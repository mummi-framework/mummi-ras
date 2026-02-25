#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2023, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------
# @authors
#           Loic Pottier <pottier1@llnl.gov>
# ------------------------------------------------------------------------------

import argparse
import csv
import sys
import os
import yaml
import mummi_ras
import traceback
from logging import getLogger
import numpy as np
from typing import List
import pathlib
import time
import multiprocessing

import mummi_core
from mummi_ras import Naming
from mummi_ras.ml.autoencoders import FullDenseAutoencoder, FullMiniDenseAutoencoder, HierarchicalCGAutoencoder
from mummi_ras.ml.validators import CGValidator

from mummi_ras.ml import ls_point as lsp
from mummi_ras.ml.samplers import ls_sampler
from mummi_ras.ml.samplers.interpolators import ot_interpolator as ot
from mummi_ras.ml.mlserver import create_sampling_db, update_sampling_db, checking_sampling_db, create_sets, write_patches
from mummi_ras.ml.samplers.interpolators.feedback_interpolator import FeedbackOTInterpolator
from mummi_ras.ml.feedback_frames import FeedbackFrames

LOGGER = getLogger(__name__)

# Note that this code is compatible with MuMMI C2

# ------------------------------------------------------------------------------
# read config file and command line arguments to get necessary settings
# ------------------------------------------------------------------------------
def read_specs(spath):
    # read the specs    
    with open(os.path.join(spath, "mlserver.yaml"), 'r') as data:
        config = yaml.load(data, Loader = yaml.FullLoader)

    if config.get("encoder") is None or config["encoder"].get("path") is None:
        raise ValueError("No encoder specified")
    if config.get("workspace") is None or config["workspace"].get("path") is None:
        raise ValueError("No workspace specified")
    if config.get("sampler") is None or config["sampler"]["feedback"].get("database") is None:
        raise ValueError("No feedback DB for sampler specified")

    # Load paths from eval strings in config
    gdict = {"mummi_ras": mummi_ras}
    for key in ["workspace", "encoder"]:
        if type(config[key]["path"]) is dict and "eval" in config[key]["path"]:
            config[key]["path"] = eval(config[key]["path"]["eval"], gdict)

    for key in ["sampler", "generator", "validator"]:
        for val in ["inpath", "outpath"]:
            if config[key].get(val) and type(config[key][val]) is dict and "eval" in config[key][val]:
                config[key][val] = eval(config[key][val]["eval"], gdict)

    return config

def start_feedback_frames(config: dict, database: str, path: str, model_name: str):
    # We have to duplicate model because multiprocessing does not support CUDA/Pytorch
    # or we should use spawn which is tricky
    model = FullDenseAutoencoder(
        model_path = model_name,
        device = "cpu"
    )

    # generator = HierarchicalCGAutoencoder(
    #     model_path = model_path
    # )

    # Process that watches frames being outputed by ddcmd
    feedback_frames = FeedbackFrames(
        database = database,
        path = path,
        model = model,
        model_name = model_name
    )
    feedback_frames.start()

# -----------------------------------------------------------------------------
def main(nsamples: int, strict: bool = False, do_feedback_frames: bool = False, mini_mummi: bool = False):
    if mini_mummi:
        config = read_specs(os.path.join(Naming.MUMMI_SPECS, "mini-mummi"))
    else:
        config = read_specs(os.path.join(Naming.MUMMI_SPECS, "workflow"))
    mummi_core.init_logger(config = config["config"])

    if mini_mummi:
        LOGGER.info("Run in mini MuMMI mode")

    model_path = os.path.join(config["encoder"]["path"], config["encoder"]["model"])
    outpath = config["generator"]["outpath"]

    training_dir = os.path.join(model_path, "training")
    training_data = list(pathlib.Path(training_dir).glob("validation_data_all.npz"))

    LOGGER.info(f"Training data {training_dir} => {len(training_data)} files")
    states = create_sets(training_data)
    # states 4 and 5 are not useful for know
    if mini_mummi:
        # in mini mummi we consider only 2 states 
        states = states[:2]
    else:
        states = states[:3]

    # Must be set by ML Server
    iteration_id = 0
    sub_sample_frac = float(config["sampler"]["sub_sample_frac"])
    assert 0.001 < sub_sample_frac <= 1

    # Needed for OTInterpolator because size of sets must be equals
    min_size = len(states[0])
    states_total_size = 0
    for s in states:
        states_total_size += len(s)
        if min_size > len(s):
            min_size = len(s)

    mask = np.random.uniform(size=(min_size,)) <= sub_sample_frac
    for i in range(len(states)):
        if len(states[i]) > min_size:
            states[i] = states[i][:min_size]
        # We sub-sample uniformly because POT is too slow for the whole states
        states[i] = list(np.array(states[i])[mask])

    new_structure_names = []
    new_ls_coords = []
    new_validation_status = []
    new_lambda = []

    iteration_id_str = "iter{:02d}".format(iteration_id)
    feedback_db_path = os.path.join(config["sampler"]["outpath"], iteration_id_str)
    feedback_db_path = os.path.join(feedback_db_path, "feedback")
    os.makedirs(feedback_db_path, exist_ok=True)
    feedback_db_path = os.path.join(feedback_db_path, config["sampler"]["feedback"]["database"])

    feedbackframe_db = os.path.join(config["sampler"]["outpath"], iteration_id_str)
    feedbackframe_db = os.path.join(feedbackframe_db, "feedback")
    os.makedirs(feedbackframe_db, exist_ok=True)
    feedbackframe_db = os.path.join(feedbackframe_db, config["sampler"]["feedback"]["frame_database"])

    if not os.path.isfile(feedback_db_path):
        create_sampling_db(feedback_db_path, model_path)

    if not checking_sampling_db(feedback_db_path, model_path):
        # raise ValueError("Wrong feedback DB")
        LOGGER.info(f"Disabling Feedback from {feedback_db_path}")
        feedback_db_path = None
        do_feedback_frames = False
    else:
        LOGGER.info(f"Feedback DB for sampling located in {feedback_db_path}")

    if do_feedback_frames:
        # Process that watches frames being outputed by simulation
        feedback_frames = multiprocessing.Process(
            target = start_feedback_frames,
            args = (config, feedbackframe_db, Naming.dir_root("feedback-cg"), model_path,)
        )
        feedback_frames.start()
        LOGGER.info(f"Started feedback DB for CG frame located in {feedbackframe_db}")
    else:
        feedbackframe_db = None

    start = time.time()
    interpolator = FeedbackOTInterpolator(
        states = states,
        kneigh = int(config["sampler"].get("kneighbors", 10)),
        lowerbound = float(config["sampler"].get("lambda_lowerbound", 0.4)),
        upperbound = float(config["sampler"].get("lambda_upperbound", 0.6)),
        num_iter_max = int(config["sampler"].get("num_iter_max", 1000000))
    )

    end = time.time() - start
    LOGGER.info(f"Interpolator created in {end:.03f} seconds for {states_total_size} LS points")
    sampler = ls_sampler.LSSampler(states, interpolator, feedback_file = feedback_db_path, feedback_frame = feedbackframe_db)

    if mini_mummi:
        generator = FullMiniDenseAutoencoder(
            model_path = model_path
        )
    else:
        generator = FullDenseAutoencoder(
            model_path = model_path
        )

    iteration_path = config["validator"]["outpath"]
    resource_name = config["validator"]["resources"]
    complex_name = config["validator"]["complex"]
    resource_path = mummi_ras.Naming.dir_res(resource_name)

    validator = CGValidator(
        iteration_path = iteration_path,
        resource_path = resource_path,
        complex_name = complex_name,
        healing = config["validator"]["healing"],
        cleanup = config["validator"]["cleanup"],
        mini_mummi = mini_mummi
    )

    total_valid_files = []
    num_new_sample = 0
    while len(total_valid_files) < nsamples:
        ls_coords = sampler.get_new_ls_points(nsamples)
        for pts in ls_coords:
            new_ls_coords.append(pts.get_coordinates())
            new_lambda.append(pts.get_lamda())
        num_new_sample += len(ls_coords)

        new_positions = generator.decode(ls_coords)

        names_array, positions_array = write_patches(
            outpath = outpath,
            iteration_id = iteration_id,
            new_positions = new_positions
        )

        new_structure_names.append(names_array)
        return_array = validator.validateArrayThreaded(iteration_id, names_array, positions_array)
        # return_array = validator.validateArray(iteration_id, names_array, positions_array)

        if len(return_array) > 0:
            new_validation_status.append(return_array[:,0])

        _, valid_files = validator.write_validation_info(
            iteration_id = iteration_id,
            return_array = return_array,
            all_structure_names = names_array
        )

        if do_feedback_frames:
            # Update sampling DB for feedback
            update_sampling_db(
                database = feedback_db_path,
                num_new_sample = len(ls_coords),
                model_name = model_path,
                structure_names = new_structure_names,
                ls_coords = new_ls_coords,
                lambda_values = new_lambda,
                validation_status = new_validation_status
            )

        LOGGER.debug(f"valid structures = {[os.path.join(validator.current_rpath, f.split('.gro')[0]) for f in valid_files]}")
        total_valid_files += list(valid_files)
        if not strict:
            break

    LOGGER.debug([os.path.join(validator.current_rpath, f.split('.gro')[0]) for f in total_valid_files])
    if do_feedback_frames:
        feedback_frames.join()
    return total_valid_files

# -----------------------------------------------------------------------------
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script to test sampler, generator and validator routines.")
    parser.add_argument("-n", "--nsample", help="Number of samples requested", required=True, type=int)
    parser.add_argument("-s", "--strict",
        help="Generate strictly the number of samples required (e.g., if set will return n valid samples)",
        action="store_true"
    )
    parser.add_argument("-f", "--feedback", help="Activate feedback generation", action="store_true")
    parser.add_argument("--mini", help="Activate Mini MuMMI", action="store_true")

    args = parser.parse_args()

    mummi_core.init()
    mummi_core.create_root()

    mini_mummi = args.mini
    if Naming.MUMMI_MINI:
        mini_mummi = True

    main(
        nsamples = args.nsample,
        strict = args.strict,
        do_feedback_frames = args.feedback,
        mini_mummi = mini_mummi
    )
