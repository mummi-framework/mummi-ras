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
import os
import yaml
from logging import getLogger
import numpy as np
from typing import List
import time
import pickle
import pathlib

import mummi_core
import mummi_ras
from mummi_ras import Naming

from mummi_ras.ml.mlserver import create_sets
from mummi_ras.ml.samplers.interpolators.feedback_interpolator import FeedbackOTInterpolator

LOGGER = getLogger(__name__)

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

# -----------------------------------------------------------------------------
def main(fraction: float, lambda_lower: float, lambda_upper: float):
    mummi_core.init()
    mummi_core.create_root()
    if Naming.MUMMI_MINI:
        config = read_specs(os.path.join(Naming.MUMMI_SPECS, "mini-mummi"))
    else:
        config = read_specs(os.path.join(Naming.MUMMI_SPECS, "workflow"))
    mummi_core.init_logger(config = config["config"])

    model_path = os.path.join(config["encoder"]["path"], config["encoder"]["model"])
    outpath = config["generator"]["outpath"]

    training_dir = os.path.join(model_path, "training")
    training_data = list(pathlib.Path(training_dir).glob("validation_data_all.npz"))

    LOGGER.info(f"Training data {training_dir} => {len(training_data)} files")
    states = create_sets(training_data)

    # states 4 and 5 are not useful for know
    states = states[:3]

    sub_sample_frac = fraction
    assert 0.001 < sub_sample_frac <= 1

    assert 0 <= lambda_lower <= 1
    assert 0 <= lambda_upper <= 1
    assert lambda_lower < lambda_upper

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

    LOGGER.info(f"sub_sample_frac   = {fraction}")
    LOGGER.info(f"lambda            = [{lambda_lower}, {lambda_upper}]")
    LOGGER.info(f"kneigh            = 10")
    LOGGER.info(f"num_iter_max      = 1000000")

    start = time.time()
    interpolator = FeedbackOTInterpolator(
        states = states,
        kneigh = 10,
        lowerbound = lambda_lower,
        upperbound = lambda_upper,
        num_iter_max = 1000000
    )

    end = time.time() - start
    LOGGER.info(f"Interpolator created in {end:.03f} seconds for {states_total_size} LS points")
    
    output = f"interpolator_frac-{sub_sample_frac}_ll-{lambda_lower}_lu-{lambda_upper}.pkl"

    start = time.time()
    with open(output, 'wb') as f:
        pickle.dump(interpolator, f)

    end = time.time() - start
    LOGGER.info(f"Pickle file ({output}) of Interpolator created in {end:.03f} seconds (size = {len(interpolator.all_points)})")
    
    # Test if we can re-open the pickled file
    
    start = time.time()
    with open(output, 'rb') as f:
        interpolator_pkl = pickle.load(f)

    end = time.time() - start

    LOGGER.info(f"Loading pickle file ({output}) done in {end:.03f} seconds (size = {len(interpolator_pkl.all_points)})")

# -----------------------------------------------------------------------------
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script to generate a pickled interpolator")
    parser.add_argument("-f", "--fraction", help="Fraction of states to use", required=True, type=float)
    parser.add_argument("-lu", "--lambda-upper", help="Upper bound for lambda", required=False, type=float, default=1.0)
    parser.add_argument("-ll", "--lambda-lower", help="Lower bound for lambda", required=False, type=float, default=0.0)

    args = parser.parse_args()

    main(
        fraction = args.fraction,
        lambda_lower = args.lambda_lower,
        lambda_upper = args.lambda_upper
    )
