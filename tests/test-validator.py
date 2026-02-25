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
import sys
import os
import yaml
import mummi_ras
from logging import getLogger
from typing import List
import time
import pathlib

import torch
import numpy as np

import mummi_core
from mummi_ras import Naming

from mummi_ras.ml.models import write_patches
from torch.utils.data import DataLoader
from mummi_ras.ml.validators import CGValidator_C3 as CGValidator

# from mummi_ras.ml.mlserver_c3 import create_sampling_db, update_sampling_db, checking_sampling_db


LOGGER = getLogger(__name__)

# -----------------------------------------------------------------------------
def main(input_samples: str, output: str, nthreads: int):
    spath = os.path.join(Naming.MUMMI_SPECS, "workflow")
    config = mummi_core.read_specs(os.path.join(spath, "mlserver.yaml"))
    mummi_core.init_logger(config = config["config"])

    os.makedirs(output, exist_ok=True)
    LOGGER.info(f"Created output directory {output}")

    iteration_id = 0 # Shortcut here: Must be read from wfmanager.yaml
    iteration_id_str = "iter{:02d}".format(iteration_id)

    # sampling_db = os.path.join(config["sampler"]["outpath"], iteration_id_str)
    # sampling_db = os.path.join(sampling_db, "feedback")
    # os.makedirs(sampling_db, exist_ok=True)
    # sampling_db = os.path.join(sampling_db, config["sampler"]["feedback"]["database"])

    # LOGGER.info(f"Feedback DB {sampling_db}")

    iteration_path = output
    resource_name = config["validator"]["resources"]
    complex_name = config["validator"]["complex"]
    resource_path = mummi_ras.Naming.dir_res(resource_name)

    validator = CGValidator(
        iteration_path = iteration_path,
        resource_path = resource_path,
        complex_name = complex_name,
        healing = config["validator"]["healing"],
        cleanup = config["validator"]["cleanup"],
        nthreads = nthreads,
        mini_mummi = False
    )

    new_structure_names = []
    new_validation_status = []
    # RPATH = os.path.join(iteration_path, f"{iteration_id_str}")
    # name_table_all_generations = os.path.join(os.path.join(RPATH, "table_all_generations.npz"))

    LOGGER.info(f"Reading {input_samples}")

    gen_data = np.load(input_samples, allow_pickle=True)
    pred_cg = gen_data["pred_cg"]
    LOGGER.info(f"pred_cg = {pred_cg.shape}")

    names_array, _ = write_patches(
        outpath = iteration_path,
        iteration_id = iteration_id,
        new_positions = pred_cg
    )

    # Validation
    validation_time = time.time()

    new_structure_names.append(names_array)
    ## Can return more than one structure
    return_array = validator.validateArrayThreaded(iteration_id, names_array, pred_cg)
    if len(return_array) > 0:
        new_validation_status.append(validator.get_status(return_array))

    LOGGER.info(return_array)
    LOGGER.info(return_array.shape)

    _, valid_files = validator.write_validation_info(
        iteration_id = iteration_id,
        return_array = return_array,
        all_structure_names = names_array
    )

    validation_time = time.time() - validation_time
    LOGGER.info(f"Validation: sec/sample: {validation_time/pred_cg.shape[0]} with {nthreads} threads")

    # # Save the details for all the predictions
    status = np.array([ff for struct in return_array[:,0] for ff in struct])
    gromacs_files = np.array([ff for struct in return_array[:,1] for ff in struct])
    valid_ids = np.where(status)[0]
    valid_files = gromacs_files[valid_ids]
    status_uniq_str = []
    for f in names_array:
        if np.any(np.char.find(valid_files, f) >= 0):
            status_uniq_str.append(True)
        else:
            status_uniq_str.append(False)

    arrays = {key: gen_data[key] for key in gen_data.files}
    arrays["structure_names"] = names_array
    arrays["status"] = status_uniq_str
    arrays["selected"] = np.full(len(status_uniq_str), False)

    # File we will be writing (just in case to not corrupt the file in case of interruption)
    timestr = time.strftime("%Y%m%d-%H%M%S")
    db_path = os.path.splitext(input_samples)
    input_samples_tmp = f"{db_path[0]}-{timestr}{db_path[1]}"

    np.savez_compressed(input_samples_tmp, **arrays)

    # Now we replace the files (atomic)
    os.replace(input_samples_tmp, input_samples)

    # update_sampling_db(
    #     database = sampling_db,
    #     model_name = [checkpoint_v2ucg_paths, checkpoint_ucg2cg_paths],
    #     num_new_sample = len(new_v_coords),
    #     structure_names = new_structure_names,
    #     v_coords = new_v_coords,
    #     v2ucg_coords = new_v2ucg_coords,
    #     ucg2cg_coords = new_ucg2cg_coords,
    #     validation_status = new_validation_status
    # )

    return return_array

# -----------------------------------------------------------------------------
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script to test validator routine.")
    parser.add_argument("-i", "--input-npz", help="NPZ with structures", required=True, type=str)
    parser.add_argument("-o", "--output", help="Where to write all the structures (usually $MUMMI_ROOT/ml/iter00)", required=True, type=str)
    parser.add_argument("-t", "--threads", help="Number of threads for the validation", required=False, default=100)
    args = parser.parse_args()

    mummi_core.init()
    mummi_core.create_root()

    d = main(
        input_samples = args.input_npz,
        output = args.output,
        nthreads = args.threads
    )
    LOGGER.info(f"Returned {d.shape}")
    
