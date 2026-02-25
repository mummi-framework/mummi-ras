#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------

import os
import yaml
import mummi_ras
from logging import getLogger
import numpy as np
import pathlib
import argparse
from timeit import default_timer as timer
import csv
import time

# import mummi_core
# from mummi_ras import Naming
from mummi_ras.ml.validators import CGValidator
from mummi_ras.ucg.validators import UCGValidator

LOGGER = getLogger(__name__)

# -----------------------------------------------------------------------------
def main(positions_dir: str, validator_resources: str, outpath: str):
    if os.path.isdir(outpath) != True:
        os.mkdir(outpath)

    iteration_id = 0 # We do not care about the iteration here but we need to create the folder
    iteration_id_str = "{:02d}".format(iteration_id)
    RPATH = os.path.join(outpath, f"iter{iteration_id_str}")
    if os.path.isdir(RPATH) != True:
        os.mkdir(RPATH)

    position_files = list(pathlib.Path(positions_dir+"/").glob("structure_iter*.npz"))
    if position_files == []:
        print(f"No structure files found in {positions_dir}. Abort.")
        return []
    else:
        print(f"Structure files found in {positions_dir} => {len(position_files)}")
        print(f"Output directory => {outpath}")

    validator = CGValidator(iteration_path = outpath, resource_path = validator_resources, cleanup = True)
    # validator = UCGValidator(iteration_path = outpath, resource_path = validator_resources, cleanup = True)

    names_array = []
    positions_array = []
    for file in position_files:
        names_array.append(str(file.stem))
        positions_array.append(np.load(file, allow_pickle=True)['data'])

    return_array = validator.validateArray(iteration_id, names_array, positions_array)

    # return_array = np.array(return_array)
    # _, valid_files = validator.write_validation_info(
    #     iteration_id = iteration_id,
    #     return_array = return_array,
    #     all_structure_names = names_array
    # )

    return np.array(return_array)

# -----------------------------------------------------------------------------
if __name__ == "__main__":
   # get arguments and configuration
    parser = argparse.ArgumentParser(description='Script to validate structures with a given validator (CG, UCG etc)')
    parser.add_argument('--patch-dir', '-p', required=True, type=str, dest='positions', help='Path containing structure files as npz.')
    parser.add_argument('--validator-resources', '-r', required=True, type=str, dest='validator', help='Path containing validator resources.')
    parser.add_argument('--outpath', '-o', required=True, type=str, dest='outpath', help='Path where to output files (preferably absolute paths).')
    args = parser.parse_args()

    # Example:
    #   python3 run_validation.py -o /p/vast1/mummiusr/outputs -p $MUMMI_ROOT/ml/iterXX -r $MUMMI_RESOURCES/martini3-validator-resources
    start = timer()

    results = main(
        positions_dir = args.positions,
        validator_resources = args.validator,
        outpath = os.path.abspath(args.outpath)
    )
    # position_files = list(pathlib.Path(args.positions+"/").glob("structure_iter*.npz"))
    valid_files = np.where(results == "True")[0]

    timestr = time.strftime("%Y%m%d-%H%M%S")
    output_csv = f"output-{timestr}.csv"
    with open(output_csv, mode='w+') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(["ID", "Valid", "MF", "PE"])
        for row in results:
            writer.writerow([row[1].split(".")[0], row[0], row[2], row[3]])

    nbvalids_str = "{:04d}".format(valid_files.shape[0])
    total_str = "{:04d}".format(results.shape[0])
    total_time = "{:.3f}".format(timer()-start)

    print(f"Valid patches found {nbvalids_str} / {total_str} in {total_time} seconds")
    print(f"Results written in => {output_csv}")