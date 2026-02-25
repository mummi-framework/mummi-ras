#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------

import argparse
import os
import logging
import numpy as np
from typing import List, Any
from timeit import default_timer as timer

import mummi_core
from mummi_ras import Naming
from mummi_ras.ml.validators import CGValidator

LOGGER = logging.getLogger(__name__)

# ------------------------------------------------------------------------------
# setup command line arguments
# ------------------------------------------------------------------------------
def setup_argparser():
    parser = argparse.ArgumentParser(description='')
    # inpath and outpath are not in config files
    parser.add_argument('--positions', '-p', required=True,
                        type=str, dest='positions', help='RAS-RBD-CRD 3d positions')
    parser.add_argument('--iteration-path', '-o', required=False, default = os.path.abspath('.'),
                        type=str, dest='iteration', help='Output folder default (.)')
    parser.add_argument("-d", "--loglevel", type=int, default=1,
                        help="Level of logging messages to be output:\n"
                             "5 - Critical\n"
                             "4 - Error\n"
                             "3 - Warning\n"
                             "2 - Info (Default)\n"
                             "1 - Debug")
    parser.add_argument("-l", "--logpath", type=str, default="./",
                        help="Path for log information")
    parser.add_argument("-c", "--logstdout", action="store_true",
                        help="Log to stdout in addition to a file")
    parser.add_argument("--logfile", type=str, default="validator_log",
                        help="Name of validator logger")

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    mummi_core.init()

    # ge
    args = setup_argparser()
    mummi_core.init_logger(argparser=args)
    # --------------------------------------------------------------------------
    # test structure for specified file
    # --------------------------------------------------------------------------
    positions_files = args.positions.split()

    # The validator wants a directory named iterXX
    iteration_id = 0
    iteration_str = "iter{:02d}".format(iteration_id)
    os.makedirs(iteration_str, exist_ok=True)

    # Read 3d positions
    positionsArray = []
    namesArray = []
    for pos_file in positions_files:
        positionsArray.append(np.load(pos_file, allow_pickle=True)['data'])
        namesArray.append(pos_file.split(".")[0])

    resource_path = Naming.dir_res("martini3-validator")
    complex_name = "ras-raf-14-3-3-ref-CG.gro"
    print(f"Run validator for RAS-RBD-CRD structures {namesArray} coordinates in {positions_files}")
    start = timer()
    validator = CGValidator(
        iteration_path = os.path.abspath(args.iteration),
        resource_path = resource_path,
        complex_name = complex_name,
        cleanup = False
    )

    #valid_status, gromacs_files = validator.validate(name, positions)
    returnArray = validator.validateArray(iteration_id, namesArray, positionsArray)
    # returnArray = validator.validateArrayThreaded(namesArray, positionsArray, 6)
    print(f"Validator returned {returnArray} \n len {len(returnArray)} took {timer() - start} s")
