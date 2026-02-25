#!/usr/bin/env python3
#
# Copyright (c) 2025, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------

from .autoencoders import FullDenseAutoencoder, FullMiniDenseAutoencoder, HierarchicalCGAutoencoder
from .generatives import V_Generator, V2UCG_Generator, UCG2CG_Generator

import glob
import logging
import os

import numpy as np
from typing import Union, List, Any

LOGGER = logging.getLogger(__name__)

# Utility functions

def latest_patch_id_generated(RPATH: str) -> int:
    """
    Find in the directory containing potentially already some structures
    like structure_iter00_X.gro. This function finds the highest X so we can
    write the patch X+1 safely without erasing any structures.

    This function returns X+1
    """
    nsamples = len(glob.glob(RPATH + '/structure*.npz'))
    # nsamples = nsamples - 3 # the - 3 is to remove master_iterX.npz and table_all_predictions.npz and because the structure id starts from 0
    return max(nsamples, 0) # to avoid returning a negative number if RPATH was empty

def write_patches(outpath: str, iteration_id: int, new_positions: np.array) -> Union[List[str], List[Any]]:
    """
    From positions, write a patch as a npz file and 
    return its path.
    """
    iteration_id = "{:02d}".format(iteration_id)
    RPATH = os.path.join(outpath, f"iter{iteration_id}")
    if not os.path.isdir(RPATH):
        os.makedirs(RPATH, exist_ok=True)

    LOGGER.info(f"Write patches in {RPATH}")
    npatches = new_positions.shape[0]
    files = np.empty(npatches, dtype='object')
    positionsArray = []
    namesArray = []
    all_structure_names = []

    # We want to avoid erasing previously generated structures
    offset = latest_patch_id_generated(RPATH)

    # Create an .npz for each new patch
    for sidx in range(npatches):

        positions = new_positions[sidx]
        structure_id = "{:012d}".format(offset+sidx)

        structure_name = f"structure_iter{iteration_id}_{structure_id}"
        all_structure_names.append(structure_name)
        outfile_positions = os.path.join(RPATH, f"{structure_name}.npz")
        files[sidx] = f"{structure_name}.npz"
        LOGGER.info(f"Writing {outfile_positions}")
        np.savez_compressed(outfile_positions, data = positions)

        positionsArray.append(positions)
        namesArray.append(structure_name)

    return namesArray, positionsArray
