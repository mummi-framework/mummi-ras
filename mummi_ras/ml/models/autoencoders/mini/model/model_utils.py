#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2023, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------

import os
import glob
import logging
import numpy as np

LOGGER = logging.getLogger(__name__)

# ------------------------------------------------------------------------------
# simple model utilities
# ------------------------------------------------------------------------------
def load_training_history(model_path):
    dir = model_path
    files = glob.glob(os.path.join(dir, 'model_history_r*.npz'))

    LOGGER.info(f'Loading {len(files)} files from ({model_path})')

    thistory = {}
    for i, f in enumerate(files):
        with np.load(f, allow_pickle=True) as npz:
            if i == 0:
                for k in npz.files:
                    thistory[k] = [npz[k]]
            else:
                for k in npz.files:
                    thistory[k].append(npz[k])

    thistory = {k: np.array(v) for k, v in thistory.items()}
    return thistory


# ------------------------------------------------------------------------------
def save_validation_data(tpath, glb_rank=0, nvalidation_write=-1, **kwargs):

    if nvalidation_write == 0:
        return

    # all components of data should have equal number
    n = [v.shape[0] for k, v in kwargs.items()]
    assert np.unique(n).shape[0] == 1, f'Invalid validation data given, sizes = {n}'
    n = n[0]

    outfile = os.path.join(tpath, f'validation_data_r{glb_rank}.npz')

    # write all
    if nvalidation_write == -1:
        LOGGER.info(f'Saving {n} validation data to ({outfile})')
        np.savez_compressed(outfile, **kwargs)

    # write a random subset
    else:
        nvalidation_write = min(n, nvalidation_write)
        idxs = np.random.choice(n, nvalidation_write, replace=False)
        LOGGER.info(f'Saving {idxs.shape[0]} subset of validation data to ({outfile})')
        np.savez_compressed(outfile, **{k: v[idxs] for k, v in kwargs.items()})


# ------------------------------------------------------------------------------
def load_validation_data(tpath, glb_rank=-1, keys_to_load = []):

    # need to function to decide whether to read a key or not
    read_key = lambda k: k in keys_to_load if len(keys_to_load) > 0 else True

    # first, figure out of all files need to be loaded or just for one rank
    if glb_rank != -1:
        fnames = [os.path.join(tpath, f'validation_data_r{glb_rank}.npz')]
    else:
        fnames = sorted(glob.glob(os.path.join(tpath, 'validation_data_r*.npz')))

    # read all files
    LOGGER.info(f'Loading {len(fnames)} files from ({tpath})')
    vdata = {k: [] for k in keys_to_load}
    for f in fnames:
        with np.load(f, allow_pickle=True) as npz:
            for k in np.intersect1d(npz.files, keys_to_load):
                vdata[k].append(npz[k])

    vdata = {k: np.concatenate(v) for k,v in vdata.items()}
    return vdata

# ------------------------------------------------------------------------------