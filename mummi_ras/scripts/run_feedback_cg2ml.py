#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2023, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------

import argparse
import os
import yaml
import traceback
from logging import getLogger
import numpy as np
from typing import List
import pathlib
import time

import mummi_core
import mummi_ras
from mummi_ras import Naming
from mummi_ras.feedback import FeedbackFrames


LOGGER = getLogger(__name__)

# -----------------------------------------------------------------------------
def main(mini_mummi: bool, refresh_time_sec: int = 180,  nthreads: int = None, no_encoding: bool = False):
    mummi_core.init()
    mummi_core.create_root()
    if Naming.MUMMI_MINI:
        LOGGER.warning(f"Running in mini-MuMMI mode")
        spath = os.path.join(Naming.MUMMI_SPECS, "mini-mummi")
        from mummi_ras.ml.models.autoencoders import FullMiniDenseAutoencoder as DenseAutoEncoder
    else:
        spath = os.path.join(Naming.MUMMI_SPECS, "workflow")
        # from mummi_ras.ml.models.autoencoders import FullDenseAutoencoder as DenseAutoEncoder
        # For C3
        from mummi_ras.ml.models.autoencoders import FullDenseAutoencoderUCG as DenseAutoEncoder


    config = mummi_core.read_specs(os.path.join(spath, "mlserver.yaml"))
    mummi_core.init_logger(config = config["config"])

    lockfile = os.path.join(config["config"]["lockfile"]["path"], config["config"]["lockfile"]["name"])

    # We take the first model that exists
    model_dir = None
    for path in config["models"]["ucg_autoencoder"]["paths"]:
        model_dir = path if os.path.isdir(path) else None

    if model_dir is not None:
        print(f"Error: did not find any valid paths for the UCG model")

    model_path = os.path.join(model_dir, config["models"]["ucg_autoencoder"]["model"])
    outpath = config["generator"]["outpath"]
    indice_cluster = os.path.join(config["models"]["ucg_autoencoder"]["resources"], "indices_per_cluster_dimer.npz")

    iteration_id = 0
    iteration_id_str = "iter{:02d}".format(iteration_id)

    feedbackframe_db = os.path.join(config["sampler"]["outpath"], iteration_id_str)
    feedbackframe_db = os.path.join(feedbackframe_db, "feedback")
    os.makedirs(feedbackframe_db, exist_ok=True)
    feedbackframe_db = os.path.join(feedbackframe_db, config["sampler"]["feedback"]["frame_database"])

    try:
        nthreads = nthreads or len(os.sched_getaffinity(0))
    except Exception as e:
        nthreads = 1

    model_name = config["encoder"]["model"]
    if bool(config["sampler"]["feedback"]["do_feedback"]):
        start = time.time()

        if no_encoding:
            model = None
            model_path = None
        else:
            model = DenseAutoEncoder(
                model_path = model_path,
                device = "cpu",
                indice_per_cluster = indice_cluster
            )

        # Process that watches frames being outputed by ddcmd/GROMACS
        feedback_frames = FeedbackFrames(
            database = feedbackframe_db,
            path = Naming.dir_root("feedback-cg"),
            model = model,
            model_name = model_path,
            refresh_time = refresh_time_sec,
            lockfile = lockfile,
            nthreads = nthreads
        )
        feedback_frames.start()

    else:
        LOGGER.warning(f"Feedback has been deactivated (do_feedback = False) ")

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to test the Feedback from CG to ML")
    parser.add_argument("-m", "--mini", help="Run in Mini MuMMI mode", action='store_true')
    parser.add_argument("-r", "--refresh", help="Refresh time between each pass in seconds", required=False, type=int, default=180)
    parser.add_argument("-t", "--threads", help="Number of threads to use", required=False, type=int)
    parser.add_argument("--no-encoding", help="Do not use ML model to encode frames", action='store_true')
    args = parser.parse_args()

    main(
        mini_mummi = args.mini,
        refresh_time_sec = args.refresh,
        nthreads = args.threads,
        no_encoding = args.no_encoding
    )
