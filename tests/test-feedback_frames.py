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
def main(mini_mummi: bool, refresh_time_sec: int = 180,  nthreads: int = None):
    mummi_core.init()
    mummi_core.create_root()
    if mini_mummi or Naming.MUMMI_MINI:
        LOGGER.warning(f"Running in mini-MuMMI mode")
        from mummi_ras.ml.autoencoders import FullMiniDenseAutoencoder as DenseAutoEncoder
        config = read_specs(os.path.join(Naming.MUMMI_SPECS, "mini-mummi"))
    else:
        from mummi_ras.ml.autoencoders import FullDenseAutoencoder as DenseAutoEncoder
        config = read_specs(os.path.join(Naming.MUMMI_SPECS, "workflow"))

    config["config"]["loglevel"] = 1 # Forcing debug level for tests

    mummi_core.init_logger(config = config["config"])

    # We take the first model that exists
    if isinstance(config["encoder"]["path"], list):
        for epath in config["encoder"]["path"]:
            model_path = os.path.join(epath, config["encoder"]["model"])
            if os.path.exists(model_path):
                break
    else:
        model_path = os.path.join(config["encoder"]["path"], config["encoder"]["model"])
    outpath = config["generator"]["outpath"]

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

        model = DenseAutoEncoder(
            model_path = model_path,
            device = "cpu"
        )

        # Process that watches frames being outputed by ddcmd/GROMACS
        feedback_frames = FeedbackFrames(
            database = feedbackframe_db,
            path = Naming.dir_root("feedback-cg"),
            model = model,
            model_name = model_path,
            refresh_time = refresh_time_sec,
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
    args = parser.parse_args()

    main(
        mini_mummi = args.mini,
        refresh_time_sec = args.refresh,
        nthreads = args.threads
    )
