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
from typing import List
import time
import traceback

import mummi_core
import mummi_ras
from mummi_ras import Naming
from mummi_ras.feedback import FeedbackCreatesims

LOGGER = getLogger(__name__)

# ------------------------------------------------------------------------------
# read config file and command line arguments to get necessary settings
# ------------------------------------------------------------------------------
def read_specs(spath : str) -> dict:
    # read the specs
    try:  
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

        with open(os.path.join(spath, 'wfmanager.yaml'), 'r') as data:
            wfmngr = yaml.load(data, Loader=yaml.FullLoader)

        try:
            iteration_id = int(wfmngr["wfmanager"]["config"].get("mlserver_round_id"))
        except Exception as _:
            iteration_id = 0
        config["iteration_id"] = iteration_id
        iteration_path = "iter{:02d}".format(iteration_id)
        config["iteration_path"] = iteration_path

        do_feedback = bool(config["sampler"]["feedback"]["do_feedback"])

        sampling_db = os.path.join(config["sampler"]["outpath"], iteration_path)
        sampling_db = os.path.join(sampling_db, "feedback")

        os.makedirs(sampling_db, exist_ok = True)

        sampling_db = os.path.join(sampling_db, config["sampler"]["feedback"]["database"])
        
        model_name = config["encoder"]["model"]
        config["model_name"] = model_name
        config["do_feedback"] = do_feedback
        config["sampling_db"] = sampling_db

        if "mini_mummi" in config["config"]:
            mini_mummi = bool(config["config"]["mini_mummi"])

        if mini_mummi:
            LOGGER.info("Run in mini MuMMI mode")

        LOGGER.info(f"> Initializing Createsim Feedback")
    except Exception as e:
        traceback.print_exc()
        LOGGER.error("Createsim Feedback failed during reading specs!")
        raise e

    return config

# -----------------------------------------------------------------------------
def main(mini_mummi: bool, refresh_time_sec: int = 180,  nthreads: int = None):
    mummi_core.init()
    mummi_core.create_root()
    if mini_mummi or Naming.MUMMI_MINI:
        LOGGER.warning(f"Running in mini-MuMMI mode")
        config = read_specs(os.path.join(Naming.MUMMI_SPECS, "mini-mummi"))
    else:
        config = read_specs(os.path.join(Naming.MUMMI_SPECS, "workflow"))

    config["config"]["loglevel"] = 1 # Forcing debug level for tests
    mummi_core.init_logger(config = config["config"])

    try:
        nthreads = nthreads or len(os.sched_getaffinity(0))
    except Exception as e:
        nthreads = 1

    feedback = FeedbackCreatesims(
        dir_to_watch =  Naming.dir_root("all-cg"),
        sampling_db = config["sampling_db"],
        model_name = config["model_name"],
        waittime = refresh_time_sec, # In seconds
        nthreads = nthreads,
        logger = LOGGER
    )
    if bool(config["do_feedback"]):
        feedback.launch()

# -----------------------------------------------------------------------------
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script to test the createsims LS sampling DB feature")
    parser.add_argument("-m", "--mini", help="Run in Mini MuMMI mode", action='store_true')
    parser.add_argument("-r", "--refresh", help="Refresh time between each pass in seconds", required=False, type=int, default=180)
    parser.add_argument("-t", "--threads", help="Number of threads to use", required=False, type=int)
    args = parser.parse_args()

    main(
        mini_mummi = args.mini,
        refresh_time_sec = args.refresh,
        nthreads = args.threads
    )