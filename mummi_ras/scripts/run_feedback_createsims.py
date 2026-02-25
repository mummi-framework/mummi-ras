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

from mummi_ras.ml.models.autoencoders import FullDenseAutoencoderUCG as DenseAutoEncoder

LOGGER = getLogger(__name__)

def select_mlmodels(lists_paths):
    """
    The YAML config can provided multiple paths for the same model.
    THis happens when we run on different machines. Models are located on
    different paths. We select the first valid one.
    """
    if not isinstance(lists_paths, list):
        lists_paths = [lists_paths]
    if isinstance(lists_paths, list):
        for epath in lists_paths:
            if os.path.exists(epath):
                LOGGER.info(f"Found valid path {epath}")
                return epath

# -----------------------------------------------------------------------------
def main(mini_mummi: bool, refresh_time_sec: int = 180,  nthreads: int = None, run_campaign2: bool = False):
    mummi_core.init()
    mummi_core.create_root()
    if Naming.MUMMI_MINI:
        LOGGER.warning(f"Running in mini-MuMMI mode")
        spath = os.path.join(Naming.MUMMI_SPECS, "mini-mummi")
    else:
        spath = os.path.join(Naming.MUMMI_SPECS, "workflow")

    config = mummi_core.read_specs(os.path.join(spath, "mlserver.yaml"))
    mummi_core.init_logger(config = config["config"])

    lockfile = os.path.join(config["config"]["lockfile"]["path"], config["config"]["lockfile"]["name"])

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
    config["do_feedback"] = do_feedback
    config["sampling_db"] = sampling_db

    if run_campaign2:
        config["model_name"] = model_name
    else:
        checkpoint_v2ucg_paths = select_mlmodels(config["models"]["generative"]["v2ucg"]["paths"])
        checkpoint_ucg2cg_paths = select_mlmodels(config["models"]["generative"]["ucg2cg"]["paths"])
        config["model_name"] = [checkpoint_v2ucg_paths, checkpoint_ucg2cg_paths]

    # We take the first model that exists
    ucg_model_dir = None
    ucg_model_path = None
    for path in config["models"]["ucg_autoencoder"]["paths"]:
        ucg_model_dir = path if os.path.isdir(path) else None

    if ucg_model_dir is None:
        print(f"Error: did not find any valid paths for the UCG model")
    else:
        ucg_model_path = os.path.join(ucg_model_dir, config["models"]["ucg_autoencoder"]["model"])
        indice_cluster = os.path.join(config["models"]["ucg_autoencoder"]["resources"], "indices_per_cluster_dimer.npz")
        ucg_model = DenseAutoEncoder(
            model_path = ucg_model_path,
            device = "cpu",
            indice_per_cluster = indice_cluster
        )

    try:
        nthreads = nthreads or len(os.sched_getaffinity(0))
    except Exception as e:
        nthreads = 1
    if bool(config["do_feedback"]):
        feedback = FeedbackCreatesims(
            dir_to_watch =  Naming.dir_root("all-cg"),
            sampling_db = config["sampling_db"],
            model_name = config["model_name"],
            ls_encoder = ucg_model,
            waittime = refresh_time_sec, # In seconds
            nthreads = nthreads,
            lockfile = lockfile,
            logger = LOGGER
        )
        feedback.launch()
    else:
        LOGGER.warning(f"Feedback is not activated in specs/")
        return

# -----------------------------------------------------------------------------
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script to test the createsims LS sampling DB feature")
    parser.add_argument("-m", "--mini", help="Run in Mini MuMMI mode", action='store_true')
    parser.add_argument("-r", "--refresh", help="Refresh time between each pass in seconds", required=False, type=int, default=180)
    parser.add_argument("-t", "--threads", help="Number of threads to use", required=False, type=int)
    parser.add_argument("--campaign-2", help="Run campaign 2 model (by default run Campaign 3 model)", action='store_true')

    args = parser.parse_args()

    main(
        mini_mummi = args.mini,
        refresh_time_sec = args.refresh,
        nthreads = args.threads,
        run_campaign2 = args.campaign_2
    )