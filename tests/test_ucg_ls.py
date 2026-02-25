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
import csv
import sys
import os
import yaml
import mummi_ras
import traceback
from logging import getLogger
import numpy as np
from typing import List
import pathlib
import time
import multiprocessing

import mummi_core
from mummi_core import read_specs
from mummi_ras import Naming
from mummi_ras.ml.models.autoencoders import FullDenseAutoencoderUCG

# from mummi_ras.ml.validators import CGValidator

# from mummi_ras.ml import ls_point as lsp
# from mummi_ras.ml.samplers import ls_sampler
# from mummi_ras.ml.samplers.interpolators import ot_interpolator as ot
# from mummi_ras.ml.mlserver import create_sampling_db, update_sampling_db, checking_sampling_db, create_sets, write_patches
# from mummi_ras.ml.samplers.interpolators.feedback_interpolator import FeedbackOTInterpolator
# from mummi_ras.ml.feedback_frames import FeedbackFrames


LOGGER = getLogger(__name__)

# # ------------------------------------------------------------------------------
# # read config file and command line arguments to get necessary settings
# # ------------------------------------------------------------------------------
# def read_specs(spath):
#     # read the specs    
#     with open(os.path.join(spath, "mlserver.yaml"), 'r') as data:
#         config = yaml.load(data, Loader = yaml.FullLoader)

#     if config.get("encoder") is None or config["encoder"].get("path") is None:
#         raise ValueError("No encoder specified")
#     if config.get("workspace") is None or config["workspace"].get("path") is None:
#         raise ValueError("No workspace specified")
#     if config.get("sampler") is None or config["sampler"]["feedback"].get("database") is None:
#         raise ValueError("No feedback DB for sampler specified")

#     # Load paths from eval strings in config
#     gdict = {"mummi_ras": mummi_ras}
#     for key in ["workspace", "encoder"]:
#         if type(config[key]["path"]) is dict and "eval" in config[key]["path"]:
#             config[key]["path"] = eval(config[key]["path"]["eval"], gdict)

#     for key in ["sampler", "generator", "validator"]:
#         for val in ["inpath", "outpath"]:
#             if config[key].get(val) and type(config[key][val]) is dict and "eval" in config[key][val]:
#                 config[key][val] = eval(config[key][val]["eval"], gdict)

#     return config

# def start_feedback_frames(config: dict, database: str, path: str, model_name: str):
#     # We have to duplicate model because multiprocessing does not support CUDA/Pytorch
#     # or we should use spawn which is tricky
#     model = FullDenseAutoencoderUCG(
#         model_path = model_name,
#         device = "cpu"
#     )

#     # generator = HierarchicalCGAutoencoder(
#     #     model_path = model_path
#     # )

#     # Process that watches frames being outputed by ddcmd
#     feedback_frames = FeedbackFrames(
#         database = database,
#         path = path,
#         model = model,
#         model_name = model_name
#     )
#     feedback_frames.start()

# -----------------------------------------------------------------------------
def main(input_feedback: List[str], mini_mummi: bool = False):
    mummi_core.init()
    mummi_core.create_root()

    if mini_mummi:
        spath = os.path.join(mummi_ras.Naming.MUMMI_SPECS, "mini-mummi")
    else:
        spath = os.path.join(mummi_ras.Naming.MUMMI_SPECS, "workflow")

    config = mummi_core.read_specs(os.path.join(spath, "mlserver.yaml"))

    mummi_core.init_logger(config = config["config"])

    if mini_mummi:
        LOGGER.info("Run in mini MuMMI mode")

    model_dir = None
    for path in config["models"]["ucg_autoencoder"]["paths"]:
        model_dir = path if os.path.isdir(path) else None

    if model_dir:
        print(f"Error: did not find any valid paths for the UCG model")

    model_path = os.path.join(model_dir, config["models"]["ucg_autoencoder"]["model"])
    outpath = config["generator"]["outpath"]

    indice_cluster = os.path.join(config["models"]["ucg_autoencoder"]["resources"], "indices_per_cluster_dimer.npz")

    generator = FullDenseAutoencoderUCG(
        model_path = model_path,
        indice_per_cluster = indice_cluster
    )

    ls_coords = []
    for positions_file in input_feedback:
        print(f"Positions file {positions_file}")
        d = np.load(positions_file, allow_pickle=True)
        CG_positions = d["positions"]
        print(f"CG_positions {CG_positions.shape}")

        ls_coords.append(generator.encode(CG_positions))
    ls_coords = np.array(ls_coords)
    print(ls_coords.shape)

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to test UCG LS for feedback sampling.")
    parser.add_argument("-i", "--feedback", help="Directory where the NPZ are located", nargs='+', required=True)

    args = parser.parse_args()

    main(input_feedback = args.feedback, mini_mummi = False)
