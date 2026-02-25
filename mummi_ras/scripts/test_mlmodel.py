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
from logging import getLogger
import torch
import yaml

import mummi_core
from mummi_ras.ml.models.model_ae_inf import CustomAutoencoder as CGAutoencoder

LOGGER = getLogger(__name__)
mummi_core.init()
mummi_core.init_logger(level=2)

FRAMEWORK = "torch"

class TestModel():
    def __init__(self, model_path: str, outpath: str) -> None:
        self.model_path = model_path
        self.config_fname = os.path.join(self.model_path, "config.yaml")
        with open(self.config_fname, 'r') as data:
            self.config = yaml.load(data, Loader = yaml.Loader)
        self.config["model"]["name"] = os.path.basename(self.model_path)
        
        self.outpath = outpath
        self.inshape = tuple(self.config["data"]["inshape"])
        self.outshape = tuple(self.config["model"]["z_dim"])
        self.config["model"]["type"] = "protein_latent_space"
        self.config["model"]["inshape"] = self.inshape
        self.config["outpath"] = self.outpath

        LOGGER.info(f"> Initialization of generator with:")
        LOGGER.info(f"  > model_path = {model_path}")
        LOGGER.info(f"  > outpath    = {outpath}")
        LOGGER.info(f"  > encoder    = {self.model_path}")
        LOGGER.info(f"  > dimensions = {self.outshape}")

        self.model = CGAutoencoder(
            inshape = self.inshape,
            latent_dims = self.outshape, 
            INPATH = self.model_path,
            config = self.config,
            tstamp = '0',
            create_subdirs = False
        )

        device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
        # device = torch.device('cpu')
        LOGGER.info(f" > Device used : {device}")
        if torch.cuda.is_available():
            LOGGER.info(f"   > CUDA total devices         : {torch.cuda.device_count()}")
            LOGGER.info(f"   > CUDA active device         : device:{torch.cuda.current_device()}")
            # if torch.backends.cudnn.is_available() and torch.backends.cudnn.enabled:
            #     LOGGER.info(f"   > cudnn version               : {torch.backends.cudnn.version()}")

        self.model = self.model.to(device)
        self.model.device = device   # so the model can move the data to device later
        model_name = os.path.join(self.model_path, "model.pth")
        std_fname = os.path.join(self.model_path, "CG_data_summary_pos_dis.npz")
        self.model.load_model(model_name, std_fname, map_location = device)

def main():
    test = TestModel(
        model_path="PATH_OF_THE_MODEL",
        outpath = "test"
        )

if __name__ == "__main__":
    main()
