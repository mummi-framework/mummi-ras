#!/usr/bin/env python3
#
# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------


import torch
import numpy as np
import os
import yaml
from logging import getLogger, Logger

from torch import Tensor
import numpy.typing as npt

from mummi_ras.ml import ls_point as lsp
from ..base.cgautoencoder import CGAutoencoder
# from .pc3 import CustomAutoencoder as CGAutoencoder_model
# from .c1 import CustomAutoencoder as CGAutoencoder_model
from .c2 import CustomAutoencoder as CGAutoencoder_model



# ------------------------------------------------------------------------------
LOGGER = getLogger(__name__)
FRAMEWORK = "torch"
# ------------------------------------------------------------------------------

class FullDenseAutoencoder(CGAutoencoder):
    """Class of a full dense CG autoencoder 
    
    Attributes (public):
        model_path (str): Path containing trained ML model (parameter or checkpoint files).
        encoder (PyTorch module): Encoder module.
        decoder (PyTorch module): Decoder module.
        
    Methods (public):
        encode(pos): Encodes CG positions into latent space coordinates.
        decode(l): Decodes latent space coordinates into CG positions.
    """
    
    def __init__(self, model_path: str, device: str = "cuda") -> None:
        super().__init__(model_path)
        self.framework = FRAMEWORK
        is_framework_supported = self._setup()
        assert is_framework_supported, f"Invalid framework: {self.framework}"
        self.config_fname = os.path.join(self.model_path, "config.yaml")
        with open(self.config_fname, 'r') as data:
            self.config = yaml.load(data, Loader = yaml.Loader)
        self.config["model"]["name"] = os.path.basename(self.model_path)

        self.inshape = tuple(self.config["data"]["inshape"])
        self.outshape = tuple(self.config["model"]["z_dim"])
        self.config["model"]["type"] = "protein_latent_space"
        self.config["model"]["inshape"] = self.inshape

        LOGGER.info(f"> Initialization of generator with:")
        LOGGER.info(f"  > encoder      = {self.model_path}")
        LOGGER.info(f"  > dimensions   = ({self.inshape}, {self.outshape})")

        self._load_model(device)

    def encode(self, pos: npt.NDArray) -> npt.NDArray:
        """Encodes CG positions into latent space.

        Args:
            pos (npt.NDArray): CG positions with shape [B, N, 3], where B is number of batchs,
                and N is number of CG beads.

        Returns:
            l (npt.NDArray): Latent space coordinates with shape [B, D], where B is number of batchs,
                and D is latent space dimensionality.
        """
        std_pos = self.model.standardize(pos)
        std_pos = np.reshape(std_pos, (std_pos.shape[0] * std_pos.shape[1])) # Flatten for the input if the ML model
        std_pos = np.expand_dims(std_pos, axis=0) # to have a (1, 8691)
        return self.model.encode_ae(std_pos)

    def decode(self, l: npt.NDArray[lsp]) -> npt.NDArray:
        """Decodes latent space coordinates into CG positions.

        Args:
            l (npt.NDArray[lsp]): Latent space coordinates with shape [B, D], where B is number of batchs,
                and D is latent space dimensionality.

        Returns:
            pos (v): CG positions with shape [B, N, 3], where B is number of batchs,
                and N is number of CG beads.
        """
        res = np.array([pts.get_coordinates() for pts in l])
        return self._decode(ls_coords = res)

    def _setup(self) -> bool:
        """Checks the used framework and returns True if PyTorch is used."""

        if self.framework == "torch":
            LOGGER.info(f" > Using torch v{torch.__version__} ({torch.__file__})")
            return True
        return False

    def _load_model(self, device_str: str = "cuda") -> bool:
        """Loads the autoencoder model."""
        self.model = CGAutoencoder_model(
            inshape = self.inshape,
            latent_dims = self.outshape, 
            INPATH = self.model_path,
            config = self.config,
            tstamp = '0',
            create_subdirs = False
        )

        if self.framework == "torch":
            # TODO: ideally we should remove hard-coded paths
            if device_str == "cuda":
                device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
            else:
                device = torch.device('cpu')
            LOGGER.info(f" > Device used : {device}")
            if device_str == "cuda" and torch.cuda.is_available():
                LOGGER.info(f"   > GPU total devices         : {torch.cuda.device_count()}")
                LOGGER.info(f"   > GPU active device         : device:{torch.cuda.current_device()}")
                # Could lead to some errors such as:
                # Could not load library libcudnn_ops_infer.so.8. Error: libcudnn_ops_infer.so.8:
                # if torch.backends.cudnn.is_available() and torch.backends.cudnn.enabled:
                #     LOGGER.info(f"   > cudnn version              : {torch.backends.cudnn.version()}")

            self.model = self.model.to(device)
            self.model.device = device   # so the model can move the data to device later
            model_name = os.path.join(self.model_path, "model.pth")
            self.model.load_model(model_name, self.model_path, map_location = device)
        else:
            ValueError(f"Only PyTorch is supported, {self.framework} is not supported")

        return True

    def _decode(self, ls_coords: npt.NDArray) -> npt.NDArray:
        """ Generate new positions from Latent Space coordinates """
        reconstructed_positions = self.model.reconstruct_from_ls_ae(ls_coords)
        # TODO: remove magic numbers 2253, 751 etc
        # TODO: move that to mlserver.yaml
        # Pretty dirty management between C1 and PC3 models here, move to pc3/c1/c2
        pc3 = False
        c1 = True
        c2 = True
        if pc3:
            reconstructed_positions = np.reshape(reconstructed_positions, (reconstructed_positions.shape[0], 751, 3))
        if c1 or c2:
            reconstructed_positions = np.reshape(reconstructed_positions, (reconstructed_positions.shape[0], 2897, 3))

        # Apply destandardization
        new_positions = self._destandardize(
            reconstructed_positions,
            "predicted_positions",
            self.model.datamean,
            self.model.datastd
        )

        LOGGER.debug(f"new_positions= {new_positions.shape} from ls_coords = {ls_coords.shape}")
        return new_positions

    def _destandardize(self, data, label, dmean = None, dstd = None):
        """ Destandardize an array """
        LOGGER.info(f'De-standardizing {data.shape} {label}')
        data *= dstd
        data += dmean
        return data
