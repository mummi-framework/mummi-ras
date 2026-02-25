#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2023, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------

import logging

import torch
from .model_torch import BaseModel
from .model_encoder_torch import Encoder
from .model_decoder_torch import Decoder_pos
import numpy as np
import os
LOGGER = logging.getLogger(__name__)


# ------------------------------------------------------------------------------
# A autoencoder model
# ------------------------------------------------------------------------------
class Autoencoder(BaseModel):
    """An Abstract Model."""

    # --------------------------------------------------------------------------
    def __init__(self, inshape, latent_dims, INPATH, config, tstamp, create_subdirs, **kwargs):
        super().__init__(config, tstamp, create_subdirs, **kwargs)

        assert self.type == 'protein_latent_space', f'Invalid type = {self.type}'

        self.input_shp = inshape
        self.output_shp = latent_dims
        self.show_plots = False
        self.INPATH = INPATH

        # architectures can be defined in architectures.py
        # and chosen using config file
        # The mini model has 763 beads
        self.encoder = Encoder(inshape, latent_dims)
        self.decoder = Decoder_pos(latent_dims, inshape)
        self.loss_fn = torch.nn.MSELoss()

    def _encode(self, data, training=False):

        z = self.encoder(data)

        if training:
            return z
        else:
            return z.detach()


    def _reconstruct(self, data, training=False):
        pos, dis = self.forward(data)
        if training:
            return pos, dis
        else:
            pos = pos.detach()
            dis = dis.detach()
            return pos, dis

# ------------------------------------------------------------------------------
# a custom autoencoder model
# ------------------------------------------------------------------------------
class CustomAutoencoder(Autoencoder):

    def __init__(self, inshape, latent_dims, INPATH, config, tstamp, create_subdirs):
        super(CustomAutoencoder,self).__init__(inshape, latent_dims,  INPATH, config, tstamp, create_subdirs)

    def forward(self, x):
        if not isinstance(x, torch.Tensor):
            x = torch.Tensor(x)

        x = x.to(self.device)
        z = self.encoder(x)
        # d = self.decoder_dist(z)
        # r = torch.cat((z, d), 1)
        return self.decoder(z)
