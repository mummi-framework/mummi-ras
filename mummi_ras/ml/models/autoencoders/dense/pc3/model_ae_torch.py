#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------

import logging

import torch
from .model_torch import BaseModel
from .model_encoder_torch import Encoder
from .model_decoder_torch import Decoder

LOGGER = logging.getLogger(__name__)


# ------------------------------------------------------------------------------
# A autoencoder model
# ------------------------------------------------------------------------------
class Autoencoder(BaseModel):
    """An Abstract Model."""

    # --------------------------------------------------------------------------
    def __init__(self, inshape, latent_dims, config, tstamp, create_subdirs, **kwargs):
        super().__init__(config, tstamp, create_subdirs, **kwargs)

        assert self.type == 'protein_latent_space', f'Invalid type = {self.type}'

        self.input_shp = inshape
        self.output_shp = latent_dims
        self.show_plots = False

        # architectures can be defined in architectures.py
        # and chosen using config file
        #arch = getattr(architectures, self.arch_name)
        #self.encoder = arch(self.input_shp, self.z_dim)
        self.encoder = Encoder(inshape, latent_dims)
        self.decoder = Decoder(inshape, latent_dims)

        self.loss_fn = torch.nn.MSELoss()

    def _compute_losses(self, data, training=True):

        x = data[0]
        x_pred = self._encode(x, training=training)

        loss_value = self.loss_fn(x, x_pred)
        loss_report = {'loss': loss_value}
        return loss_value, loss_report

    def _encode(self, data, training=False):
        _ = self.forward(data)
        return _ if training else _.detach()

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# a custom autoencoder model
# ------------------------------------------------------------------------------
class CustomAutoencoder(Autoencoder):

    def __init__(self, inshape, latent_dims, config, tstamp, create_subdirs):
        super(CustomAutoencoder,self).__init__(inshape, latent_dims, config, tstamp, create_subdirs)


    def forward(self, x):
        if not isinstance(x, torch.Tensor):
            x = torch.Tensor(x)

        x = x.to(self.device)
        z = self.encoder(x)
        return self.decoder(z)
    