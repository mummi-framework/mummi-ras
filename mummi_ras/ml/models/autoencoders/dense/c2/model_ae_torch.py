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
        #arch = getattr(architectures, self.arch_name)
        #self.encoder = arch(self.input_shp, self.z_dim)
        self.encoder = Encoder(inshape, latent_dims)
        # 2897 beads * 3 dimensions = 8691
        self.decoder = Decoder_pos(latent_dims, 8691)
        # self.decoder_dist = Decoder_dist(latent_dims, [15328])

        std_fname_pos = os.path.join(self.INPATH, 'CG_pos_data_summary_pos_dis_C2.npz')

        with np.load(std_fname_pos, allow_pickle=True) as npz:
            self.mean_pos = torch.from_numpy(npz['mean_pos'])
            self.std_pos = torch.from_numpy(npz['std_pos'])

        # norm_fname_pos = os.path.join(INPATH, 'CG_pos_data_summary_pos_min_max_v1.npz')

        # with np.load(norm_fname_pos, allow_pickle=True) as npz:
        #     self.min = torch.from_numpy(npz['min'])
        #     self.max = torch.from_numpy(npz['max'])

        indices_fname = os.path.join(INPATH, 'indices_per_cluster.npz')
        # Load the CG indices per UCG site
        with np.load(indices_fname, allow_pickle=True) as npz:
            self.ids_per_cluster = npz['indices']

        bbids_fname_ = os.path.join(INPATH, 'bb_ids.npz')
        with np.load(bbids_fname_) as npz:
            self.bb_ids = torch.from_numpy(npz['bb_ids'])

        self.loss_fn = torch.nn.MSELoss()

    def _compute_losses(self, data, training=True):

        x = data[0]
        pos_pred, pos_ucg_pred = self._reconstruct(x, training=training)

        pos_real = x[:,:8691] # 2897 beads * 3 dimensions = 8691

        pos_real_3d = pos_real.reshape(pos_real.shape[0], -1, 3)
        pos_bb_real = pos_real_3d[:,self.bb_ids,:]
        pos_bb_real = torch.flatten(pos_bb_real, start_dim=1)

        pos_pred_3d = pos_pred.reshape(pos_pred.shape[0], -1, 3)
        pos_bb_pred = pos_pred_3d[:,self.bb_ids,:]
        pos_bb_pred = torch.flatten(pos_bb_pred, start_dim=1)
        loss_pos_bb = torch.sqrt(self.loss_fn(pos_bb_real, pos_bb_pred))

        loss_pos = torch.sqrt(self.loss_fn(pos_real, pos_pred))
        loss_value = (loss_pos) + (2*loss_pos_bb)
        loss_report = {'loss-total': loss_value, 'loss-pos': loss_pos,
                        'loss-loss_pos_bb': loss_pos_bb}
        # loss_report = {'loss-total': loss_value, 'loss-pos': loss_pos}
        return loss_value, loss_report

    def _encode(self, data, training=False):

        z = self.encoder(data)

        if training:
            return z
        else:
            return z.detach()


    def _reconstruct(self, data, training=False):
        # _ = self.forward(data)
        # return _ if training else _.detach()

        pos, dis = self.forward(data)
        # return _ if training else _.detach()
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
