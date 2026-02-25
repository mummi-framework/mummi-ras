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

import logging

import torch
from .model_torch import BaseModel
# from model_ae.model_encoder_torch_cons import Encoder
from .model_encoder_torch import Encoder
from .model_decoder_torch_cons_pos import Decoder_pos
# from model_ae.model_decoder_torch_cons_dist import Decoder_dist
# from utils.dataloaders_admirral import produce_hierarchical_distances
import numpy as np
import os
# from model_ae.potentials import potential_loss

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
        self.decoder_pos = Decoder_pos(latent_dims, 2253)
        # self.decoder_dist = Decoder_dist(latent_dims, [15328])

        std_fname_pos = os.path.join(self.INPATH, 'CG_data_summary_pos_dis.npz')

        with np.load(std_fname_pos) as npz:
            self.mean_pos = torch.from_numpy(npz['mean_pos'])
            self.std_pos = torch.from_numpy(npz['std_pos'])

        self.loss_fn = torch.nn.MSELoss()

    def _compute_losses(self, data, training=True):

        x = data[0]
        pos_pred, dis_pred = self._reconstruct(x, training=training)

        pos_real = x[:,:2253]

        # LJ loss
        pos_pred_3d = pos_pred.reshape(pos_pred.shape[0], -1, 3)
        dis_from_pos = torch.cdist(pos_pred_3d, pos_pred_3d)
        loss_potential = potential_loss(dis_from_pos)
        
        dis_real = x[:,2253:]

        loss_dis = torch.sqrt(self.loss_fn(dis_real, dis_pred))

        loss_pos = torch.sqrt(self.loss_fn(pos_real, pos_pred))

        # For consistency loss of distances
        # # pos_pred = np.array(pos_pred.tolist())
        # # pos_pred = np.reshape(pos_pred, (pos_pred.shape[0], 751, 3))
        # # # destandardize positions
        # # std_fname_pos = os.path.join(self.INPATH, 'CG_data_summary_pos_dis.npz')
    
        # # with np.load(std_fname_pos) as npz:
        # #     mean_pos = npz['mean_pos']
        # #     std_pos = npz['std_pos']
        # #     pos_pred *= std_pos
        # #     pos_pred += mean_pos

        # # pos_pred = np.array(pos_pred.tolist())
        # # pos_pred = np.reshape(pos_pred, (pos_pred.shape[0], 751, 3))
        # pos_pred = pos_pred.reshape(pos_pred.shape[0], 751, 3)
        # # destandardize positions
        # pos_pred = pos_pred.cpu() * self.std_pos
        # pos_pred = pos_pred.cpu() + self.mean_pos

        # pos_pred = np.array(pos_pred.tolist())
        # dist_by_pred_pos = produce_hierarchical_distances(self.INPATH, pos_pred)

        # dist_by_pred_pos = torch.tensor(dist_by_pred_pos, dtype=torch.float32).to(self.device) 
        # loss_cons = self.loss_fn(dis_pred, dist_by_pred_pos) 

        # loss_value = loss_dis + loss_pos + 0.001 * loss_cons

        # loss_value = loss_dis + 2 * loss_pos
        # loss_value = loss_dis + loss_pos
        # print(f"kaka be {loss_value} -- {loss_potential}")
        if torch.isfinite(loss_potential):
            loss_potential = 1 * loss_potential #0.1 is arbitrary weighing term since potential loss is somewhat larger even with log
        # loss_value = loss_dis + 2*loss_pos + 1*loss_potential
        loss_value = loss_dis + (2 * loss_pos) + loss_potential
        # loss_value = loss_potential
        # loss_value = loss_pos
        # print(f"kaka aft {loss_value} -- {loss_potential}")
        # loss_report = {'loss-total': loss_value, 'loss-pos': loss_pos,
        #                 'loss-dis': loss_dis, 'loss-cons':loss_cons}

        # print(f"{loss_value} -- {loss_potential} {loss_dis} {loss_pos}")
        # loss_report = {'loss-total': loss_value, 'loss-pos': loss_pos,
        #                 'loss-dis': loss_dis}
        loss_report = {'loss-total': loss_value, 'loss-pos': loss_pos,
                        'loss-dis': loss_dis, 'loss-potential': loss_potential}
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
        return self.decoder_pos(z)
    