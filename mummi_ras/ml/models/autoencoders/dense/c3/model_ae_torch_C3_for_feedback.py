# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
import logging

import torch
import torch.nn as nn
# from model_base.model_torch_C3_for_feedback import BaseModel
from .model_torch import BaseModel
import numpy as np
import os
# from model_ae.rmsd_loss import rmsd
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
        self.loss_fn = torch.nn.MSELoss()

        self.encoder = nn.Sequential(
            nn.Linear(31125, 512),  # 256 common part has sub_beads, each with 3 dimension
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Linear(512, 256),  # 256 common part has sub_beads, each with 3 dimension
            nn.BatchNorm1d(256),
            nn.ReLU(), 
             
            nn.Linear(256, latent_dims[0]), #49 80
        )

        self.decoder = nn.Sequential(
            nn.Linear(latent_dims[0],256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Linear(256, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Linear(512, 31125)
        )

    def normalize_distances(self, distances):
        min_val = distances.min()
        max_val = distances.max()
        return (distances - min_val) / (max_val - min_val)

    def _compute_losses(self, data, training=True):

        structures = data[0] 
        pairwise_distances = torch.cdist(structures, structures)
        
        # Get the upper triangular part (excluding diagonal)
        upper_tri_indices = torch.triu_indices(250, 250, offset=1)
        x = pairwise_distances[:, upper_tri_indices[0], upper_tri_indices[1]]  # Shape: (10, num_upper_tri)
        x_pred = self._reconstruct(x, training=training)
        loss_reconstruction = torch.sqrt(self.loss_fn(x,x_pred))

        loss_value = loss_reconstruction
        loss_report = {'loss-total': loss_value, 'loss_reconstruction': loss_reconstruction}
        return loss_value, loss_report

    def _encode(self, data, training=False):

        z = self.encoder(data)

        if training:
            return z
        else:
            return z.detach()


    def _reconstruct(self, data, training=False):

        dis = self.forward(data)
        if training:
            return dis
        else:
            dis = dis.detach()
            return dis

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
        return self.decoder(z)
