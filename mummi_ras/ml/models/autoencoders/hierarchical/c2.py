#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2023, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------

import torch
import numpy as np

# Typing
from torch import Tensor
import numpy.typing as npt
from typing import List, Tuple, Optional, Union, Any
from mummi_ras.ml import ls_point as lsp

from .base import HierarchicalAutoencoder
from ..base.cgautoencoder import CGAutoencoder

class HierarchicalCGAutoencoder(CGAutoencoder):
    def __init__(self, model_path: str) -> None:
        """Instantiates the autoencoder.

        Args:
            model_path (str): Path containing the trained ML model (parameter or checkpoint files).
        """
        from pathlib import Path
        self.model_path = Path(model_path)

        # Instantiate model
        model_config_fname = list(self.model_path.glob('*.yaml'))[0]
        ucg_fname          = self.model_path/'ucg-idx.txt'
        bb_bond_fname      = self.model_path/'bb_bonds_dis_raw.npz'
        self.autoencoder = HierarchicalAutoencoder(
            model_config_fname = model_config_fname,
            ucg_fname          = ucg_fname,
            bb_bond_index      = np.load(bb_bond_fname)['indices'].astype(int).T
        )

        # Load model parameters
        ckpt_path = list(self.model_path.glob('*.ckpt'))[0]
        ckpt_dict = torch.load(ckpt_path)
        self.autoencoder.load_state_dict(ckpt_dict['state_dict'])
        self.autoencoder.to('cpu')
        self.autoencoder.eval()

    def encode(self, pos: npt.NDArray) -> npt.NDArray[lsp]:
        """Encodes CG positions into latent space.

        Args:
            pos (npt.NDArray): CG positions with shape [B, N, 3], where B is number of batchs,
                and N is number of CG beads.

        Returns:
            l (npt.NDArray): Latent space coordinates with shape [B, D], where B is number of batchs,
                and D is latent space dimensionality.
        """
        with torch.no_grad():
            pos = torch.tensor(pos, dtype=torch.float)
            x   = self.autoencoder._pos2x(pos)
            l, z = self.autoencoder.encode(x)
            return l.numpy()

    def decode(self, l: npt.NDArray[lsp]) -> npt.NDArray:
        """Decodes latent space coordinates into CG positions.

        Args:
            l (npt.NDArray[lsp]): Latent space coordinates with shape [B, D], where B is number of batchs,
                and D is latent space dimensionality.

        Returns:
            pos (npt.NDArray): CG positions with shape [B, N, 3], where B is number of batchs,
                and N is number of CG beads.
        """
        with torch.no_grad():
            l = torch.tensor([pts.get_coordinates() for pts in l], dtype=torch.float)
            x_, z_ = self.autoencoder.decode(l)
            pos_ = self.autoencoder._x2pos(x_)
            return pos_.numpy()