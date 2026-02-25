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
from numpy import ndarray
from torch import Tensor
from typing import Tuple, List, Optional

####################### Base modules #######################

from torch import nn
import torch_sparse

class MLP(nn.Module):
    """Multi-layer perceptron with custom activation functions.
    """
    def __init__(self, dims: List[int], act=None) -> None:
        """
        Args:
            dims (List[int]): Input, hidden, and output dimensions.
            act (torch activation function, optional): Activation function that
                applies to all but the output layer. For example, 'nn.ReLU()'.
                If None, no activation function is applied.
        """
        super().__init__()
        self.dims = dims
        self.act = act
        
        num_layers = len(dims)

        layers = []
        for i in range(num_layers-1):
            layers += [nn.Linear(dims[i], dims[i+1])]
            if (act is not None) and (i < num_layers-2):
                layers += [act]

        self.mlp = nn.Sequential(*layers)

    def forward(self, x: Tensor) -> Tensor:
        return self.mlp(x)
    
    def __repr__(self) -> str:
        return f'{self.__class__.__name__}(dims={self.dims}, act={self.act})'


class BlockDiagonalLinear(nn.Module):
    """Sparse linear layer with block diagonal matrix.
    """
    def __init__(self, in_dims: List[int], out_dims: List[int]) -> None:
        """
        Args:
            in_dims (List[int]): List of input block dimensions.
            out_dims (List[int]): List of output block dimensions.
        """
        super().__init__()
        
        assert len(in_dims) == len(out_dims)

        self.in_dims  = in_dims
        self.out_dims = out_dims

        blocks = []
        for nrows, ncols in zip(out_dims, in_dims):
            blocks.append(torch.ones(nrows, ncols))

        sparse_weight  = torch.block_diag(*blocks).to_sparse()
        self.m, self.n = sparse_weight.shape
        self.indices   = sparse_weight.indices()
        self.values    = nn.Parameter(sparse_weight.values()) 
        self.bias      = nn.Parameter(torch.randn(sum(out_dims)))

        self.reset_parameters()

    def reset_parameters(self) -> None:
        torch.nn.init.normal_(self.values)

    def forward(self, input: Tensor) -> Tensor:
        return torch_sparse.spmm(self.indices.to(input.device), self.values, self.m, self.n, input.T).T + self.bias

    def extra_repr(self) -> str:
        return f'in_dims={self.in_dims}, out_dims={self.out_dims}, bias={self.bias is not None}'


class BlockDiagonalMLP(nn.Module):
    """MLP made of block diagonal linear matrices.
    """
    def __init__(self, dims: ndarray, act=None) -> None:
        """
        Args:
            dims (ndarray): 2D int array of shape (N, M), where N is number of blocks
                and M is number of (block diagonal) linear layers.
            act (torch activation function, or None): activation function that
                applies to all but the output layer. For example, 'torch.nn.ReLU()'.
                If None, no activation function is applied.
        """
        super().__init__()

        assert dims.ndim == 2

        self.dims = dims
        self.act  = act

        num_layers = dims.shape[1]

        layers = []
        for i in range(num_layers-1):
            in_dims, out_dims = dims[:, i:i+2].T
            layers.append(BlockDiagonalLinear(in_dims, out_dims))
            if (act is not None) and (i < num_layers-2):
                layers.append(act)

        self.mlp = nn.Sequential(*layers)

    def forward(self, x: Tensor) -> Tensor:
        return self.mlp(x)

####################### Model #######################

from torch import nn
from torch_scatter import scatter
import yaml
import csv

class HierarchicalAutoencoder(nn.Module):
    def __init__(self, model_config_fname: str, ucg_fname: str = None, scale: float = 0.1, bb_bond_index: ndarray = None) -> None:
        super().__init__()

        # Read model config
        with open(model_config_fname, 'r') as f:
            model_config = yaml.safe_load(f)
        local_encoder_dims  = model_config['ucg_encoder_dims']
        local_encoder_dims  = np.array(local_encoder_dims)
        local_decoder_dims  = np.flip(local_encoder_dims, axis=1)
        global_encoder_dims = model_config['global_encoder_dims']
        global_decoder_dims = global_encoder_dims[::-1]

        # Local autoencoders
        self.local_encoder = BlockDiagonalMLP(local_encoder_dims, act=nn.SiLU())
        self.local_decoder = BlockDiagonalMLP(local_decoder_dims, act=nn.SiLU())

        # Global autoencoder
        self.global_encoder = MLP(global_encoder_dims, act=nn.SiLU())
        self.global_decoder = MLP(global_decoder_dims, act=nn.SiLU())
        
        # Parameters stored for later use
        self.scale = scale
        self.ucg_x_dim = local_encoder_dims[0, 0].item() // 3
        self.bb_bond_index = bb_bond_index

        self.ucg_fname = ucg_fname
        if self.ucg_fname is not None:
            with open(ucg_fname, 'r') as f:
                ucg_idx = [np.array(indices, dtype=int) for indices in csv.reader(f)]
                self.scatter_idx = self._ucg_idx_to_scatter_idx(ucg_idx)
                self.flattened_ucg_idx = np.concatenate(ucg_idx)
                self.inverse_idx = np.argsort(self.flattened_ucg_idx)
        
    def _ucg_idx_to_scatter_idx(self, ucg_idx: List[List[int]]) -> Tensor:
        N = len(np.concatenate(ucg_idx))
        scatter_idx = torch.zeros(N, dtype=torch.long)
        for i, idx in enumerate(ucg_idx):
            scatter_idx[idx] = i
        return scatter_idx

    def _x2pos(self, x: Tensor) -> Tensor:
        bsize = x.size(0)
        x = x.view(bsize, -1, 3)
                
        # Convert flattened 1D array 'x' to 3D coordinates 'pos'
        ucg_pos = x[:, :self.ucg_x_dim, :]
        pos     = x[:, self.ucg_x_dim:, :]
        pos     = pos[:, self.inverse_idx, :]
        pos    += ucg_pos[:, self.scatter_idx, :]
        return pos / self.scale

    def _pos2x(self, pos: Tensor) -> Tensor:
        bsize = pos.size(0)

        # Center and scale the positions
        mean_pos = pos.mean(dim=1, keepdim=True)
        mean_pos[:, :, -1] = 0.0 # do not center z-axis
        pos -= mean_pos
        pos *= self.scale

        # Convert 3D coordinates 'pos' to flattened 1D array 'x'
        ucg_pos = scatter(src=pos, index=self.scatter_idx, dim=1, reduce='mean')
        pos -= ucg_pos[:, self.scatter_idx, :]
        pos  = pos[:, self.flattened_ucg_idx, :]
        # print(ucg_pos.shape, pos.shape)
        return torch.concatenate([ucg_pos, pos], dim=1).view(bsize, -1)

    def encode(self, x: Tensor) -> Tuple[Tensor, Tensor]:
        z   = self.local_encoder(x)
        l   = self.global_encoder(z)
        return l, z

    def decode(self, l: Tensor) -> Tuple[Tensor, Tensor]:
        z_ = self.global_decoder(l)
        x_ = self.local_decoder(z_)
        return x_, z_