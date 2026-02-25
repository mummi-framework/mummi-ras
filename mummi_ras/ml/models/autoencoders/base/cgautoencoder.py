#!/usr/bin/env python3
#
# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# @authors
#           Loic Pottier <pottier1@llnl.gov>
#           Tim Hsu <hsu16@llnl.gov>
#           Konstantia Georgouli <georgouli1@llnl.gov>
# ------------------------------------------------------------------------------

import torch
import numpy as np

# Typing
from torch import Tensor
import numpy.typing as npt
from typing import List, Tuple, Optional, Union, Any
from mummi_ras.ml import ls_point as lsp

from abc import ABC, abstractmethod

class CGAutoencoder(ABC):
    """Class of an abstract CG autoencoder 
    
    Attributes (public):
        model_path (str): Path containing trained ML model (parameter or checkpoint files).

    Methods (public):
        encode(pos): Encodes CG positions into latent space coordinates represented as LSPoint.
        decode(l): Decodes latent space coordinates (LSPoint) into CG positions.
    """

    def __init__(self, model_path: str) -> None:
        """Instantiates the autoencoder.

        Args:
            model_path (str): Path containing the trained ML model (parameter or checkpoint files).
        """
        self.model_path  = model_path

    @classmethod
    def __subclasshook__(cls, subclass):
        """Ensure subclass implement all the abstract method 
        defined in the interface. Errors will be raised 
        if all methods aren't overridden.
        """
        return (hasattr(subclass, "encode") and
                callable(subclass.encode) and
                hasattr(subclass, "decode") and
                callable(subclass.decode) or
                NotImplemented)

    @abstractmethod
    def encode(self, pos: npt.NDArray) -> npt.NDArray:
        """Encodes CG positions into latent space.

        Args:
            pos (npt.NDArray): CG positions with shape [B, N, 3], where B is number of batchs,
                and N is number of CG beads.

        Returns:
            l (npt.NDArray): Latent space coordinates with shape [B, D], where B is number of batchs,
                and D is latent space dimensionality.
        """
        raise NotImplementedError

    @abstractmethod
    def decode(self, l: npt.NDArray[lsp]) -> npt.NDArray:
        """Decodes latent space coordinates into CG positions.

        Args:
            l (npt.NDArray[lsp]): Latent space coordinates with shape [B, D], where B is number of batchs,
                and D is latent space dimensionality.

        Returns:
            pos (npt.NDArray): CG positions with shape [B, N, 3], where B is number of batchs,
                and N is number of CG beads.
        """
        raise NotImplementedError
