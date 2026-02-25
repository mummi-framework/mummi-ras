#!/usr/bin/env python3
#
# Copyright (c) 2025, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# @authors
#           Tim Hsu <hsu16@llnl.gov>
#           Loic Pottier <pottier1@llnl.gov>
# ------------------------------------------------------------------------------

# Typing
import numpy.typing as npt
from typing import List, Tuple, Optional, Union, Any

from abc import ABC, abstractmethod

class GenModel(ABC):
    """Class of an abstract Generative ML model
    
    Attributes (public):
        model_path (str): Path containing trained ML model (parameter or checkpoint files).

    Methods (public):
        setup(device): Setup the model on the given device.
        generate(inputs, iteration): Transforms the inputs into outputs
    """

    def __init__(self, model_path: str) -> None:
        """Instantiates the model.

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
        return (hasattr(subclass, "setup") and
                callable(subclass.encode) and
                hasattr(subclass, "generate") and
                callable(subclass.decode) or
                NotImplemented)

    @abstractmethod
    def setup(self, device: str) -> bool:
        """This method is in charge of setting up the model,
        make sure the model is on the right device etc

        Args:
            device (str): Where to run the inference ('cpu' or 'cuda')

        Returns:
            status (bool): Return True if the setup was succesful, False otherwise
        """
        raise NotImplementedError

    @abstractmethod
    def generate(self, inputs: npt.NDArray, iterations: int) -> npt.NDArray:
        """Takes some input data in the form of batches and generates the corresponding outputs

        Args:
            inputs (npt.NDArray): Tensor of shape [B, N, D], where B is number of batchs,
                and N is number of UCG/CG beads (or the size of V) and D is some dimensions (3 for CG/UCG or 1 for v).
            iterations (int): Number of iterations to use.

        Returns:
            outputs (npt.NDArray): Return a tensor of shape [B, N', D'] where N' and D' are the new dimensions of the outputs 
        """
        raise NotImplementedError

