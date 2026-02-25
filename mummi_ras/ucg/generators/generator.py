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

from abc import ABC, abstractmethod
import numpy as np
from typing import List, Union, Any

class Generator(ABC):
    """
    Abstract class to generate protein structures from Latent Space samples.
    """
    @classmethod
    def __subclasshook__(cls, subclass):
        """
        Ensure subclass implement all the abstrac method 
        defined in the interface. Errors will be raised 
        if all methods aren't overridden.
        """
        return (hasattr(subclass, "generate") and
                callable(subclass.generate) or
                NotImplemented)

    @abstractmethod
    def new_positions(self, coords: np.array) -> np.array:
        """ Generate new positions from Latent Space / UCG coordinates """
        raise NotImplementedError

    @abstractmethod
    def write_patches(self, iteration_id: int, new_positions: np.array) -> Union[List[str], List[Any]]:
        """
        From positions, write a patch as a npz file and 
        return its path.
        """
        raise NotImplementedError