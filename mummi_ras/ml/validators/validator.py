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
from typing import List, Any

class Validator(ABC):
    """
    Abstract class to validate the generated CG protein structures.
    """

    @classmethod
    def __subclasshook__(cls, subclass):
        """
        Ensure subclass implement all the abstrac method 
        defined in the interface. Errors will be raised 
        if all methods aren't overridden.
        """
        return (hasattr(subclass, "validate") and
                callable(subclass.validate) and
                hasattr(subclass, "validateArray") and
                callable(subclass.validateArray) and
                hasattr(subclass, "validateArrayThreaded") and
                callable(subclass.validateArrayThreaded) or
                NotImplemented)

    @abstractmethod
    def validate(self, name, positions) -> List:
        """
        Validate one structure. Return the name of the Gromacs file
        and whether or not that file has been validated.
        """
        raise NotImplementedError

    @abstractmethod
    def validateArray(self, nameArray, positionsArray) -> np.array:
        """
        Run validate on a collection of structures by calling
        validate() on each of them.
        """
        raise NotImplementedError

    @abstractmethod
    def validateArrayThreaded(self, nameArray, positionsArray, np):
        """
        Run validate on a collection of structures by calling
        validate() on each of them in a multithreaded fashion.
        """
        raise NotImplementedError

