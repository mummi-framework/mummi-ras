#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------

from enum import Enum
from .simulation import Simulation, NoSimulation, ddcMD, GROMACS, GROMACS_PARTS

class MDEngine(Enum):
    NOSIM = NoSimulation
    DDCMD = ddcMD
    GROMACS = GROMACS
    GROMACS_PARTS = GROMACS_PARTS

    def __str__(self):
        return self.name

    @staticmethod
    def from_string(s: str):
        try:
            return MDEngine[s]
        except KeyError:
            raise ValueError()
