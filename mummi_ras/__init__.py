# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
from mummi_core.utils.timer import Timer
from mummi_core.utils.logger import init_logger
from mummi_core.utils.naming import MuMMI_NamingUtils as Naming

from .interfaces import get_interfaces, get_io


# ------------------------------------------------------------------------------
# lipid order
# ------------------------------------------------------------------------------
LIPID_TYPES = ['POPC', 'PAPC', 'POPE', 'DIPE', 'DPSM', 'PAPS', 'PAP6', 'CHOL']
NLIPIDS = 14

# -----------------------------------------------------------------------------
