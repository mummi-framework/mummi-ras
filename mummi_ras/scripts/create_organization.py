# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
import os
import sys

from mummi_ras import Naming

if __name__ == '__main__':

    nargs = len(sys.argv)
    if nargs < 1 or nargs > 3:
        print('Usage:', sys.argv[0], 'root (default: MUMMI_ROOT) resources (default: MUMMI_RESOURCES)')
        exit()

    MUMMI_ROOT = ''
    MUMMI_RES = ''

    if nargs == 2 or nargs == 3:
        MUMMI_ROOT = sys.argv[1]

    if nargs == 3:
        MUMMI_RES = sys.argv[2]

    Naming.init(MUMMI_ROOT, MUMMI_RES)
    Naming.create_root()
