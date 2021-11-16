# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# This is IO interface used in MuMMI
# It provides a single interface with change-able backend
# ------------------------------------------------------------------------------


KNOWN_INTERFACES = ['mummi']


def get_interfaces():
    from mummi_core import get_interfaces
    return get_interfaces() + KNOWN_INTERFACES


def get_io(_):
    if _ == 'mummi':
        from .mummi import IO_Mummi
        interface = IO_Mummi
        interface.check_environment()
    else:
        from mummi_core import get_io
        interface = get_io(_)

    return interface

# ------------------------------------------------------------------------------
