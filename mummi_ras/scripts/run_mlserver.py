#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------

import os
import yaml
import mummi_ras
import traceback
from logging import getLogger

import mummi_core
from mummi_ras import Naming
from mummi_ras.ml import MLServerRPC

LOGGER = getLogger(__name__)

# -----------------------------------------------------------------------------
def main():
    mummi_core.init()
    mummi_core.create_root()
    if Naming.MUMMI_MINI:
        spath = os.path.join(Naming.MUMMI_SPECS, "mini-mummi")
    else:
        spath = os.path.join(Naming.MUMMI_SPECS, "workflow")

    config = mummi_core.read_specs(os.path.join(spath, "mlserver.yaml"))
    mummi_core.init_logger(config = config["config"])

    lockfile = os.path.join(config["config"]["lockfile"]["path"], config["config"]["lockfile"]["name"])

    try:
        server = MLServerRPC(config = config, lockfile = lockfile, logger = LOGGER)
        server.setup()
        server.run()
    except Exception as e:
        LOGGER.error(f"> Exiting ML server due to error ({e})")
        traceback.print_exc()
        exit(1)

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()