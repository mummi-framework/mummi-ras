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
from mummi_ras.ucg import UCGServerRPC

LOGGER = getLogger(__name__)

# ------------------------------------------------------------------------------
# read config file and command line arguments to get necessary settings
# ------------------------------------------------------------------------------
def read_specs(spath):
    # read the specs    
    with open(os.path.join(spath, "ucgserver.yaml"), 'r') as data:
        config = yaml.load(data, Loader = yaml.FullLoader)

    if config.get("encoder") is None or config["encoder"].get("path") is None:
        raise ValueError("No encoder specified")
    if config.get("workspace") is None or config["workspace"].get("path") is None:
        raise ValueError("No workspace specified")

    # Load paths from eval strings in config
    gdict = {"mummi_ras": mummi_ras}
    for key in ["workspace", "encoder"]:
        if type(config[key]["path"]) is dict and "eval" in config[key]["path"]:
            config[key]["path"] = eval(config[key]["path"]["eval"], gdict)

    for key in ["generator", "validator"]:
        for val in ["inpath", "outpath"]:
            if config[key].get(val) and type(config[key][val]) is dict and "eval" in config[key][val]:
                config[key][val] = eval(config[key][val]["eval"], gdict)

    return config

# -----------------------------------------------------------------------------
def main():
    mummi_core.init()
    mummi_core.create_root()
    if Naming.MUMMI_MINI:
        config = read_specs(os.path.join(Naming.MUMMI_SPECS, "mini-mummi"))
    else:
        config = read_specs(os.path.join(Naming.MUMMI_SPECS, "workflow"))
    mummi_core.init_logger(config = config['config'])

    try:
        server = UCGServerRPC(config = config, logger = LOGGER)
        server.setup()
        server.run()
    except Exception as e:
        LOGGER.error(f"> Exiting UCG server due to error ({e})")
        traceback.print_exc()
        exit(1)

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()