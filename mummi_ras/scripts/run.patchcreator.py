# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
import os
import yaml
from queue import Queue

import mummi_core
import mummi_ras
from mummi_core.utils import Naming
from mummi_ras.datastructures.patch import PatchCreatorManager
from core.workflow import PatchCreatorManager
from mummi_core.utils.logger import init_logger

from logging import getLogger
LOGGER = getLogger(__name__)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if __name__ == '__main__':

    hostname = mummi_ras.get_hostname(contract_hostname=False)
    print ('hostname =', hostname)

    mummi_core.init()
    mummi_core.create_root()

    # --------------------------------------------------------------------------
    # read the specs
    spath = os.path.join(Naming.MUMMI_SPECS, 'workflow/pf2pfpatches.yaml')
    with open(spath, 'r') as data:
        specs = yaml.load(data, Loader=yaml.FullLoader)

    # --------------------------------------------------------------------------
    # create logger
    mummi_core.init_logger(config = specs['macro_patch_creator']['config'])

    # --------------------------------------------------------------------------
    # start the job

    env = None
    exceptions = Queue()
    job = PatchCreatorManager('MacroPatchCreator', exceptions)
    job.setup(env, specs)
    job.run()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
