# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
import sys
import time
import numpy as np
import random
import multiprocessing
import logging
from multiprocessing import Pool

import mummi_core, mummi_ras
from mummi_core.utils import timeout, Naming
from mummi_ras.datastructures.rdf import RDF
from mummi_ras.feedback import MacroFeedbackManager, FeedbackManagerType

LOGGER = logging.getLogger(__name__)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if __name__ == '__main__':

    # --------------------------------------------------------------------------
    mummi_core.init()
    mummi_core.create_root()
    mummi_core.init_logger(level=1)

    LOGGER.info(f"Starting macro feedback on host({mummi_core.get_hostname()})")

    interface_io = mummi_ras.get_io('mummi')
    interface_fb = mummi_ras.get_io('mummi')
    wspace_fb = Naming.dir_root('feedback-cg')
    nspace_fb = Naming.dir_root('feedback-cg')

    fb_cg2macro = MacroFeedbackManager(FeedbackManagerType.Manager,
                                                    'test_manager', wspace_fb, wspace_fb, nspace_fb,
                                                    interface_io, interface_fb,
                                                    None)
    fb_cg2macro.restore()
    fb_cg2macro.aggregate()
    fb_cg2macro.report()
    fb_cg2macro.checkpoint()
    

    # --------------------------------------------------------------------------
# ------------------------------------------------------------------------------
