#!/usr/bin/env python3

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

def run_loop(i):
    
    name = Naming.pfpatch(multiprocessing.current_process().pid)
    interface_io = mummi_ras.get_io('mummi')
    interface_fb = mummi_ras.get_io('mummi')
    wspace_fb = Naming.dir_root('feedback-cg')
    nspace_fb = Naming.dir_root('feedback-cg')

    fb = MacroFeedbackManager(FeedbackManagerType.Worker, name, wspace_fb,
                              wspace_fb, nspace_fb, interface_io, interface_fb,
                              None)

    rdfs = np.ones(RDF.SHAPE, dtype=RDF.DTYPE)
    LOGGER.info(f"Start parallel job {name}")
    nSecPerFeedback = 2
    i = 0
    while True:
        LOGGER.info(f"Fake aggregate {i}")
        fb.aggregate(rdfs)
        fb.report()
        fb.reset()

        time.sleep(nSecPerFeedback + random.uniform(0, 1))
        i += 1

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if __name__ == '__main__':

    # --------------------------------------------------------------------------
    mummi_core.init()
    mummi_core.create_root()
    mummi_core.init_logger(level=1)

    LOGGER.info(f"Starting fake cg analysis on host({mummi_core.get_hostname()})")

    # --------------------------------------------------------------------------
    # perform fake feedback
    nprocs = 100

    pool = Pool(processes=nprocs)
    pool.map(run_loop, range(nprocs))
    pool.close()
    pool.join()

    # --------------------------------------------------------------------------
# ------------------------------------------------------------------------------
