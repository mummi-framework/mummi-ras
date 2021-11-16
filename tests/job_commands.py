# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
#!/usr/bin/env python

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
import os
import glob
import numpy as np
import yaml
from queue import Queue
import logging


import mummi_core
import mummi_ras
from mummi_core import Naming

LOGGER = logging.getLogger(__name__)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def do_createsim(spath):

    from mummi_ras.workflow.jobTracker_createsim import CreateSimTracker

    rpath = os.path.join(spath, 'createsim.yaml')
    with open(rpath, 'r') as data:
        csim = yaml.load(data, Loader=yaml.FullLoader)
        csim['createsim']['config']['total_nnodes'] = 100

    jt = CreateSimTracker(csim['createsim'], wpath, None, 'Flux', {}, Queue())
    return jt.command(['pfpatch_0000000000'])


def do_ddcmd(spath):

    from mummi_ras.workflow.jobTracker_cg import CGTracker

    rpath = os.path.join(spath, 'cg.yaml')
    with open(rpath, 'r') as data:
        cg = yaml.load(data, Loader=yaml.FullLoader)
        cg['cg']['config']['total_nnodes'] = 100

    jt = CGTracker(cg['cg'], wpath,  None, 'Flux', {}, Queue())
    s = jt.command(['pfpatch_0000000000']) #,'pfpatch_0000000001','pfpatch_0000000002','pfpatch_0000000003'])

    return s


def do_aa(spath):

    from mummi_ras.workflow.jobTracker_aa import AATracker

    rpath = os.path.join(spath, 'aa.yaml')
    with open(rpath, 'r') as data:
        cg = yaml.load(data, Loader=yaml.FullLoader)
        cg['aa']['config']['total_nnodes'] = 100

    jt = AATracker(cg['aa'], wpath,  None, 'Flux', {}, Queue())
    s = jt.command([Naming.cgframe('pfpatch_0000000000', 0)])

    return s


def do_backmapping(spath):

    from mummi_ras.workflow.jobTracker_backmapping import BackmappingTracker

    rpath = os.path.join(spath, 'backmapping.yaml')
    with open(rpath, 'r') as data:
        csim = yaml.load(data, Loader=yaml.FullLoader)
        csim['backmapping']['config']['total_nnodes'] = 100

    patch = 'pfpatch_0000000000'
    frame = 99

    jt = BackmappingTracker( csim['backmapping'], wpath,  None, 'Flux', {}, Queue())
    return jt.command([Naming.cgframe(patch, frame)])


# ------------------------------------------------------------------------------
if __name__ == '__main__':

    mummi_core.init()
    mummi_core.create_root()
    mummi_core.init_logger(level=1)

    wpath = './'
    spath = os.path.join(Naming.MUMMI_SPECS, 'workflow')

    #s = do_createsim(spath)
    #s = do_backmapping(spath)
    #s = do_ddcmd(spath)
    s = do_aa(spath)

    print ('\n')
    print ('----------------------------------------')
    print (s)
    print ('----------------------------------------')
    print ('\n')

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
