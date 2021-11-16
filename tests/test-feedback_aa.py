# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
#!/usr/bin/env python

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

import os
import sys
import yaml
import time
import logging

#sys.path=['/ccs/home/fikret/mummi/CODE_ROOT'] + sys.path

import mummi_core
import mummi_ras
from mummi_core import Naming
from mummi_ras.feedback import FeedbackManager_AA2CG, FeedbackManagerType

LOGGER = logging.getLogger(__name__)

# ------------------------------------------------------------------------------
if __name__ == '__main__':

    print (f'importing mummi_core from ({mummi_core.__file__})')
    print (f'importing mummi_ras from ({mummi_ras.__file__})')

    mummi_core.init()
    mummi_core.create_root()
    mummi_core.init_logger(level=1)

    # --------------------------------------------------------------------------
    iointerface = mummi_ras.get_io('mummi')
    fbinterface = mummi_ras.get_io('mummi')
    workspace = Naming.dir_root('feedback-aa')
    outpath = Naming.dir_root('feedback-aa')
    #fbpath = fbinterface.namespace('feedback-aa')
    fbpath = Naming.dir_root('feedback-aa')

    # overwrite the paths for testing
    #workspace = '/gpfs/alpine/proj-shared/lrn005/Campaign3/roots/feedback_speedup_tests/feedback-aa2cg'
    #workspace = '/g/g92/aydin1/mummi/CODE_ROOT/mummi/feedback-aa2cg'
    workspace = '/g/g92/aydin1/mummi_new/feedback-aa2cg'
    outpath = workspace
    fbpath = workspace

    hvr_th = 0.15
    crd_th = 1.15
    nframes_fb_threshold = 3


    # file names to load
    #traj_name = os.path.join(inpath, 'CUT.xtc')
    #top_name  = os.path.join(inpath, 'CUT.gro')

    # itp files
    #main_itp = os.path.join(dpath, 'KRAS-GTP-04-HVR-best-guess-M22-CYFpos_main.itp')
    #revised_itp = os.path.join(inpath, 'Protein_temporary.itp')
    #revised_itp = os.path.join('Protein_temporary.itp')
    #revised_itp = os.path.join('Protein_reduced.itp')

    #itpkey = "KRAS-GTP-04-HVR-best-guess-M22-CYFpos_main.itp"
    #itpkey = "RAS_RAF_TERNARY_scfix-CYFpos.itp"

    #pdbkey_hvr = "protein_chainA.pdb"
    #pdbkey_crd = "protein_chainB.pdb"

    #protein_structure =os.path.join(dpath, 'KRAS-GTP-protein.pdb')


    # start with one itp file
    #multiple_itps = False

    load_checkpoint = False


    # --------------------------------------------------------------------------
    # this is a fake workflow
    fb = FeedbackManager_AA2CG(FeedbackManagerType.Manager,
                               workspace, outpath, fbpath,
                               iointerface, fbinterface, 
                               hvr_th, crd_th, nframes_fb_threshold)
                                       
    num = 0

    for it in range(0,1): #18):

        # copy the same and rename it

        #single_result = fb.load_analysis(top_name, traj_name)
        #overall_result = fb.aggregate_analysis(single_result)
        if it%10 == 0:
            fb.aggregate()

        if it%4 == 0:
            fb.report()

        #overall_result = fb.get_aggregate()
        #fb.do_feedback(overall_result)

        #if it==0:
        fb.checkpoint() #'checkpoint', overall_result)

        #if it%10 == 2:

        #    for i in range(0,5):
        #        os.system('cp' + ' ' + outpath_new + '/pfpatch_01_cgf01_aaf00000000_analysis.npz' + ' ' + outpath + '/pfpatch_01_cgf01_aaf00000000{0}_analysis.npz'.format(num))
        #        num = num + 1

        #sleep(100)


    if load_checkpoint == True:
        print("loading checkpoint now")
        fb.restore() #'checkpoint_2020-08-05_12-20-44')



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
