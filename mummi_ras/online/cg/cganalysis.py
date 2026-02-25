#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
#
# @authors  Harsh Bhatia  <bhatia4@llnl.gov>
#           Helgi I. Ingolfsson <ingolfsson1@llnl.gov>
#           Tim Carpenter <carpenter36@llnl.gov>
#           Christopher A. Neale <cneale@lanl.gov>
#           Francesco Di Natale <dinatale3@llnl.gov>
#           Jeremy O. B. Tempkin <tempkin1@llnl.gov>
#           Loic Pottier <pottier1@llnl.gov>
# -----------------------------------------------------------------------------

import argparse
import multiprocessing as mp
import numpy as np
import os
import time
import traceback
from MDAnalysis.analysis import distances 

import mummi_core
from mummi_core.utils.logger import init_logger
from mummi_core.utils import Naming
from mummi_core.utils.utilities import read_arg, sys_call
from mummi_core.workflow.feedback_manager import FeedbackManagerType

import mummi_ras
from mummi_ras.online import Simulation, ddcMD, GROMACS, GROMACS_PARTS, MDEngine

import mummi_ras.online.cg.cnlibs as cnlibs
import mummi_ras.online.cg.fast_get_rdf
import mummi_ras.online.cg.get_tiltrot_z_state_multiple

from mummi_ras.online.cg.cgoptions import CgAnalyzerOptions
from mummi_ras.datastructures.rdf import RDF
from mummi_ras.parsers.cg import CGParser

from mummi_ras.ml.selectors import CGSelector, CGSelectorType
from mummi_ras.feedback.manager.feedback_manager_cg2ml import FeedbackManager_cg2ml

from typing import List

sim = None # global var for ddcMD sim

from logging import getLogger
LOGGER = getLogger(__name__)

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
## Hide:
### cganalysis.py:599: DeprecationWarning: Conversion of an array with ndim > 0 to a scalar is deprecated, and will error in future. Ensure you extract a single element from your array before performing this operation. (Deprecated NumPy 1.25.)

def processOutput(options : CgAnalyzerOptions, cOutputList, iointerface):
    """Parse the cOutoutList dictionary and send feedback messages and write to analysis files."""

    # Sort output list by frame number - so they are processed in order
    cOutputList.sort(key=lambda x: x["frameNum"])

    # loop over frames 
    # this loops only processes data for the feedback mechanism. Rest of cOutput gets written later
    for cOutput in cOutputList:

        if len(cOutput) <= 1:
            LOGGER.warning(f"Frame {cOutput['frameNum']} invalid, skipping for feedback")
            continue

        # Only do any feedback if some protein in patch (i.e. lipid only patches are currently not used for feedback)
        # if len(cOutput["states"]) > 0:
        # In C4 Only do any feedback if one protein in the patch - not for lipids only and not for dimers

        # For C4 CG to AA selection is done for all monomers (of all x4 protein types)
        # turning off for ADMIRRAL to atomistic NOTE: Can only work on frames saved in positions.tar
        #if (len(cOutput["states"]) == 1) and ((cOutput["frameNum"] % options.counterSaveIncr) == 0):
        #state = cOutput["states"][0][0][5]
        #tilt = cOutput["states"][0][0][2]
        #rot = cOutput["states"][0][0][3]
        # @WARNING only do if not RAS-RAF  - check
        #depth = cOutput["states"][0][0][4]
        # LOGGER.info(f"cOutput = {cOutput.keys()}")

        #if depth == 'n/a':
        #    depth = np.nan
        LOGGER.info(f'Add feedback for frameNum {cOutput["frameNum"]} - state {cOutput["states"]}')
        # we are removing this! as no atomistic in PC. Normally put in ROOT/ml/cg
        #options.cgs.add_snapshot(options.simname, cOutput["frameNum"], tilt, rot, depth)
        if options.fb != None:
            options.fb.aggregate(options.simname, cOutput["frameNum"], cOutput["prot_positions"], cOutput["states"])

    #LOGGER.info("out of the loop")
    ## this is outside the "for cOutput in cOutputList" loop
    #if len(cOutputList[0]["states"]) > 0:
    # In C4 Only do any feedback if one protein in the patch - not for lipids only and not for dimers
    # LOGGER.info('report now')
    if options.fb != None:
        options.fb.report()
    # LOGGER.info('reported')

    # Update how many frames we have proccessed
    options.fcounter = cOutputList[-1]["frameNum"]
    # @TODO check for missing frames - in threory a frame before this one
    # could still be in processing.

    LOGGER.info('processOutput: {}'.format(options.fcounter))

    # get a list of all frameNums being processed
    listAll = list([x["frameNum"] for x in cOutputList[:]])

    cSavePath = options.path
    if options.pathremote != "":
        # Copy other done files to remote (e.g. restart)
        copyRemoteRestart(options, listAll)
        cSavePath = options.pathremote

    # Consolidate and save all analysis
    iointerface.save_analysis_data(cOutputList, cSavePath)

    # Consolidate and save all processed snapshots
    saveList = []
    for cFramNum in listAll:
        if (cFramNum % options.counterSaveIncr == 0):
            saveList.append(cFramNum)
    #iointerface.save_snapshots(saveList, options.path, cSavePath)

    # Copy over all logs - @WARNING these can be ok big and very vastfull to to this here - can remove if needed.
    if options.pathremote != "":
        copyRemoteLogs(options)

    # @Warning snapshots are always cleared - restart only if copied above
    # (using pathremote)
    cleanupLocal(options, listAll)

def run_main_analysis(options: CgAnalyzerOptions, cg: CGParser, iointerface, lipid_types: dict, simulation: Simulation, mini_mummi: bool):
    cfcounter = options.fcounter
    # only used for GROMACS as frames only read/loaded in first round
    reading = True
    cFrameList = []
    cOutputList = []    # List of frames that have alreaddy been analyzed
    cProcList = np.empty((options.nt), dtype=object)
    outputQueue = mp.Queue()

    startup_time = 5
    wait_time_s = 2

    # keep looking for new frames
    # check that the simulation is stil running
    restarts_left = 0
    if not options.noninteractive:
        restarts_left = simulation.restart_limit

    # in noninteractive set frame count to the full trajectory
    if options.noninteractive:
        cFrameList = list(range(0, cg.syst.trajectory.n_frames))

    while True:
        if not options.noninteractive:
            # Check if the simulation is successful.
            if simulation.successful and simulation.running:
                LOGGER.info(f"Simulation successful")
                # We should stop the simulation
                simulation.stop()
            elif not simulation.running and not simulation.successful:
                # Check if the simulation is still running.
                if restarts_left > 0:
                    simulation.setup_restart()
                    restarts_left -= 1
                    simulation.run()
                    LOGGER.info(f"Started simulation, waiting {startup_time} secondes...")
                    time.sleep(startup_time)

            # -----------------------------------------------------
            # Read in all new files (currently ready in file system)
            if options.MDEngine == MDEngine.GROMACS_PARTS:
                # Only add new frames if current list has been processed (as old frame would be in a diffrent .xtc)
                if cFrameList == []:
                    cFrameList = simulation.get_new_frames()
            else:
                cFrameList = cFrameList + simulation.get_new_frames()

        # -----------------------------------------------------
        # Process a single frame
        if len(cFrameList) > 0:
            if options.nt == 1:
                cFrameNum = cFrameList.pop(0)
                LOGGER.info(f"Process frame {cFrameNum}")
                cOut = processFrame(cFrameNum, cg, options, lipid_types, mini_mummi)
                outputQueue.put(cOut)
            elif options.nt > 1:
                # Find empty or finished precessor
                for i in range(0, len(cProcList)):
                    if cProcList[i] is None or not cProcList[i].is_alive():
                        cFrameNum = cFrameList.pop(0)
                        cProcList[i] = mp.Process(target=threadProcessFrame, args=(cFrameNum, cg, options, outputQueue, lipid_types, mini_mummi))
                        LOGGER.info(f"Process frame {cFrameNum} in subprocess num / name {i} / {cProcList[i].name}")
                        cProcList[i].start()
                        # break for loop; only start one frame per outer loop
                        break

        # -----------------------------------------------------
        # Get output from already processed frames - as soon as it's ready
        # @WARNING output that is not retrieved results in waiting processes
        # (they dont' terminate unless you emtpy the queue)
        while not outputQueue.empty():
            cOutputList.append(outputQueue.get())

        if len(cOutputList) >= options.frameProcessBatchSize:
            processOutput(options, cOutputList, iointerface)
            # This will empty the cOutputList
            # cOutputList.clear() # Empty the list after processing -
            # change back to this for python3
            del cOutputList[:]  # Empty the list after processing

        if not options.noninteractive and simulation.running:
            # -----------------------------------------------------
            # Wait before next run
            # @NOTE Adjust depending on need - this is a busy loop
            LOGGER.debug(f"Sleeping for {wait_time_s} seconds")
            time.sleep(wait_time_s)

        # Check simulation status
        success_flag, _, stop_flag = Naming.status_flags("cg")

        # -----------------------------------------------------
        # Stop criteria is max number of frames
        # If reached, stop sim, process last frames and exit
        if not options.noninteractive and options.fcounter >= options.counterSimStop:
            LOGGER.info("Sim {} is past {} steps, will now be stopped.".format(options.simname, options.counterSimStop))
            # Stop simulation run
            simulation.stop()
            # Copy everyhing to remote
            if len(cOutputList) > 0:
                processOutput(options, cOutputList, iointerface)

            # Mark simulations as converged - will be read by workflow - and
            # copy other logs to remote

            if options.pathremote != "":
                iointerface.send_signal(options.pathremote, success_flag)
                copyRemoteLogs(options)
            else:
                iointerface.send_signal(options.path, success_flag)

            # We are now done
            LOGGER.info("EXITING: %s", options.simname)
            exit(0)
        elif not options.noninteractive and options.fcounter >= options.counterSimStopMinLimit:
            if options.pathremote != "":
                _stop = iointerface.test_signal(options.pathremote, stop_flag)
            else:
                _stop = iointerface.test_signal(options.path, stop_flag)
            if _stop:
                LOGGER.info(
                    "Sim {} is past {} minimum required steps and it has been found uninteresting, will now be stopped.".format(options.simname, options.counterSimStopMinLimit))

                # Stop simulation run 
                simulation.stop()

                # Copy everything to remote storage
                if len(cOutputList) > 0:
                    processOutput(options, cOutputList, iointerface)

                if options.pathremote != "":
                    iointerface.send_signal(options.pathremote, success_flag)
                    copyRemoteLogs(options)
                else:
                    iointerface.send_signal(options.path, success_flag)

                # We are now done
                LOGGER.info("EXITING: %s", options.simname)
                exit(0)

        # -----------------------------------------------------
        # Stop the function here if simply testing the analysis routines
        if options.test:
            return

        # If list to process is empty and not in interactive mode - exit while loop and wait for last snapshots.
        if len(cFrameList) <= 0 and (options.noninteractive or not simulation.running):
            LOGGER.info(f"Simulation is not running and no frames to process, stopping")
            break

    # -------------------------------------------
    # Exited main while - now finalize wait for sub-process to finish and
    # process last data
    LOGGER.info(f"CG analysis is stopping for {options.simname}.")

    if options.nt > 1:
        while True:
            while not outputQueue.empty():
                cOutputList.append(outputQueue.get())
            cAive = 0
            for p in cProcList:
                if p is not None and cProcList[i].is_alive():
                    cAive += 1
            if cAive == 0:
                break
            LOGGER.info(
                "Waiting for {} subprocess there are still alive.".format(cAive))
            time.sleep(2)

    if len(cOutputList) > 0:
        LOGGER.info("Process ouput from last frames")
        processOutput(options, cOutputList, iointerface)

    # Copy other logs to remove
    if options.pathremote != "":
        copyRemoteLogs(options)

    outputQueue.close()
    return

def copyRemoteLogs(options: CgAnalyzerOptions):
    """ PLACEHOLDER """

    # Copy over all log files - warning these can be big so don't do all the
    # time. @TDOO we want to do this in a smarter way !!!
    try:
        LOGGER.info("Copy over log files")
        copyDirList = ""
        if options.MDEngine == MDEngine.DDCMD:
            copyDirList = copyDirList + " " + os.path.join(options.path, "data")
            copyDirList = copyDirList + " " + os.path.join(options.path, "ddcMD.header")
            copyDirList = copyDirList + " " + os.path.join(options.path, "hpm.data")
            copyDirList = copyDirList + " " + os.path.join(options.path, "ddcmd.out")
        else:
            copyDirList = copyDirList + " " + os.path.join(options.path, "cg_analysis.out")
        # @TODO don't copy the python log if it gets big
        copyDirList = copyDirList + " " + os.path.join(options.path, "*.log")
        sys_call("cp -r "+copyDirList+" "+options.pathremote+"/ || true")
        LOGGER.info(f"copyRemoteLogs = {copyDirList}")
    except Exception as e:
        LOGGER.error(f"cganalysis - copy to remote error {e}")
    return 

def cleanupLocal(options: CgAnalyzerOptions, frameList: List):
    """Cleanup local files"""

    # @WARNING currently this only works for ddcMD file structure
    if options.MDEngine != MDEngine.DDCMD:
        LOGGER.warning("Copy cleanup local only works for ddcMD files")
        return

    for cFrameIndex in frameList:
        if (cFrameIndex % options.counterCheckpointIncr == 0):
            #and options.pathremote == "":

            # Keep dir but remove snapshot file
            cFrameName = options.get_sim_file_name(cFrameIndex, options.fileformat)
            sys_call("rm "+cFrameName)
            LOGGER.debug('Remove snapshot {} num {}'.format(cFrameName, cFrameIndex))
            # If folder is empty - restat has been previusly removed - remove folder
            cFolderName = os.path.join(options.path, options.get_frame_folder_name(cFrameIndex))
            if len(os.listdir(cFolderName)) == 0:
                sys_call("rm -r "+cFolderName)
                LOGGER.debug('Empty remove frame {} num {}'.format(cFrameName, cFrameIndex))
        else:
            # Remove folder
            cFrameName = os.path.join(
                options.path, options.get_frame_folder_name(cFrameIndex))
            sys_call("rm -r " + cFrameName)
            LOGGER.debug("Remove frame {} num {}".format(cFrameName, cFrameIndex))

    return 


def copyRemoteRestart(options: CgAnalyzerOptions, frameList: List):
    """Copy over ddcMD restart files + any analysis (in folders)"""

    if options.MDEngine == MDEngine.DDCMD:
        try:
            LOGGER.info("Copy to remote frames {}".format(frameList))
            copyDirList = ""
            newCheckpointIndex = -1

            for copyFrameIndex in frameList:
                # Save snapshot for storage
                #if (copyFrameIndex % options.counterSaveIncr == 0):
                #    copyDirList = copyDirList + " " + os.path.join(options.path, options.get_frame_folder_name(copyFrameIndex))
                # Not used now as the files are saved/aggregated in processOutput

                # We have a candidat checkpoint to use as a restart file
                if (copyFrameIndex % options.counterCheckpointIncr == 0):
                    newCheckpointIndex = copyFrameIndex
                    copyDirList = copyDirList + " " + os.path.join(options.path, options.get_frame_folder_name(copyFrameIndex))

            # If no checkpoint file found there is nothing todo ``
            if copyDirList == "":
                return

            # Dont copy restart link- see below # Copy over ddcMD restart link;
            # ddcMD link to nevest restart file
            # copyDirList = copyDirList + " " + \
            #   os.path.join(options.path, "restart")

            # Copy over python log (for this script) - Name hard coded below
            # @WARNING cp of log is not nessisary can be left for end
            # copyDirListAll = copyDirList + " " + \
            #    os.path.join(options.path, "analysis_logger.log")
            # sys_call('cp -r '+copyDirList+' '+options.pathremote+'/')
            # @TODO change to to for loops
            sys_call("rsync -va --no-g " + copyDirList + " " + options.pathremote + "/ --exclude subset#000000*")
            # skip setup files and their .pdb's

            # After copy remove local data files - so not to fillup local dir
            #sys_call("rm -r "+copyDirList)

            # After copy if needed create a ddcMD restart link in remote dir
            # (@NODE the link ddcMD creates points to dir that are not copied yet)
            if (newCheckpointIndex > options.currentCheckpointIndex):
                sys_call("ln -sf ./{}/restart {}/restart".format(options.get_frame_folder_name(newCheckpointIndex), options.pathremote))
                options.currentCheckpointIndex = newCheckpointIndex

        except Exception as e:
            LOGGER.error(f"cganalysis - copy to remote error {e}")
    elif options.MDEngine == MDEngine.GROMACS_PARTS:
        try:
            LOGGER.info("Copy to remote frames {}".format(frameList))
            # if traj_comp.xtc is not there GROMACS will complain (even if we do not write to traj_comp.xtc)
            copyDirList = "traj_comp.xtc topol.tpr state.cpt state_prev.cpt ener.edr md.log mdrun.log"
            newCheckpointIndex = -1

            for copyFrameIndex in frameList:
                # We have a candidat checkpoint to use as a restart file
                # We have to convert the frame number into a timestep as the XTC filenames are in timesteps like traj_comp.step000002750000.xtc for ts = 2750000
                if (int(copyFrameIndex * options.counterIncr / options.bundlingsize) % options.counterCheckpointIncr == 0):
                    newCheckpointIndex = copyFrameIndex
                    copyDirList = copyDirList + " " + os.path.join(options.path, options.get_sim_analysis_name(copyFrameIndex))

            # If no checkpoint file found there is nothing todo
            if copyDirList == "":
                return

            sys_call(f"rsync -va --no-g {copyDirList} {options.pathremote} || true")
            # skip setup files and their .pdb's

            # We create a symlink for the last xtc analyzed. We can restart from that XTC
            if (newCheckpointIndex > options.currentCheckpointIndex):
                sys_call(f"ln -sf ./{options.get_sim_analysis_name(newCheckpointIndex)} {options.pathremote}/restart || true")
                options.currentCheckpointIndex = newCheckpointIndex

        except Exception as e:
            LOGGER.error(f"cganalysis - copy to remote error {e}")
    return

def get_protein_contacts_from_z(ppos, atomids_dict, cutoff, farnesyl_id):
    """Return an array with protein counts"""

    domains = atomids_dict.keys()

    ncontacts = {}
    ref_z = ppos[farnesyl_id,2]  # get farnesly height 

    for dom in domains:
        bead_zs = ppos[atomids_dict[dom].indices,2]
        diff = np.abs(bead_zs - ref_z)
        bary = diff <= cutoff
        ncontacts[dom] = bary.sum()

    return ncontacts

def get_state_from_protein_coordinates(ppos,
                                       atomids_dict,
                                       codes_to_state,
                                       cutoff,
                                       contacts_threshold,
                                       farnesyl_id,
                                       COM_1,
                                       COM_2):
    """
    Return state code for this protein structure.

    checks if we are in state A or B

    cutoff is the z distance to the CYF BB required and contacts_threshold is the
    number of beads within cutoff required for the domain to be "touching"
    """
    M1_domains = ["KD1", "CRD1", "loop11", "loop21", "b1433", "RAS1"]
    M2_domains = ["KD2", "CRD2", "loop12", "loop22", "RAS2"]

    _M1_dict = {key : atomids_dict[key] for key in M1_domains}
    _M2_dict = {key : atomids_dict[key] for key in M2_domains}

    # do M1 contacts
    contacts = get_protein_contacts_from_z(ppos, _M1_dict, cutoff, farnesyl_id[0])

    # make an array of all contacts
    code1 = ""
    # now accumulate all the codes into the dictionaries  
    for dom in M1_domains:

        # we are always assigning RAS to 1
        if dom == "RAS1":
            code1 += "1"
            continue

        if contacts[dom] >= contacts_threshold:
            code1 += "1"
        else:
            code1 += "0"

    # tag this with KD orientation 
    kdpos = [atomids_dict['KD1'].positions[:,0].mean(),
             atomids_dict['KD1'].positions[:,1].mean()]
    crdpos = [atomids_dict['CRD1'].positions[:,0].mean(),
              atomids_dict['CRD1'].positions[:,1].mean()]

    kdvec = np.array([kdpos[0] - COM_1[0], kdpos[1] - COM_1[1]])
    crdang = np.arctan2(crdpos[0] - COM_1[0], crdpos[1] - COM_1[1])

    kdrot = np.array([kdvec[0] * np.cos(-crdang) - kdvec[1] * np.sin(-crdang),
                      kdvec[0] * np.sin(-crdang) + kdvec[1] * np.cos(-crdang)]) 

    ymean = np.arctan2(kdrot[0], kdrot[1])

    if ymean >= 0:
        ytag = "ypos"
    else:
        ytag = "yneg"

    if code1 in codes_to_state.keys():
        state = code1
    else:
        state = "-".join([code1,ytag])

    M1_state = codes_to_state[state]

    # do M2 contacts
    contacts = get_protein_contacts_from_z(ppos, _M2_dict, cutoff, farnesyl_id[1])

    code2 = ""
    # now accumulate all the codes into the dictionaries  
    for dom in M2_domains:

        # we are always assigning RAS to 1
        if dom == "RAS2":
            code2 += "1"
            continue

        if contacts[dom] >= contacts_threshold:
            code2 += "1"
        else:
            code2 += "0"

    # tag this with KD orientation 
    kdpos = [atomids_dict['KD2'].positions[:,0].mean(),
             atomids_dict['KD2'].positions[:,1].mean()]
    crdpos = [atomids_dict['CRD2'].positions[:,0].mean(),
              atomids_dict['CRD2'].positions[:,1].mean()]

    kdvec = np.array([kdpos[0] - COM_2[0], kdpos[1] - COM_2[1]])
    crdang = np.arctan2(crdpos[0] - COM_2[0], crdpos[1] - COM_2[1])

    kdrot = np.array([kdvec[0] * np.cos(-crdang) - kdvec[1] * np.sin(-crdang),
                      kdvec[0] * np.sin(-crdang) + kdvec[1] * np.cos(-crdang)]) 

    ymean = np.arctan2(kdrot[0], kdrot[1])

    if ymean >= 0:
        ytag = "ypos"
    else:
        ytag = "yneg"

    if code2 in codes_to_state.keys():
        state = code2
    else:
        state = "-".join([code2,ytag])

    M2_state = codes_to_state[state]

    return ("_".join(["state",M1_state]), "_".join(["state", M2_state]))


def processFrame(cFrameNum: int, cg: CGParser, options: CgAnalyzerOptions, lipid_types: dict, mini_mummi = False):
    """
    Run all analysis returns a dictionary of all results.
    It returns a dictionary contaning the following:
    """

    lastTime = time.time()

    # add frame number
    return_dic = {"frameNum": cFrameNum}

    if options.MDEngine == MDEngine.DDCMD:
        cFrameName = options.get_sim_file_name(cFrameNum, options.fileformat)
        # We need to convert ddcMD frame to pdb and use the pdb name
        cFrameName = cFrameName + '.pdb'
    else:
        # used in GROMACS and in noninteractive mode also uses this variable 
        cFrameName = cg.gmx_trajfile

    if not cg.isSetup():
        # This will not happen if initlized first
        # cg not initilzed yet, then this is the first time called and we
        # have to create
        cg.setup(cFrameName, options.fnameTopo)
        # Set all cg MDAnalysis selections - that are used by all frames
        options.set_cgSelections(cg)
    else:
        LOGGER.info(f"Load {len(cg.syst.trajectory)} frames: current frame ID {cg.syst.trajectory.filename}")
        if options.MDEngine == MDEngine.GROMACS_PARTS:
            # For simulations split in parts - we need the subindex for that part
            # Note cg.fnametraj_ts should be the frame step of the current xtc file - updated in the simulation
            cPartFRameNum = int((options.counterIncr * cFrameNum - cg.fnametraj_ts) / options.counterIncr)
            LOGGER.debug(f"Convert cFrameNum {cFrameNum} to local frame {cPartFRameNum} from xtc file number {cg.fnametraj_ts}")
            cg.syst.trajectory[cPartFRameNum]
        elif len(cg.syst.trajectory) > cFrameNum:
            # Trajectory alreaddy loaded, set MDAnalysis to correct frame
            LOGGER.debug(f"Trajectory already loaded, setting MDA universe to Frame {cFrameNum}")
            # This modify the MDA universe to set it to cFrameNum
            cg.syst.trajectory[cFrameNum]
        else:
            # Normal for a ddcMD run load each rame into CGTrajt
            # cg parser shold alreaddy be initilized - we need to load in new
            # frame and update leaflets

            # @TODO should not happen with GROMACS or GROMACS_PARTS need to update for ddcMD 
            LOGGER.error("ERROR accessing frame not in loaded trajectory - should not happen")
            
            ## @todo ddcMDFrameName is not defined? what are we suposed to do here?
            # cg.update_frame(ddcMDFrameName)

    # get box size and lipid layer selections
    bbox = [0, cg.syst.dimensions[0], 0, cg.syst.dimensions[1]]
    tp, ce, bt = cg.get_layers()
    LOGGER.info(
        "Frame: {0:5d}, Time: {1:8.3f} ps, Bbox = {2:s}"
        .format(cFrameNum, cg.syst.trajectory.time, str(bbox)))

    LOGGER.info(f"Lengths of layers: {len(tp)} {len(ce)} {len(bt)}")
   
    # aligner should center the protein so that we
    # 1) avoid PBC issues  
    # 2) prep it for ML processing
    cg.align_frame()

    # get all the protein positions
    if mini_mummi:
        prot_positions = cg.get_selection("index 0 to 762")  # revised selection for mini mummi
    else:
        prot_positions = cg.get_selection("index 0 to 4715")  # dimer all protein indexes (including ATP, ZN, etc..)
        
    
    return_dic["prot_positions"] = prot_positions.positions 

    # FA: Do we want to keep native fractions analysis? Jeremy may need to check this part to see if it's still applicable here.   
    if not mini_mummi:
        # compute the contact distances
        native_fractions = []
        for int_id, interface in enumerate(options.get_interface_selections()):
            seli, selj = interface

            # compute the contacts with MDAnalysis call
            _all_distances = distances.distance_array(cg.syst.atoms[seli],
                                                    cg.syst.atoms[selj],
                                                    box=cg.syst.dimensions)
    
            # get the atom selections for each interface pair
            return_dic[f"contact_distances_interface_{int_id+1}"] = _all_distances
            _ncontacts = (_all_distances < 11.0).sum()  #  
            _nf = (_ncontacts, _ncontacts / options.params['interfaces'][f'interface_{int_id+1}']['num_native'])
            native_fractions.append(_nf)
        return_dic["native_fractions"] = native_fractions

    results = cg.syst.dimensions[0:3]
    return_dic["box"] = results

    # include dummy states array, placeholding for further implementation.
    #return_dic["states"] = np.array([[[-1,-1,-1,"n/a","n/a","n/a"]]])

    # additional distance analysis
    if mini_mummi:
        KRAS_ASP38 = cg.get_selection('name BB and resname ASP and index 81')
        KRAS_GLN25 = cg.get_selection('name BB and resname GLN and index 48')
        BRAF_ARG166 = cg.get_selection('name BB and resname ARG and index 476')
        BRAF_MET187 = cg.get_selection('name BB and resname MET and index 522')
    else:
        m1_KRAS_ASP38 = cg.get_selection('name BB and resname ASP and index 81')
        m1_KRAS_GLN25 = cg.get_selection('name BB and resname GLN and index 48')
        m1_BRAF_ARG166 = cg.get_selection('name BB and resname ARG and index 476')
        m1_BRAF_MET187 = cg.get_selection('name BB and resname MET and index 522')

        m2_KRAS_ASP38 = cg.get_selection('name BB and resname ASP and index 2978')
        m2_KRAS_GLN25 = cg.get_selection('name BB and resname GLN and index 2945')
        m2_BRAF_ARG166 = cg.get_selection('name BB and resname ARG and index 3373')
        m2_BRAF_MET187 = cg.get_selection('name BB and resname MET and index 3419')

    # commenting out for mini mummi
    if not mini_mummi:
        m1_SEP365 = cg.get_selection('resname SEP and index 942 and name SC1')
        m1_SEP729 = cg.get_selection('resname SEP and index 1783 and name SC1')
        #m1_SEP365_ARG56 = cg.get_selection('resname ARG and index 1943 and name SC2')# and resid 56 
        #m1_SEP365_ARG127 = cg.get_selection('resname ARG and index 2116 and name SC2')#and resid 127 
        m2_SEP729_ARG56 = cg.get_selection('resname ARG and index 2482 and name SC2')# and resid 56 m2 14-3-3
        m2_SEP729_ARG127 = cg.get_selection('resname ARG and index 2655 and name SC2')#and resid 127 m2  14-3-3
        m1_ZN740 = cg.get_selection('name ZN and index 1805')#and resid 740
        m1_CYM280 = cg.get_selection('resname CYM and name SC1 and index 750')#and resid 280 FA:changed CYS to CYM 

        m2_SEP365 = cg.get_selection('resname SEP and index 3839 and name SC1')
        m2_SEP729 = cg.get_selection('resname SEP and index 4680 and name SC1')
        #m2_SEP365_ARG56 = cg.get_selection('resname ARG and index 2482 and name SC2')
        #m2_SEP365_ARG127 = cg.get_selection('resname ARG and index 2655 and name SC2')
        m1_SEP729_ARG56 = cg.get_selection('resname ARG and index 1943 and name SC2') # m1 14-3-3
        m1_SEP729_ARG127 = cg.get_selection('resname ARG and index 2116 and name SC2') # m1 14-3-3
        m2_ZN740 = cg.get_selection('name ZN and index 4702')#and resid 740
        m2_CYM280 = cg.get_selection('resname CYM and name SC1 and index 3647')#and resid 280 FA:changed CYS to CYM

        # Residues at the KD-KD interface
        m1_KD_ASP565 = cg.get_selection('name BB and index 1389')
        m1_KD_PHE516 = cg.get_selection('name BB and index 1270')
        m1_KD_GLN530 = cg.get_selection('name BB and index 1303')

        m2_KD_HSE477 = cg.get_selection('name BB and index 4080')
        m2_KD_ARG509 = cg.get_selection('name BB and index 4150')
        m2_KD_HSD510 = cg.get_selection('name BB and index 4153')

        # Residues at 1433 - KD interface
        m1_1433_ARG127 = cg.get_selection('name BB and index 2114')
        m2_KD_SEP729 = cg.get_selection('name BB and index 4679')


    if mini_mummi:
        RAS = cg.get_selection('name BB and index 0 to 437')
        RBD = cg.get_selection('name BB and index 452 to 618')
        CRD = cg.get_selection('name BB and index 637 to 740')
        CRD_loops = cg.get_selection('name BB and index 649 to 660 685 to 694')
    else:
        m1_RAS = cg.get_selection('name BB and index 0 to 437')
        m1_RBD = cg.get_selection('name BB and index 452 to 618')
        m1_CRD = cg.get_selection('name BB and index 637 to 740')
        m1_CRD_loops = cg.get_selection('name BB and index 649 to 660 685 to 694')

        m2_RAS = cg.get_selection('name BB and index 2897 to 3334')
        m2_RBD = cg.get_selection('name BB and index 3349 to 3515')
        m2_CRD = cg.get_selection('name BB and index 3534 to 3637')
        m2_CRD_loops = cg.get_selection('name BB and index 3546 to 3557 3582 to 3591')

    if mini_mummi:
        BRAF = cg.get_selection('name BB and index 452 to 762') # revised range for mini mummi
        domains = [RAS,RBD,CRD,CRD_loops,BRAF]
        domains_COM = [RAS.center_of_geometry(),RBD.center_of_geometry(),CRD.center_of_geometry(),CRD_loops.center_of_geometry(), BRAF.center_of_geometry()]
        domains_COM = np.array(domains_COM)
    else:
        m1_loop1 = cg.get_selection('name BB and index 753 to 928')
        m1_loop2 = cg.get_selection('name BB and index 955 to 1135')
        m1_KD = cg.get_selection('name BB and index 1136 to 1768')
        m1_1433 = cg.get_selection('name BB and index 1819 to 2357')
        m1_BRAF = cg.get_selection('name BB and index 452 to 1818') # FA: I removed the 14-3-3 from this, is it reasonable?

        m2_loop1 = cg.get_selection('name BB and index 3650 to 3825')
        m2_loop2 = cg.get_selection('name BB and index 3852 to 4032')
        m2_KD = cg.get_selection('name BB and index 4033 to 4665')
        m2_1433 = cg.get_selection('name BB and index 2358 to 2896')
        m2_BRAF = cg.get_selection('name BB and index 3349 to 4715') # No 14-3-3


        m1_domains = [m1_RAS,m1_RBD,m1_CRD,m1_CRD_loops,m1_loop1,m1_loop2,m1_KD,m1_1433,m1_BRAF]
        m1_domains_COM = [m1_RAS.center_of_geometry(),m1_RBD.center_of_geometry(),m1_CRD.center_of_geometry(),m1_CRD_loops.center_of_geometry(),m1_loop1.center_of_geometry(),
                    m1_loop2.center_of_geometry(),m1_KD.center_of_geometry(), m1_1433.center_of_geometry(), m1_BRAF.center_of_geometry()]

        m2_domains = [m2_RAS,m2_RBD,m2_CRD,m2_CRD_loops,m2_loop1,m2_loop2,m2_KD,m2_1433,m2_BRAF]
        m2_domains_COM = [m2_RAS.center_of_geometry(),m2_RBD.center_of_geometry(),m2_CRD.center_of_geometry(),m2_CRD_loops.center_of_geometry(),m2_loop1.center_of_geometry(),
                    m2_loop2.center_of_geometry(),m2_KD.center_of_geometry(), m2_1433.center_of_geometry(), m2_BRAF.center_of_geometry()]
    

    if mini_mummi:
        domains_COM = np.array(domains_COM)
        PO4 = cg.get_selection('name PO4')

        Domains_membrane_dist = np.zeros([1,len(domains)])
        
        RASRAF_distances = np.zeros([1,3])
    
        for d in range(len(domains_COM)):
            temp = PO4.atoms.center_of_geometry()[2]-domains_COM[d][2]
            if temp < 0:
                temp = temp + cg.syst.dimensions[2]
            Domains_membrane_dist[0,d] = temp
    else:
        m1_domains_COM = np.array(m1_domains_COM)
        m2_domains_COM = np.array(m2_domains_COM)
        PO4 = cg.get_selection('name PO4')

        m1_Domains_membrane_dist = np.zeros([1,len(m1_domains)])
        m2_Domains_membrane_dist = np.zeros([1,len(m2_domains)])

        m1_RASRAF_distances = np.zeros([1,3])
        m2_RASRAF_distances = np.zeros([1,3])

        for d in range(len(m1_domains_COM)):
            temp_m1 = PO4.atoms.center_of_geometry()[2]-m1_domains_COM[d][2]
            if temp_m1 < 0:
                temp_m1 = temp_m1 + cg.syst.dimensions[2]
            m1_Domains_membrane_dist[0,d] = temp_m1

        for d in range(len(m2_domains_COM)):
            temp_m2 = PO4.atoms.center_of_geometry()[2]-m2_domains_COM[d][2]
            if temp_m2 < 0:
                temp_m2 = temp_m2 + cg.syst.dimensions[2]
            m2_Domains_membrane_dist[0,d] = temp_m2

    # Not for mini MuMMI
    if not mini_mummi:
        xy_threshold = 40
        z_threshold = 80
        BB = cg.get_selection("name BB")
        if np.max(BB.atoms.positions[:,0])-np.min(BB.atoms.positions[:,0]) > (cg.syst.dimensions[0]-xy_threshold):
            x_pbc = True
        else:
            x_pbc = False
        if np.max(BB.atoms.positions[:,1])-np.min(BB.atoms.positions[:,1]) > (cg.syst.dimensions[1]-xy_threshold):
            y_pbc = True
        else:
            y_pbc = False
        if np.max(BB.atoms.positions[:,2])-np.min(BB.atoms.positions[:,2]) > (cg.syst.dimensions[2]-z_threshold):
            z_pbc = True
        else:
            z_pbc = False
        pbc_checks = [x_pbc, y_pbc, z_pbc]
    
    # DeprecationWarning: Conversion of an array with ndim > 0 to a scalar is deprecated, and will error in future.
    # Ensure you extract a single element from your array before performing this operation. (Deprecated NumPy 1.25.)
    if mini_mummi:
        RASRAF_distances[0,0] = distances.distance_array(KRAS_ASP38,BRAF_ARG166)
        RASRAF_distances[0,1] = distances.distance_array(KRAS_GLN25,BRAF_MET187)
        RASRAF_distances[0,2] = distances.distance_array(RAS.center_of_geometry(),RBD.center_of_geometry())
    else:
        m1_RASRAF_distances[0,0] = distances.distance_array(m1_KRAS_ASP38,m1_BRAF_ARG166)
        m1_RASRAF_distances[0,1] = distances.distance_array(m1_KRAS_GLN25,m1_BRAF_MET187)
        m1_RASRAF_distances[0,2] = distances.distance_array(m1_RAS.center_of_geometry(),m1_RBD.center_of_geometry())

        m2_RASRAF_distances[0,0] = distances.distance_array(m2_KRAS_ASP38,m2_BRAF_ARG166)
        m2_RASRAF_distances[0,1] = distances.distance_array(m2_KRAS_GLN25,m2_BRAF_MET187)
        m2_RASRAF_distances[0,2] = distances.distance_array(m2_RAS.center_of_geometry(),m2_RBD.center_of_geometry())

    if not mini_mummi:
        m1_SEP_distances = np.zeros([1,4])
        m1_SEP_distances[0,0] = distances.distance_array(m1_SEP365,m2_SEP729_ARG56)
        m1_SEP_distances[0,1] = distances.distance_array(m1_SEP365,m2_SEP729_ARG127)
        m1_SEP_distances[0,2] = distances.distance_array(m1_SEP729,m2_SEP729_ARG56)
        m1_SEP_distances[0,3] = distances.distance_array(m1_SEP729,m2_SEP729_ARG127)
        m1_ZN_distances = distances.distance_array(m1_ZN740,m1_CYM280)

        m2_SEP_distances = np.zeros([1,4])
        m2_SEP_distances[0,0] = distances.distance_array(m2_SEP365,m1_SEP729_ARG56)
        m2_SEP_distances[0,1] = distances.distance_array(m2_SEP365,m1_SEP729_ARG127)
        m2_SEP_distances[0,2] = distances.distance_array(m2_SEP729,m1_SEP729_ARG56)
        m2_SEP_distances[0,3] = distances.distance_array(m2_SEP729,m1_SEP729_ARG127)
        m2_ZN_distances = distances.distance_array(m2_ZN740,m2_CYM280)

        KD_KD_distances = np.zeros([1,3])
        KD_KD_distances[0,0] = distances.distance_array(m1_KD_ASP565,m2_KD_HSE477)
        KD_KD_distances[0,1] = distances.distance_array(m1_KD_PHE516,m2_KD_ARG509)
        KD_KD_distances[0,2] = distances.distance_array(m1_KD_GLN530,m2_KD_HSD510)
        
        KD_1433_distances = distances.distance_array(m1_1433_ARG127,m2_KD_SEP729)

        # FA: Do we still want to keep these checks? I'm not touching this part for now
        # below are the checks for the ensemble B status:

        '''
        sep_cutoff = 8 #XXXXXX
        if SEP_distances[0,0] > sep_cutoff or SEP_distances[0,1] > sep_cutoff:
            SEP365_unbound = True
        else:
            SEP365_unbound = False
        
        zn_cutoff = 8 #XXXXXX
        if ZN_distances < zn_cutoff:
            zn_cys_bound = True
        else:
            zn_cys_bound = False
        '''

        crd_loops_cutoff = 30 #need to check this...
        if m1_Domains_membrane_dist[0,3] < crd_loops_cutoff:
            m1_crd_loops_bound = True
        else:
            m1_crd_loops_bound = False
       
        if m2_Domains_membrane_dist[0,3] < crd_loops_cutoff:
            m2_crd_loops_bound = True
        else:
            m2_crd_loops_bound = False

        '''
        KD_membrane_cutoff = 68
        if Domains_membrane_dist[0,6] < KD_membrane_cutoff:
            KD_membrane_distance = True
        else:
            KD_membrane_distance = False 
        '''

        # b_checks = [SEP365_unbound, zn_cys_bound, crd_loops_bound, KD_membrane_distance]
        b_checks = [m1_crd_loops_bound, m2_crd_loops_bound]

        return_dic["m1_sep_distances"] = m1_SEP_distances
        return_dic["m1_zn_distances"] = m1_ZN_distances
        return_dic["m2_sep_distances"] = m2_SEP_distances
        return_dic["m2_zn_distances"] = m2_ZN_distances

        return_dic["b_checks"] = b_checks
        return_dic["pbc_checks"] = pbc_checks

        return_dic["kd_kd_distances"] = KD_KD_distances
        return_dic["kd_1433_distances"] = KD_1433_distances

    if mini_mummi:
        Domain_Domain_dist = distances.distance_array(domains_COM,domains_COM)
    
        return_dic["domains_membrane_dist"] = Domains_membrane_dist
        return_dic["domain_domain_dist"] = Domain_Domain_dist
        return_dic["rasraf_distances"] = RASRAF_distances
    else:
        m1_Domain_Domain_dist = distances.distance_array(m1_domains_COM,m1_domains_COM)
        m2_Domain_Domain_dist = distances.distance_array(m2_domains_COM,m2_domains_COM)

        return_dic["m1_domains_membrane_dist"] = m1_Domains_membrane_dist
        return_dic["m2_domains_membrane_dist"] = m2_Domains_membrane_dist
        return_dic["m1_domain_domain_dist"] = m1_Domain_Domain_dist
        return_dic["m2_domain_domain_dist"] = m2_Domain_Domain_dist
        return_dic["m1_rasraf_distances"] = m1_RASRAF_distances
        return_dic["m2_rasraf_distances"] = m2_RASRAF_distances

    # FA: I haven't changed this part as I'm not sure if this is still valid for ensemble C 
    if not mini_mummi:

        domains_dict = {"KD1": prot_positions.select_atoms('index 1136 to 1768'),
                        "CRD1": prot_positions.select_atoms('index 637 to 740'),
                        "loop11": prot_positions.select_atoms('index 753 to 928'),
                        "loop21": prot_positions.select_atoms('index 955 to 1135'),
                        "b1433": prot_positions.select_atoms('index 2358 to 2896'),
                        "RAS1": prot_positions.select_atoms('index 0 to 437'),
                        "KD2": prot_positions.select_atoms('index 4034 to 4664'),
                        "CRD2": prot_positions.select_atoms('index 3535 to 3638'),
                        "loop12": prot_positions.select_atoms('index 3651 to 3826'),
                        "loop22": prot_positions.select_atoms('index 3853 to 4033'),
                        "RAS2": prot_positions.select_atoms('index 2898 to 3335')}

        # this is for B state classification only
        codes_to_state = {"000001": "a",
                          "000011": "a",
                          "000101": "a",
                          "000111": "a",
                          "010001": "b",
                          "010011": "b",
                          "010101": "b",
                          "010111": "b",
                          "100001-ypos": "c",
                          "100011-ypos": "c",
                          "100101-ypos": "c",
                          "100111-ypos": "c",
                          "110001-ypos": "d",
                          "110011-ypos": "d",
                          "110101-ypos": "d",
                          "110111-ypos": "d",
                          "100001-yneg": "e",
                          "100011-yneg": "e",
                          "100101-yneg": "e",
                          "100111-yneg": "e",
                          "110001-yneg": "f",
                          "110011-yneg": "f",
                          "110101-yneg": "f",
                          "110111-yneg": "f",
                          "001001": "g",
                          "001011": "g",
                          "001101": "g",
                          "001111": "g",
                          "011001": "b",
                          "011011": "b",
                          "011101": "b",
                          "011111": "b",
                          "101001-ypos": "i",
                          "101011-ypos": "i",
                          "101101-ypos": "i",
                          "101111-ypos": "i",
                          "111001-ypos": "d",
                          "111011-ypos": "d",
                          "111101-ypos": "d",
                          "111111-ypos": "d",
                          "101001-yneg": "k",
                          "101011-yneg": "k",
                          "101101-yneg": "k",
                          "101111-yneg": "k",
                          "111001-yneg": "f",
                          "111011-yneg": "f",
                          "111101-yneg": "f",
                          "111111-yneg": "f",
                          "00001": "a_Bstar",
                          "00011": "a_Bstar",
                          "01001": "b_Bstar",
                          "01011": "b_Bstar",
                          "10001-ypos": "c_Bstar",
                          "10011-ypos": "c_Bstar",
                          "11001-ypos": "d_Bstar",
                          "11011-ypos": "d_Bstar",
                          "10001-yneg": "e_Bstar",
                          "10011-yneg": "e_Bstar",
                          "11001-yneg": "f_Bstar",
                          "11011-yneg": "f_Bstar",
                          "00101": "g_Bstar",
                          "00111": "g_Bstar",
                          "01101": "b_Bstar",
                          "01111": "b_Bstar",
                          "10101-ypos": "i_Bstar",
                          "10111-ypos": "i_Bstar",
                          "11101-ypos": "d_Bstar",
                          "11111-ypos": "d_Bstar",
                          "10101-yneg": "k_Bstar",
                          "10111-yneg": "k_Bstar",
                          "11101-yneg": "f_Bstar",
                          "11111-yneg": "f_Bstar"}

        com_m1_GDOMAIN = prot_positions.select_atoms('index 0 to 380 and name BB').center_of_mass() # should be KRAS G-domain_
        com_m2_GDOMAIN = prot_positions.select_atoms('index 2897 to 3277 and name BB').center_of_mass()

        # protein atom indexes go as M1 is 0 to 2897 inclusive, then M2 is 2898 to 4716
        states_M12 = get_state_from_protein_coordinates(prot_positions.atoms.positions, 
                                                        domains_dict,
                                                        codes_to_state,
                                                        15,
                                                        5,
                                                        [430,3327],
                                                        com_m1_GDOMAIN,  # TODO 
                                                        com_m2_GDOMAIN)  # farnesyl id is hardcoded as BB

        return_dic["states"] = states_M12
    else:
        # We create an empty states to not break the code in case of Mini MuMMI
        return_dic["states"] = []

    # Remove all-to-all dist as to big 11MB per frame 
    # adding in the all to all matrix 
    #all_bbs = cg.get_selection("name BB")
    #a2a_distance = distances.distance_array(all_bbs, all_bbs)
    #return_dic["all_to_all_bb_dist"] = a2a_distance

    LOGGER.info(
        "Time for analysis - frameNum %s, time %s s",
        cFrameNum, time.time() - lastTime)

    # Return all data as a dictonary
    return return_dic


def threadProcessFrame(cFrameNum: int, cg: CGParser, options: CgAnalyzerOptions, outputQueue: mp.Queue, lipid_types: dict, mini_mummi: bool = False):
    try:
        outputQueue.put(processFrame(cFrameNum, cg, options, lipid_types, mini_mummi))
        LOGGER.info("Frame done {}".format(str(cFrameNum)))
    except Exception as e:
        outputQueue.put({"frameNum": cFrameNum})
        LOGGER.error("cganalysis - frame {} - sup process error {}".format(str(cFrameNum), e))
    return


def setup_argparser():
    """Set up the program's argument parser."""
    parser = argparse.ArgumentParser(description='Analyze CG simulations individually.')

    parser.add_argument("--simname", required=True,
                        help="Name of the simulation")

    parser.add_argument("-f", "--fcount", type=int, default=1,
                       help="Frame number for inital frame (default: 1)")

    parser.add_argument("--noninteractive", action="store_true",
                       help="Set if running on a previusly run folder - not with live simulations")

    parser.add_argument("-p", "--path", required=True, type=str,
                       help="Path for input and output")

    parser.add_argument("-pr", "--pathremote", type=str, default="",
                       help="If set periotically copy files to remote WARNING this will also remove everything in path on exit")

    parser.add_argument("--fstype", type=str, default="mummi", choices=mummi_ras.get_interfaces(),
                        help="Chose how to write files/data (default: mummi)")

    parser.add_argument("--fbio", type=str, default="mummi", choices=mummi_ras.get_interfaces(),
                        help="Choose how to write files/data (default: mummi)")

    parser.add_argument("--nprocs", type=int, default=1,
                        help="Numer of processes to use for analysis (default: 1)")

    parser.add_argument("-d", "--loglevel", type=int, default=2,
                        help="Level of logging messages to be output:\n"
                             "5 - Critical\n"
                             "4 - Error\n"
                             "3 - Warning\n"
                             "2 - Info (Default)\n"
                             "1 - Debug")

    parser.add_argument("-l", "--logpath", type=str, default="./",
                        help="Path for log information.")

    parser.add_argument("-c", "--logstdout", action="store_true",
                        help="Log to stdout in addition to a file.")

    parser.add_argument("--test", action="store_true",
                        help="Analyze a single structure file and exit.")

    parser.add_argument("--logfile", type=str, default="cg_analysis",
                        help="Name of cg analysis logger")

    parser.add_argument("--maxsimtime", type=int, default=250000000,
                        help="Stop simulation and analysis after this time step is reached")

    parser.add_argument("--stopsimtime", type=int, default=600,
                        help="Stop simulation and analysis after this simulated time is reached (in nanosecond)")

    parser.add_argument("--minsimtime", type=int, default=100000000,
                        help="Don't stop simulation and analysis before this time step is reached")

    parser.add_argument("--frameProcessBatchSize", type=int, default=10,
                        help="Show many frames are processed together")

    # Simulation specific arguments
    parser.add_argument("--siminputs", type=str, required=True,
                        help="Path to collection of simulation inputs.")

    parser.add_argument("--simbin", type=str, required=True,
                        help="Path to simulation binary to use for simulation. Empty string is don't run simulation")

    parser.add_argument("--simruntime", type=str, required=False, default=-1,
                        help="Max time ot run simulation (mostly for GROMACS -maxh).")

    parser.add_argument("--simcores", type=int, required=False, default=1,
                        help="Number of cores for the simulation (mostly for GROMACS).")

    parser.add_argument("--nofeedback", action="store_true",
                        help="Turn off the feedback calculation.")

    parser.add_argument("--backend", "-b", type=MDEngine.from_string, choices=list(MDEngine), default="NOSIM",
                        help="MD backend to use (default NOSIM, meaning just analysis no simulations)")

    parser.add_argument("--mini", action='store_true',
                        help="If set, CGAnalysis will use mini-MuMMI settings")

    args = parser.parse_args()

    if args.noninteractive:
        LOGGER.info(f"Changed backend from {args.backend} to NOSIM because noninteractive flag is set")
        args.backend = MDEngine.NOSIM

    if args.backend != MDEngine.NOSIM:
        if args.simbin == '':
            raise ValueError(f"--simbin='' while the backend is {args.backend}. A valid path is required.")

    return args


def run_cganalysis(args: argparse.Namespace):
    # Show options
    LOGGER.info("Running CG analysis, options are:")
    LOGGER.info(args)
    LOGGER.info("local path = {}".format(args.path))
    if args.pathremote != '':
        LOGGER.info("remote path = {}".format(args.pathremote))
    else:
        LOGGER.info("No remote path")

    # Initilize ouput DBR/files etc
    iointerface = mummi_ras.get_io(args.fstype)
    if args.nofeedback:
        fbinterface = None
    else:
        fbinterface = mummi_ras.get_io(args.fbio)
    LOGGER.info(f"Using IO_Interface [{iointerface}] for data output")
    LOGGER.info(f"Using IO_Interface [{fbinterface}] for feedback output")

    # Read in all argumetns and store in options objectselt
    options = CgAnalyzerOptions(args, iointerface, fbinterface)

    # Setup FeadbackManager
    # This will fail because Harsh commened Feedback manager's import
    # on Mar 19, 2020
    # until the new one is written!
    cPath = options.path
    if options.pathremote != "":
        cPath = options.pathremote

    # Create inital MDAnalysis object - make up front index definitions
    # Check .tpr file (only one per simulation but it has to be there)
    if not os.path.isfile(options.fnameTopo):
        error = "Simulation .tpr file {} not found".format(options.fnameTopo)
        LOGGER.error(error)
        raise ValueError(error)

    if options.MDEngine == MDEngine.DDCMD:
        # Check init .gro file (only one per simulation but it has to be there)
        # fname = "lipids-water-eq6.gro"
        fname = os.path.join(options.path, "lipids-water-eq6.gro")
    elif options.MDEngine == MDEngine.GROMACS:
        fname = os.path.join(options.path, "lipids-water-eq6.gro")
    elif options.MDEngine == MDEngine.GROMACS_PARTS:
        fname = os.path.join(options.path, "lipids-water-eq6.gro")
    else:
        # If not ddcMD we assume this is a GROMACS trajectory - .xtc file
        # provied as the sim name
        # fname = options.simname
        fname = os.path.join(options.path, "traj_comp.xtc")

    if not os.path.isfile(fname):
        error = "Simulation file {} not found".format(fname)
        LOGGER.error(error)
        raise ValueError(error)

    # Load Inital MDAnalysis universe (each frame after this will add to this
    # universe)
    # @NOTE Here we are lucky the sub-processes copy their own python stack
    # (new copy of this object)
    #  so they are all independent. If we later change to threadeds - we need
    # to change to array of cg
    #  objects with lazy load (commented out)
    # @NOTE this is slow ~15 sec each CGPrase creation
    cg = CGParser(fname, options.fnameTopo, options.MDEngine, empty = False, gmx_trajfile = os.path.join(options.path, "traj_comp.xtc"))
    # Lazy load code not used - so we make this a lazy load (so each sub-p can
    # do it independanly)
    # self.cgParseList = np.empty((self.nt), dtype=object)
    # for i in range(len(self.cgParseList)):
    #    #self.cgParseList[i] = CGParser(fname, self.fnameTopo, 1, empty=False)
    #    self.cgParseList[i] = CGParser(fname, self.fnameTopo, 1, empty=True)
    #    # Not implemented ???
    #    #self.cgParseList[i] = copy.deepcopy(self.cgParseList[0])
    # IF added change below from single to array

    # Set all cg MDAnalysis selections - that are used by all frames - keep
    # options class
    if cg.isSetup():
        options.set_cgSelections(cg)

    # Hardcoding lipid types here for now. We will likely want to parse this out of an input.
    lipid_types = {
        "nInner": ["POPC", "PAPC", "POPE", "DLPE", "SSM", "PAPS", "SAP6", "CHOL"],
        "nOuter": ["POPX", "PAPC", "POPE", "DLPE", "SSM", "CHOL"],
    }

    # Simulation setup
    # NOTE: We're initializing here, as it looks like the run_main_analysis
    # function has the main while True loop in it.
    # TODO: We need to refactor to pull out the while true loop to an actual
    # main method. It's bad practice to bury a while true loop.
    global sim
    if not options.noninteractive:
        success_flag, fail_flag, _ = Naming.status_flags('cg')
        out_path = options.pathremote if options.pathremote else options.path
        _pass = iointerface.test_signal(out_path, success_flag)
        _fail = iointerface.test_signal(out_path, fail_flag)
        LOGGER.info("Testing flags: [Success='%s', Fail='%s']", success_flag, fail_flag)
        LOGGER.info("Checking for exit condition [%s] (pass=%s fail=%s)", out_path, _pass, _fail)
        if _pass or _fail:
            LOGGER.info("DETECTED EXIT CONDITION: Pass=%s Fail=%s", _pass, _fail)
            exit(0)

        sim = run_simulation(args, iointerface, cg)

    # Find frames and run analysis
    # @WARNING for ddcMD jobs this is an infinit loop and never returns.
    # @TODO add ddcMD process id as an argument and in run_main_analysis check if still running if not report, copy back and terminate.
    run_main_analysis(options, cg, iointerface, lipid_types, sim, mini_mummi = args.mini)

    return


def run_simulation(args: argparse.Namespace, iointerface, cg: CGParser):
    if args.backend == MDEngine.DDCMD:
        sim = ddcMD(
            sim_name = args.simname,
            input_dir = args.siminputs,
            out_dir = args.path,
            remote_dir = args.pathremote,
            bin_path = args.simbin,
            cg = cg
        )
    elif args.backend == MDEngine.GROMACS:
        sim = GROMACS(
            sim_name = args.simname,
            input_dir = args.siminputs,
            out_dir = args.path,
            remote_dir = args.pathremote,
            bin_path = args.simbin,
            cg = cg,
            iointerface = iointerface,
            ncores = args.simcores,
            stopsimtime = args.stopsimtime,
            fcounter = args.fcount,
            runtime_hours = args.simruntime, # By default we run GROMACS with runtime -1 (no limit) 
            r_limit = 2,
            xtc_file = cg.gmx_trajfile #TODO: change that to not duplicate GMX input
        )
    elif args.backend == MDEngine.GROMACS_PARTS:
        sim = GROMACS_PARTS( 
            sim_name = args.simname,
            input_dir = args.siminputs,
            out_dir = args.path,
            remote_dir = args.pathremote,
            bin_path = args.simbin,
            cg = cg,
            iointerface = iointerface,
            ncores = args.simcores,
            stopsimtime = args.stopsimtime,
            fcounter = args.fcount,
            runtime_hours = args.simruntime, # By default we run GROMACS with runtime -1 (no limit) 
            r_limit = 2,
            xtc_file = cg.gmx_trajfile #TODO: change that to not duplicate GMX input
        )
    else:
        return None
    sim.initialize()
    sim.run()
    return sim

def main():
    # Set up the necessary base data structures to begin study set up.
    args = setup_argparser()

    # Sets up MuMMI and the logger
    mummi_core.init()
    mummi_core.create_root()
    mummi_core.init_logger(argparser = args)


    iointerface = mummi_ras.get_io(args.fstype)
    success_flag, fail_flag, stop_flag = Naming.status_flags('cg')

    try:
        run_cganalysis(args)
    except Exception as e:
        if "Uneven leaflets" in e:
            LOGGER.warning(f"Detected Uneven leaflets, setting flag")
            iointerface.send_signal(args.path, fail_flag)
        LOGGER.error(e)
        traceback.print_exc()
        raise e
    finally:
        # Stop simulation if running
        global sim
        if sim and sim.running:
            sim.stop()
            time.sleep(2)

        # Copy to remote and clean local files
        if args.pathremote != "":
            sys_call("cp {}/cg_analysis.log {}/ || true".format(args.path, args.pathremote))
            sys_call("cp {}/cg_analysis.out {}/ || true".format(args.path, args.pathremote))
            if args.backend == MDEngine.DDCMD:
                sys_call("cp {}/ddcmd.out {}/ || true".format(args.path, args.pathremote))
                #sys_call("cp {}/data {}/ || true".format(args.path, args.pathremote))
                sys_call("cp {}/ddcMD.header {}/ || true".format(args.path, args.pathremote))
            elif args.backend == MDEngine.GROMACS or args.backend == MDEngine.GROMACS_PARTS:
                sys_call("cp {}/md.log {}/ || true".format(args.path, args.pathremote))
                sys_call("cp {}/mdrun.log {}/ || true".format(args.path, args.pathremote))

            # @TODO ideally here we want to copy over the latest GROMACS files so we can restart from current state
            # currently we don't do this so we can lose 10-20 min, issue here is if a big job stops all jobs will sync and 
            # potentially crash luster and could lead to bad copies
            #if args.backend == MDEngine.GROMACS_PARTS:
            # alt copyRemoteRestart that copies in progress latest .xtc    

            # @WARNING for all testing remember to comment this out                       
            # @TODO remove this don't cleanup locally - for now add this back in Helgi @TODO add back in in production
            # sys_call("rm -r "+args.path)

if __name__ == '__main__':
    main()
