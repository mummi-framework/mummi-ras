###############################################################################
# @todo add Pilot2-splash-app disclaimer
# @authors  Harsh Bhatia  <bhatia4@llnl.gov>
#           Helgi I. Ingolfsson <ingolfsson1@llnl.gov>
#           Tim Carpenter <carpenter36@llnl.gov>
#           Christopher A. Neale <cneale@lanl.gov>
#           Francesco Di Natale <dinatale3@llnl.gov>
###############################################################################

"""Functions to analyse CG particle simulations."""

import argparse
from mummi_core.utils.logger import init_logger
import multiprocessing as mp
import numpy as np
import os
import time
import traceback

import mummi_core
from mummi_core.utils import Naming
from mummi_core.utils.utilities import read_arg, sys_call
import mummi_ras
from mummi_ras.online import ddcMD
# these require mdanalysis
import mummi_ras.online.cg.cnlibs as cnlibs
import mummi_ras.online.cg.fast_get_rdf
import mummi_ras.online.cg.get_tiltrot_z_state_multiple
from mummi_ras.datastructures.rdf import RDF
from mummi_ras.parsers.cg import CGParser

from mummi_ras.ml import CGSelector, CGSelectorType
from mummi_ras.feedback.feedback_manager_macro import MacroFeedbackManager
from mummi_core.workflow.feedback_manager import FeedbackManagerType

sim = None # global var for ddcMD sim

from logging import getLogger
LOGGER = getLogger(__name__)

# ----------------------------------------------------------------------------
'''
## TODO: this is a temporary function that Harsh added on Sep 17, 2020
## in lieu of the code that actually computes the rot-tilt-depth
## we will just read the existing data and randomly sample from the set!
def read_trajectories():

    file = '/p/gpfs1/helgi/campaign1x_data/cesar_analysis_init2/RAS-RAF-TILT-ROT.xvg'
    print ('reading ({})'.format(file))

    data = np.loadtxt(file)
    print ('read {} data'.format(data.shape))
    return data
# ----------------------------------------------------------------------------
'''

"""
# ----------------------------------------------------------------------------
    This module analysis cg data (ddcMD or GROMACS) and copies back to file
    system

    For testing use e.g.
    # @TODO add a working test in a shared test dir

    workon sierrapysplash/
    source ~/pilot2-splash-app/setup/setup.env.sh
    export MUMMI_ROOT=/p/gpfs1/helgi/tests-splashApp/cganalysis-01/small-01
    cd /p/gpfs1/helgi/tests-splashApp/cganalysis-01/
    cp -r small-clean-remote/ small-05-remote/
    cp -r small-clean/ small-05/
    cd small-05/
    cganalysis --fstype simple --nscgpatches ns:cgpatches \
    --nscgrdfs ns:cg-rdfs --nspecies 14 --ngrid 5 \
    --nscg2bkilled ns:cg2bkilled \
    --path /p/gpfs1/helgi/tests-splashApp/cganalysis-01/small-05  \
    --pathremote /p/gpfs1/helgi/tests-splashApp/cganalysis-01/small-05-remote \
    -l /p/gpfs1/helgi/tests-splashApp/cganalysis-01/small-05 \
    --simname small-01 -c -d 1 --fcount 100000 --nprocs 3 \
    >> analysis-test7.log 2>&1 &

# ----------------------------------------------------------------------------
"""


# This class as predefined CG MDAnalysis selections use later
class CgSelections(object):
    """PLACEHOLDER"""

    def __init__(self, cg):

        # @TODO this should probably also be in a yaml file?
        # Constants for analysis
        # cutoff for protein-protein contacts (in Angstroms)
        self.contactcutoff = 6.0
        # cutoff for protein-lipid contacts (in Angstroms)
        self.contactcutoffLong = 12.0
        self.numResPerProtein = 184
        self.numResPerProteinGTPMg = 186
        # self.numAtomPerProtein = 418
        # self.numAtomPerProteinGTPMg = 428
        self.numResGdomain = 165
        self.proteinIdentResidueString = "THR"
        self.numProteinIdentResiduePerProtein = 26

        # Can fail on first call
        try:
            self.pdbRef = \
                mummi_ras.online.cg.get_tiltrot_z_state_multiple \
                .kras_ref_universe
            # This is a MDAnalysis universe with a reference KRAS only loaded
            # x1 time
        except ImportError as e:
            print(e)
            print(
                "> ERROR: Splash installation not found. Please handle "
                "pdbRef explicitly!")

        # Define the number of proteins in the system, has been taken out of the getKRASstates function so only has to be done once
        self.kras_indices = mummi_ras.online.cg.get_tiltrot_z_state_multiple.get_protein_info(cg.syst,'resname ACE1 and name BB')

        # find the first protein residue and number of proteins
        self.firstProteinStartResidue, self.numProteins = \
            cnlibs.findFirstProteinResidue(
                cg.syst, self.proteinIdentResidueString,
                self.numProteinIdentResiduePerProtein)

        # groupContactProt, groupPlContactLipid = \
        #   cnlibs.getSel_findProtLipidContacts(
        #       cg.syst, numProteins, firstProteinStartResidue,
        #       numResPerProteinGTPMg)
        # groupContactProt, groupContactProtBB, groupPlContactLipid = \
        #   cnlibs.getSel_findProtLipidContacts(
        #        cg.syst, numProteins, firstProteinStartResidue,
        #        numResPerProteinGTPMg)
        self.groupContactProt, self.groupContactG, self.groupContactHVR, \
            self.groupContactProtBB, self.groupPlContactLipid = \
            cnlibs.getSel_findProtLipidContacts(
                cg.syst, self.numProteins, self.firstProteinStartResidue,
                self.numResPerProteinGTPMg, self.numResPerProtein,
                self.numResGdomain)
        # could have getSel_findProtLipidContacts() return
        # groupContactAllProt but Chris' other branch is ahead of develop so
        # will leave this as a separate routine for now
        # groupContactAllProt = cnlibs.getSel_allProteinFromAllRas(
        #    pdbPatch, numProteins, firstProteinStartResidue,
        #    numResPerProteinGTPMg)

        self.groupLipidHeadSelection, self.groupLipidTailASelection, \
            self.groupLipidTailBSelection = \
            cnlibs.getSel_findLipidLeaflet(cg.syst)

        self.groupGdomainBB = cnlibs.getSel_findGdomainBB(
            cg.syst, self.numProteins, self.firstProteinStartResidue,
            self.numResPerProteinGTPMg, self.numResGdomain)
        self.groupRefGdomainBB = cnlibs.getSel_findGdomainBB(
            self.pdbRef, 1, 1, self.numResPerProteinGTPMg, self.numResGdomain)

        # options are 'serial' and 'OpenMP'
        self.distanceBackendType = 'Serial'


# Class for neede CgAnalyser options and helper functions
# @WARNING this one should be light as it's passed to the worker threads
class CgAnalyzerOptions:
    """PLACEHOLDER"""

    def __init__(self, args, iointerface, fbinterface):

        self.simname = args.simname
        self.path = args.path
        self.test = args.test
        self.noninteractive = args.noninteractive
        self.fcounter = read_arg(args.fcount)  # Warning this has to be correct first frame
        self.pathremote = args.pathremote
        self.fstype = args.fstype
        self.fbio = args.fbio
        self.nt = args.nprocs

        self.iointerface = iointerface
        self.fbinterface = fbinterface

        self.currentCheckpointIndex = 0

        # @TODO this should be in resorce yaml files - and read in here.
        # @WARNING this needs to match what ddcMD is doing - use 100000 for
        # teting
        self.counterIncr = 25000
        # @WARNING this needs to match what ddcMD is doing - e.g. !!!
        # don't !!! remove restart files
        self.counterSaveIncr = 100000
        # @WARNING this needs to match what ddcMD is doing - set in
        # object.data - checkpointrate
        self.counterCheckpointIncr = 500000
        # Stop Simulations and analysis after this max frame is reached
        # for C3 we expect to use 250000000 or 5us - about 5 days of running
        self.counterSimStop = args.maxsimtime
        # Sets how many frames (snapshots, analysis and feedback) are processed together
        # This affects the GPFS load default is 5 for large runs something larger like 20 is suggested
        self.frameProcessBatchSize = args.frameProcessBatchSize

        # If using GROMACS or testing make this False if ddcMD real analysis
        # use True
        if self.test:
            # Be very verbose about testing
            LOGGER.info("WARNING! THIS IS A TEST RUN ONLY")
            self.ddcMD = False
        else:
            self.ddcMD = True

        if self.ddcMD:
            self.fileformat = ""
        else:
            # "pdb" is default can also do "gro"
            self.fileformat = "gro"

        # Set .tpr file - use for MDAnalysis first read
        self.fnameTopo = "topol.tpr"
        self.fnameTopo = os.path.join(self.path, self.fnameTopo)

        # CG selections
        self.cgSel = None

        # --------------------------------------------------------------------------
        # temporary code for feedback and cg selection

        cg_dir = Naming.dir_root('feedback-cg')

        self.fb = MacroFeedbackManager(FeedbackManagerType.Worker, self.simname,
                                       cg_dir, cg_dir, cg_dir,
                                       self.iointerface, self.fbinterface, None)

        self.cgs = CGSelector(CGSelectorType.Worker, self.simname,
                              os.path.join(Naming.dir_root('ml'), 'cg'), {})

        # --------------------------------------------------------------------------

    def set_cgSelections(self, cg):
        """PLACEHOLDER"""

        self.cgSel = CgSelections(cg)

    def get_frame_folder_name(self, fcounter):
        """PLACEHOLDER"""

        return 'snapshot.{:012d}'.format(fcounter)

    def get_sim_file_name(self, fcounter, extension):
        """PLACEHOLDER"""

        if self.ddcMD:
            # fname = 'snapshot.{:08d}/atoms#000000'.format(fcounter)
            fname = "{}/subset#000000" \
                .format(self.get_frame_folder_name(fcounter))
            if len(extension) > 0:
                fname = fname + '.{}'.format(extension)
        else:
            fname = self.get_sim_analysis_name(fcounter, extension)
        fname = os.path.join(self.path, fname)
        return fname

    def get_sim_analysis_name(self, fcounter, extension):
        """PLACEHOLDER"""

        return "{}_{}_f{}.{}" \
            .format(self.simname, 'cg000', fcounter, extension)


def processOutput(options, cOutputList, iointerface):
    """PLACEHOLDER"""

    # Sort output list by frame number - so they are processed in order
    cOutputList.sort(key=lambda x: x["frameNum"])

    '''  @TODO need to add back feedback update
    for cOutput in cOutputList:
        # Update current RDF
        rdfs = np.zeros(RDF.SHAPE, dtype=RDF.DTYPE)
        cState = cOutput["states"][0][0][-1]
        # Note need to skip first line of rdf (these are the bins)
        if cState == 1:
            rdfs[0] = cOutput["KRAS_lipid_rdfs"][1:, :]
        elif cState == 2:
            rdfs[1] = cOutput["KRAS_lipid_rdfs"][1:, :]
        else:
            error = \
                "Analysis output gave state {} only sates 1 and 2 supported" \
                .format(cState)

            LOGGER.error(error)
            raise ValueError(error)
        #fb.accumulate_frame(rdfs)
    '''

    for cOutput in cOutputList:
        # temp code for cg selection
        ##Time, Tilting, Rotation, CRD distance - skip first line as it's a comment
        #random_trd = options.trd_data[np.random.randint(1, options.trd_data.shape[0])]
        #tilt = np.deg2rad(random_trd[1])
        #rot = np.deg2rad(random_trd[2])
        #depth = random_trd[3]
        #LOGGER.info('Add fake feedback for frameNum {} - tilt {} / rot {} / depth {} from {} '.format(cOutput["frameNum"], tilt, rot, depth, random_trd))
        #options.cgs.add_snapshot(options.simname, cOutput["frameNum"], tilt, rot, depth)

        # Only do any feedback if some protein in patch (i.e. lipid only patches are currently not used for feedback)
        # if len(cOutput["states"]) > 0:
        # In C4 Only do any feedback if one protein in the patch - not for lipids only and not for dimers
        if len(cOutput["states"]) == 1:

            # Add possible snapshot to CG to AA consideration - only do this if this patch is x1 RAS-RAF - nothing else
            LOGGER.debug('should we add this to ml? {} {} {}'.format(len(cOutput["states"]) == 1, cOutput["states"][0][0][0] == "RAS-RAF", cOutput["frameNum"] % options.counterSaveIncr == 0))
            LOGGER.debug('should we add this to ml? num proteins {}, type first protein {}, frame num {}'.format(len(cOutput["states"]), cOutput["states"][0][0][0], cOutput["frameNum"]))
            #if (len(cOutput["states"]) == 1) and (cOutput["states"][0][0][0] == "RAS-RAF") and ((cOutput["frameNum"] % options.counterSaveIncr) == 0):
            # For C4 CG to AA selection is done for all monomers (of all x4 protein types)
            if (len(cOutput["states"]) == 1) and ((cOutput["frameNum"] % options.counterSaveIncr) == 0):
                state = cOutput["states"][0][0][5]
                tilt = np.deg2rad(cOutput["states"][0][0][2])
                rot = np.deg2rad(cOutput["states"][0][0][3])
                # @WARNING only do if not RAS-RAF  - check
                depth = cOutput["states"][0][0][4]
                if depth == 'n/a':
                    depth = np.nan
                LOGGER.info('Add feedback for frameNum {} - protein {} / state {} / tilt {} / rot {} / depth {}'.format(cOutput["frameNum"], cOutput["states"][0][0][0], state, tilt, rot, depth))
                options.cgs.add_snapshot(options.simname, cOutput["frameNum"], tilt, rot, depth)

            ''' Refactor - use code form Tomas
            HardCoded = False
            if HardCoded:
                ## Map state designations to indices: rasmap,rafmap
                rasmap = { "a":0 , "b":1 , "c":2 , "za":3 , "ma":4 , "mb":5 }
                rafmap = { "za":6 , "ma":7 , "mb":8 }
            else:
                ## Compute ra{s,f}map from STATE_NAMES, to avoid hard-coding.
                rasmap = {}
                rafmap = {}
                idx = 0
                for x in STATE_NAMES:
                    if x.find("RAS") >= 0:
                        rasmap[x.replace("RAS","")] = idx
                    else:
                        rafmap[x.replace("RAF","")] = idx
                    idx = idx + 1
            ras_rdf = cOutput["RAS-ONLY_lipid_rdfs"][1:]  # Remove first linee ...
            idx = rasmap[ cOutput["states"][0][0][5] ]
            options.rdfs_total[idx] = options.rdfs_total[idx] + ras_rdf
            options.rdfs_counts[idx] += 1
            if cOutput["states"][0][0][0] == "RAS-RAF":
                raf_rdf = cOutput["RAS-RAF_lipid_rdfs"][1:]  # Remove first line ...
                idx = rafmap[ cOutput["states"][0][0][5] ]
                options.rdfs_total[6] = options.rdfs_total[6] + raf_rdf
                options.rdfs_counts[6] += 1
            '''
            ''' old code
            rdfs = np.zeros(RDF.SHAPE, dtype=RDF.DTYPE)
            ## RAS_RAF rdfs
            ##@TODO: need to fix back to RAS-RAF for actual Campaign 4
            if cOutput["states"][0][0][0] == "RAS4A-RAF":
                ras_rdf = cOutput["RAS-ONLY_lipid_rdfs"][1:]  # Remove first line that is the x axes numbers
                if cOutput["states"][0][0][5] == 'za':
                    rdfs[3] = ras_rdf
                elif cOutput["states"][0][0][5] == 'ma':
                    rdfs[4] = ras_rdf
                elif cOutput["states"][0][0][5] == 'mb':
                    rdfs[5] = ras_rdf
                rasraf_rdf = cOutput["RAS-RAF_lipid_rdfs"][1:]  # Remove first line that is the x axes numbers
                if cOutput["states"][0][0][5] == 'za':
                    rdfs[6] = rasraf_rdf
                elif cOutput["states"][0][0][5] == 'ma':
                    rdfs[7] = rasraf_rdf
                elif cOutput["states"][0][0][5] == 'mb':
                    rdfs[8] = rasraf_rdf
            ## RAS_Only rdfs
            else:
                ras_rdf = cOutput["RAS-ONLY_lipid_rdfs"][1:]  # Remove first line that is the x axes numbers
                if cOutput["states"][0][0][5] == 'a':
                    rdfs[0] = ras_rdf
                elif cOutput["states"][0][0][5] == 'b':
                    rdfs[1] = ras_rdf
                elif cOutput["states"][0][0][5] == 'c':
                    rdfs[2] = ras_rdf
            '''
            rdfs = np.zeros(RDF.SHAPE, dtype=RDF.DTYPE)
            statemapRAS = {"RAS-ONLY"  : { "a" : 0 , "b" : 1 , "c" : 2  },
                           "RAS-RAF"   : { "za": 3 , "ma": 4 , "mb": 5  },
                           "RAS4A-ONLY": { "a" : 9 , "b" : 10, "c" : 11 },
                           "RAS4A-RAF" : { "za": 12, "ma": 13, "mb": 14 }}
            statemapRAF = {"RAS-RAF"   : { "za": 6 , "ma": 7 , "mb": 8  },
                           "RAS4A-RAF" : { "za": 15, "ma": 16, "mb": 17 }}
            ras_rdf = cOutput["RAS-ONLY_lipid_rdfs"][1:]  # Remove first line that is the x axes numbers
            idx = statemapRAS[ cOutput["states"][0][0][0] ][ cOutput["states"][0][0][5] ]
            rdfs[idx] = ras_rdf
            if "RAF" in cOutput["states"][0][0][0]: # RAS-RAF for 4A or 4B
                rasraf_rdf = cOutput["RAS-RAF_lipid_rdfs"][1:]  # Remove first line that is the x axes numbers
                idx = statemapRAF[ cOutput["states"][0][0][0] ][ cOutput["states"][0][0][5] ]
                rdfs[idx] = rasraf_rdf
            options.fb.aggregate(rdfs=rdfs)

    LOGGER.info("out of the loop")
    ## this is outside the "for cOutput in cOutputList" loop
    #if len(cOutputList[0]["states"]) > 0:
    # In C4 Only do any feedback if one protein in the patch - not for lipids only and not for dimers
    if len(cOutputList[0]["states"]) == 1:
        LOGGER.info('report now')
        options.fb.report()
        options.fb.reset()
        LOGGER.info('reported and reset')

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
    iointerface.save_snapshots(saveList, options.path, cSavePath)

    # Copy over all logs - @WARNING these can be ok big and very vastfull to to this here - can remove if needed.
    if options.pathremote != "":
        copyRemoteLogs(options)

    # @Warning snapshots are always cleared - restart only if copied above
    # (using pathremote)
    cleanupLocal(options, listAll)


#def run_main_analysis(options, cg, fb, cgs, iointerface, lipid_types, simulation):
def run_main_analysis(options, cg, iointerface, lipid_types, simulation):

    cfcounter = options.fcounter
    # only used for GROMACS as frames only read/loaded in first round
    reading = True
    cFrameList = []
    cOutputList = []    # List of frames that have alreaddy been analyzed
    cProcList = np.empty((options.nt), dtype=object)
    outputQueue = mp.Queue()

    # keep looking for new frames
    # check that the simulation is stil running
    restarts_left = 0
    if not options.noninteractive:
        restarts_left = simulation.restart_limit

    while True:
        if not options.noninteractive:
            # Check if the simulation is successful.
            if simulation.successful:
                # We should stop the simulation
                simulation.stop()
            elif not simulation.running:
                # Check if the simulation is still running.
                if restarts_left > 0:
                    simulation.setup_restart()
                    restarts_left -= 1
                    simulation.run()

        # -----------------------------------------------------
        # Read in all new files (currently readdy in file system)
        if options.ddcMD:
            while True:
                # Get ddcMD frame
                fnameTraj = \
                    options.get_sim_file_name(cfcounter, options.fileformat)

                LOGGER.info("Checking for frame {}".format(fnameTraj))
                if os.path.isfile(fnameTraj):  # file not found try later
                    cFrameList.append(cfcounter)
                    cfcounter = cfcounter + options.counterIncr
                else:
                    break
        elif cFrameList == [] and reading:
            # Assumign this is a GROMACS trajectory so initiliza as all frames
            cFrameList = list(range(0, len(cg.syst.trajectory)))
            reading = False
            LOGGER.info(
                "Getting all frame numbers in trajectory from 0 to {}".format(len(cg.syst.trajectory)))

        # -----------------------------------------------------
        # Process a single frame
        if len(cFrameList) > 0:
            if options.nt == 1:
                cFrameNum = cFrameList.pop(0)
                LOGGER.info('Process frame {}'.format(cFrameNum))
                cOut = processFrame(cFrameNum, cg, options, lipid_types)
                outputQueue.put(cOut)

                time.sleep(0.5)
            elif options.nt > 1:
                # Find empty or finished precessor
                for i in range(0, len(cProcList)):
                    if cProcList[i] is None or not cProcList[i].is_alive():
                        cFrameNum = cFrameList.pop(0)
                        cProcList[i] = mp.Process(target=threadProcessFrame, args=(cFrameNum, cg, options, outputQueue, lipid_types))
                        LOGGER.info('Process frame {} in subprocess num / name {} / {} '.format(cFrameNum, i, cProcList[i].name))
                        cProcList[i].start()
                        # break for loop; only start one frame per outer loop
                        break

        # -----------------------------------------------------
        # Get output from alreaddy processed frames - as soon as it's ready
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

        # -----------------------------------------------------
        # Wait before next run
        time.sleep(1) # @TODO Addjsut depending on need - @NOTE this is a bussy loop

        # -----------------------------------------------------
        # Stop criteria is max number of frames
        # If reached, stop sim, process last frames and exit
        if not options.noninteractive and options.fcounter >= options.counterSimStop:

            LOGGER.info(
                "Sim {} is past {} steps, will now be stopped.".format(options.simname, options.counterSimStop))

            # Stop DDCMD run
            simulation.stop()

            # Copy everyhing to remote
            if len(cOutputList) > 0:
                processOutput(options, cOutputList, iointerface)

            # Mark simulations as converged - will be read by workflow - and
            # copy other logs to remote
            success_flag, _ = Naming.status_flags('cg')
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
        # if cfcounter > 12 and len(cFrameList) == 0 \
        #    and outputQueue.qsize() <= batchSize:
        #    break

        # If list to process is empty and not in interactive mode - exit while loop and wait for last snapshots.
        if len(cFrameList) <= 0 and options.noninteractive:
            break

    # -------------------------------------------
    # Exited main while - now finalize wait for sub-process to finish and
    # process last data
    LOGGER.info("CG analysis is stopping.")

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


def copyRemoteLogs(options):
    """ PLACEHOLDER """

    # Copy over all log files - warning these can be big so don't do all the
    # time. @TDOO we want to do this in a smarter way !!!
    try:
        LOGGER.info("Copy over log files")
        copyDirList = ""
        copyDirList = copyDirList + " " + os.path.join(options.path, "data")
        copyDirList = copyDirList + " " + os.path.join(options.path, "ddcMD.header")
        copyDirList = copyDirList + " " + os.path.join(options.path, "hpm.data")
        copyDirList = copyDirList + " " + os.path.join(options.path, "ddcmd.out")
        copyDirList = copyDirList + " " + os.path.join(options.path, "cg_analysis.out")
        # @TODO don't copy the python log if it gets big
        copyDirList = copyDirList + " " + os.path.join(options.path, "*.log")
        sys_call('cp -r '+copyDirList+' '+options.pathremote+'/ || true')
    except Exception as e:
        LOGGER.error("cganalysis - copy to remote error {}".format(e))


def cleanupLocal(options, frameList):
    """Cleanup local files"""

    # @WARNING currently this only works for ddcMD file structure
    if not options.ddcMD:
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


def copyRemoteRestart(options, frameList):
    """Copy over ddcMD restart files + any analysis (in folders)"""

    # @WARNING currently this only works for ddcMD file structure
    if not options.ddcMD:
        LOGGER.warning("Copy to remote only works for ddcMD files")
        return

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
        sys_call('rsync -va --no-g '+copyDirList+' '+options.pathremote+'/ --exclude subset#000000*')
        # skip setup files and their .pdb's

        # After copy remove local data files - so not to fillup local dir
        #sys_call("rm -r "+copyDirList)

        # After copy if needed create a ddcMD restart link in remote dir
        # (@NODE the link ddcMD creates points to dir that are not copied yet)
        if (newCheckpointIndex > options.currentCheckpointIndex):
            sys_call('ln -sf ./{}/restart {}/restart'.format(options.get_frame_folder_name(newCheckpointIndex), options.pathremote))
            options.currentCheckpointIndex = newCheckpointIndex

    except Exception as e:
        LOGGER.error("cganalysis - copy to remote error %s", e.message)


def processFrame(cFrameNum, cg, options, lipid_types):
    """Run all analysis returns a dictionary of all results."""

    lastTime = time.time()
    return_dic = {"frameNum": cFrameNum}

    cFrameName = options.get_sim_file_name(cFrameNum, options.fileformat)
    ddcMDFrameName = cFrameName

    # if ddcMD: @TODO add back in to support GROMACS format - XXX needed
    # We need to convert ddcMD frame to pdb and use the pdb name
    cFrameName = cFrameName+'.pdb'
    # @TODO fix - bad hack to reduce load on nfs
    #sys_call('martiniobj2pdb -i '+ddcMDFrameName+' -o '+cFrameName+' -t /usr/gapps/kras/templates/ddcMD/proItpList')
    #sys_call('python /p/gpfs1/splash/programs/martiniobj2pdb -i '+ddcMDFrameName+' -o '+cFrameName+' -t '+options.path+'/proItpList')

    if cg.empty:
        # This will not happen if initlized first
        # cg not initilzed yet, then this is the first time called and we
        # have to create
        cg.setup(cFrameName, options.fnameTopo)
    else:
        LOGGER.info("Load frame len {} framre {}".format(len(cg.syst.trajectory), cFrameName))
        if len(cg.syst.trajectory) > cFrameNum:
            # Trajectory alreaddy loaded, set MDAnalysis to correct frame
            cg.syst.trajectory[cFrameNum]
        else:
            # Normal for a ddcMD run load each rame into CGTrajt
            # cg parser shold alreaddy be initilized - we need to load in new
            # frame and update leaflets
            cg.update_frame(ddcMDFrameName)   # (cFrameName)

    bbox = [0, cg.syst.dimensions[0], 0, cg.syst.dimensions[1]]
    tp, ce, bt = cg.get_layers()
    LOGGER.info(
        "Frame: {0:5d}, Time: {1:8.3f} ps, Bbox = {2:s}"
        .format(cFrameNum, cg.syst.trajectory.time, str(bbox)))
    LOGGER.info("Lengths of layers: {} {} {}".format(len(tp), len(ce), len(bt)))
    # NOTE: if more information is needed for anlaysis,
    # we can add functionalities to fetch it from CGParser

    # @TODO Removed for now - due to instlal issues with patch - do we want this down the line
    '''
    #####################################################
    # A1:    patch creation
    # NOTES: simtime is going to be zero from individual pdb files!
    #        this will create only a single patch
    pc = PatchCreator(options.ngrid, options.nspecies)
    patches = pc.create_cg(
        cFrameName, options.fnameTopo, options.simname, float(cFrameNum),
        cFrameNum, bbox, (tp, ce, bt), lipid_types)
    LOGGER.debug(
        "Created %s patch of version %s", len(patches), cgPatch.version)
    assert(len(patches) == 1)
    p = patches[0]
    return_dic["patch"] = p
    '''
    return_dic["patch"] = np.zeros(1) # All return value shoudl be numpy arrays


    # ####################################################
    # A2:    KRAS states
    # For each KRAS return (KRAS number, t_1 angle, t_2 angle, z dist from
    # bilayer, full state[1-4], main state[1-2])
    # output file is saved as numpy matix (number of KRAS x 6)
    statesList = \
        mummi_ras.online.cg.get_tiltrot_z_state_multiple.getKRASstates(cg.syst, options.cgSel.kras_indices)
    # Update state number to state name - maybe shoudl be part of the files ???
    #print("HII-before " + str(statesList))

    ## Refactor - use code form Tomas
    statemap = {"RAS-ONLY" : { 1:"a"  , 2:"b"  , 3:"c"  },
                "RAS-RAF"  : { 1:"ma" , 2:"mb" , 3:"za" },
                "RAS4A-ONLY": {1: "a", 2: "b", 3: "c"},
                "RAS4A-RAF": {1: "ma", 2: "mb", 3: "za"}}

    for pState in statesList:
        for innerState in pState:
            try:
                innerState[5] = statemap[innerState[0]][innerState[5]]
            except:
                LOGGER.error("-get, States, protein type {} state {} not supported".format(innerState[0], innerState[5]))
                raise Exception("Error in protein type definition {}, state {} not supported".format(innerState[0], innerState[5]))
    return_dic["states"] = statesList
    '''
    for pState in statesList:
        for innerSate in pState:
            if innerSate[0] == "RAS-ONLY":
                if innerSate[5] == 1:
                    innerSate[5] = "a"
                elif innerSate[5] == 2:
                    innerSate[5] = "b"
                elif innerSate[5] == 3:
                    innerSate[5] = "c"
                else:
                    LOGGER.error("-get, States, protein type {} state {} not supported".format(innerSate[0], innerSate[5]))
                    raise Exception("Error in protein type definition {}, state {} not supported".format(innerSate[0], innerSate[5]))
            elif innerSate[0] == "RAS-RAF":
                if innerSate[5] == 1:
                    innerSate[5] = "ma"
                elif innerSate[5] == 2:
                    innerSate[5] = "mb"
                elif innerSate[5] == 3:
                    innerSate[5] = "za"
                else:
                    LOGGER.error("-get, States, protein type {} state {} not supported".format(innerSate[0], innerSate[5]))
                    raise Exception("Error in protein type definition {}, state {} not supported".format(innerSate[0], innerSate[5]))
            else:
                LOGGER.error("-get, States, protein type {} not supported".format(innerSate[0]))
                raise Exception("Error in protein type definition, type {} not supported".format(innerSate[0]))
    return_dic["states"] = statesList
    '''
    #print("HII-after " + str(statesList))

    if len(statesList) > 0:
        LOGGER.info(
            "-get, protein states for shape {}, first is {}".format(statesList[0].shape, statesList[0][0]))
    else:
        # This happens if there are no proteins so a lipid only patch
        LOGGER.info("-get, we have a lipid only patch")

    # @TODO add RDF's back in after fixing MDanalysis issue with baseFastSingleFrame.py and fast_rdf2d.py
    # ####################################################
    # A3: Get RDFs like counts from RKAS (central only) to all lipids in outer leaflet
    if len(statesList) > 0:  # Skip if no proteins in patch
        cRDF = mummi_ras.online.cg.fast_get_rdf.RAS_lipid_RDF(cg.syst, bt)
        return_dic["RAS-ONLY_lipid_rdfs"] = cRDF
        LOGGER.debug("-get, RAS RDFs size {} count {}".format(cRDF.shape, np.sum(cRDF[1:])))

        # A3b: If central protein is RAS-RAF also get RDFs like counts from RAF-CRD (central only) to all lipids
        if return_dic["states"][0][0][0] == "RAS-RAF" or return_dic["states"][0][0][0] == "RAS4A-RAF":
            cRDF = mummi_ras.online.cg.fast_get_rdf.RAF_lipid_RDF(cg.syst, bt)
            return_dic["RAS-RAF_lipid_rdfs"] = cRDF
            LOGGER.debug("-get, RAF(CRD) RDFs size {} count {}".format(cRDF.shape, np.sum(cRDF[1:])))

    '''
    # #################################################### -->
    # Now 0.2 seconds with 5 ras (contacts based on all-to-all beads) --
    # vs 0.064 seconds with old routine that was less informative
    # A4:    Get protein-lipid contacts
    results = cnlibs.findProtLipidContacts(
        cg.syst, options.cgSel.groupContactProtBB,
        options.cgSel.groupPlContactLipid, options.cgSel.numProteins,
        options.cgSel.firstProteinStartResidue,
        options.cgSel.numResPerProteinGTPMg, options.cgSel.contactcutoffLong,
        options.cgSel.distanceBackendType)
    return_dic["plContacts"] = results

    # ####################################################
    # A5: Get G domain - HVR intramolecular contacts
    results = cnlibs.findGHVRContacts(
        cg.syst, options.cgSel.groupContactG, options.cgSel.groupContactHVR,
        options.cgSel.numProteins, options.cgSel.firstProteinStartResidue,
        options.cgSel.numResPerProteinGTPMg, options.cgSel.contactcutoff,
        options.cgSel.distanceBackendType)
    return_dic["ghContacts"] = results

    # ########################################### --> 0.027 seconds with 5 Ras
    # A6:    Get protein-protein intermolecular contacts
    results = cnlibs.findProtProtContacts(
        cg.syst, options.cgSel.groupContactProt, options.cgSel.numProteins,
        options.cgSel.firstProteinStartResidue,
        options.cgSel.numResPerProteinGTPMg, options.cgSel.contactcutoff,
        options.cgSel.distanceBackendType)
    return_dic["ppContacts"] = results

    # ######################################### --> 0.0015 seconds with 5 Ras
    # A7:    Get RMSD of G domain
    # find RMSD surrogate of G domain to reference structure: this is actually
    # a RMSD of the intermolecular pairwise distance matrix (much faster to
    # compute)
    # output file has one column per Ras, with the RMSD to the reference for
    # selected residues (units=A)
    results = cnlibs.findGdomainRMSDsurrogate(
        options.cgSel.groupGdomainBB, options.cgSel.groupRefGdomainBB,
        options.cgSel.numProteins, options.cgSel.distanceBackendType)
    return_dic["rmsdGdomain"] = results

    # ######################################### --> 0.0006 seconds  with 5 Ras
    # A8: Get Rgyr of G domain (residues 1-165)
    # comment out three lines below since we do not need this (the RMSD
    # surrogate is a better sanity check)
    #     results = cnlibs.findGdomainRgyr(groupGdomainBB, numProteins)
    #     return_dic["rgyrGdomain"] = results

    # ########################################## --> 0.0004 seconds with 5 Ras
    # A9:    Get center of mass of G domain X,Y,Z (res Thr1 to His165)
    results = cnlibs.findGdomainCOM(
        options.cgSel.groupGdomainBB, options.cgSel.numProteins)
    return_dic["comGdomain"] = results

    # ########################################### --> 0.046 seconds with 5 Ras
    # A10: Get what leaflet each lipid is in. This is not expected to be exact
    # frame-by-frame, but a smoothed version should show if there are
    # flip-flops
    results = cnlibs.findLipidLeaflet(
        options.cgSel.groupLipidHeadSelection,
        options.cgSel.groupLipidTailASelection,
        options.cgSel.groupLipidTailBSelection)
    return_dic["lipidLeafletDesignation"] = results
    '''

    # #################################################### --> 0.001 seconds
    # (XTC) or 0.93 seconds (PDB) with 5 Ras
    # A10: get compressed trajectory with coordinates of one atom per lipid
    # and BB atoms of protein
    # BB is protein and F1 is farnesyl tail bead for analysis (but also keep
    # F2, F3, F4 for visualization)
    # C1A bead is in DIPE, POPC, POPE, POPX; D1A bead is in PAP6, PAPC, PAPS;
    # T1A bead is in DPSM; R1 bead is in CHOL
    # MG+ bead is magnesium; resname GTP
    # @WARNING remove for now to save space in DBR - if we put this back in -
    # uncomment groupCompressedTrajSelection above
    #   groupCompressedTrajSelection = \
    #       cnlibs.getSel_writeCompressedSelection(cg.syst)
    #   return_dic["compressed.xtc"] = results

    # #################################################### --> 0.00005 seconds
    # A11: Get box dimensions (cg.syst.dimensions has the same values as
    # cg.syst.dimensions)
    results = cg.syst.dimensions[0:3]
    return_dic["box"] = results

    LOGGER.info(
        "Time for analysis - frameNum %s, time %s s",
        cFrameNum, time.time()-lastTime)

    # Example of return values - used for feadback
    # (cFrameNum, int(statesList[0][0][-1]), cRDF, p)
    # 0 return_dic["frameNum"]
    # 1 return_dic["states"][0][0][-1]
    # 2 return_dic["KRAS_lipid_rdfs"]
    # 3 return_dic["patch"]

    # Return all data as a dictonary
    return return_dic


def threadProcessFrame(cFrameNum, cg, options, outputQueue, lipid_types):
    try:
        outputQueue.put(processFrame(cFrameNum, cg, options, lipid_types))
        LOGGER.info(".. Frame done {}".format(str(cFrameNum)))
    except Exception as e:
        outputQueue.put({"frameNum": cFrameNum})
        LOGGER.error("cganalysis - frame {} - sup process error {}".format(str(cFrameNum), e))
    return


def setup_argparser():
    """Set up the program's argument parser."""
    parser = argparse.ArgumentParser(description='Analyze CG simulations individually.')

    parser.add_argument('--simname', required=True,
                        help='Name of the simulation')

    parser.add_argument("-f", "--fcount", type=int, default=1,
                       help="Frame number for inital frame (default: 1)")

    parser.add_argument("--noninteractive", action="store_true",
                       help="Set if running on a previusly run folder - not with live simulations")

    parser.add_argument("-p", "--path", required=True, type=str,
                       help="Path for input and output")

    parser.add_argument("-pr", "--pathremote", type=str, default="",
                       help="If set periotically copy files to remote WARNING this will also remove everything in path on exit")

    parser.add_argument('--fstype', type=str, default="mummi", choices=mummi_ras.get_interfaces(),
                        help='Chose how to write files/data (default: mummi)')

    parser.add_argument('--fbio', type=str, default="mummi", choices=mummi_ras.get_interfaces(),
                        help='Chose how to write files/data (default: mummi)')

    parser.add_argument('--nprocs', type=int, default=1,
                        help='Numer of processes to use for analysis (default: 1)')

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

    parser.add_argument("--frameProcessBatchSize", type=int, default=5,
                        help="Show many frames are processed together")

    # Simulation specific arguments
    parser.add_argument("--siminputs", type=str, required=True,
                        help="Path to collection of simulation inputs.")

    parser.add_argument("--simbin", type=str, required=True,
                        help="Path to simulation binary to use for simulation. Empty string is don't run simulation")

    return parser


def run_cganalysis(args):

    # Show options
    LOGGER.info("Running CG analysis, options are:")
    LOGGER.info(args)
    LOGGER.info("path        = {}".format(args.path))
    LOGGER.info("remote - gpfs - path = {}".format(args.pathremote))

    # Initilize ouput DBR/files etc
    iointerface = mummi_ras.get_io(args.fstype)
    fbinterface = mummi_ras.get_io(args.fbio)

    LOGGER.info('Using IO_Interface [%s] for data output', iointerface)
    LOGGER.info('Using IO_Interface [%s] for feedback output', fbinterface)

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

    if options.ddcMD:
        # Check init .gro file (only one per simulation but it has to be there)
        fname = "lipids-water-eq4.gro"
    else:
        # If not ddcMD we assume this is a GROMASC trajectory - .xtc file
        # provied as the sim name
        fname = options.simname
    fname = os.path.join(options.path, fname)
    if not os.path.isfile(options.fnameTopo):
        error = "Simulation init .gro file {} not found".format(fname)
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
    cg = CGParser(fname, options.fnameTopo, empty=False)
    # Lazy laod code not used - so we make this a lazy load (so each sub-p can
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
    options.set_cgSelections(cg)

    # Hardcoding lipid types here for now. We will likely want to parse this out of an input.
    lipid_types = {
        "nInner": ["POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"],
        "nOuter": ["POPX", "PAPC", "POPE", "DIPE", "DPSM", "CHOL"],
    }

    # Simulation setup
    # NOTE: We're initializing here, as it looks like the run_main_analysis
    # function has the main while True loop in it.
    # TODO: We need to refactor to pull out the while true loop to an actual
    # main method. It's bad practice to bury a while true loop.
    global sim
    if not options.noninteractive:
        success_flag, fail_flag = Naming.status_flags('cg')
        out_path = options.pathremote if options.pathremote else options.path
        _pass = iointerface.test_signal(out_path, success_flag)
        _fail = iointerface.test_signal(out_path, fail_flag)
        LOGGER.info("Testing flags: [Success='%s', Fail='%s']", success_flag, fail_flag)
        LOGGER.info("Checking for exit condition [%s] (pass=%s fail=%s)", out_path, _pass, _fail)
        if _pass or _fail:
            LOGGER.info("DETECTED EXIT CONDITION: Pass=%s Fail=%s", _pass, _fail)
            exit(0)
        sim = run_ddcmd_simulation(args)

    # Find frames and run analysis
    # @WARNING for ddcMD jobs this is an infinit loop and never returns.
    # @TODO add ddcMD process id as an argument and in run_main_analysis check if still running if not report, copy back and terminate.
    run_main_analysis(options, cg, iointerface, lipid_types, sim)


def run_ddcmd_simulation(args):
    sim = ddcMD(args.simname, args.siminputs, args.path, args.simbin, 2)
    # @TODO 2 is how often to try restart change to e.g. 5 in C3
    sim.initialize()
    sim.run()
    return sim

def main():
    # Set up the necessary base data structures to begin study set up.
    parser = setup_argparser()
    args = parser.parse_args()

    # Sets up MuMMI and the logger
    mummi_core.init()
    mummi_core.create_root()
    mummi_core.init_logger(argparser=args)

    try:
        run_cganalysis(args)
    except Exception as e:
        LOGGER.error(e)
        traceback.print_exc()
        raise e
    finally:
        # Stop ddcMD if running
        global sim
        if sim and sim.running:
            sim.stop()
            time.sleep(2)

        # Copy to remote and clean local files
        if args.pathremote != "":
            sys_call("cp {}/cg_analysis.log {}/ || true".format(args.path, args.pathremote))
            sys_call("cp {}/cg_analysis.out {}/ || true".format(args.path, args.pathremote))
            sys_call("cp {}/ddcmd.out {}/ || true".format(args.path, args.pathremote))
            #sys_call("cp {}/data {}/ || true".format(args.path, args.pathremote))
            sys_call("cp {}/ddcMD.header {}/ || true".format(args.path, args.pathremote))
            # @WARNING for all testing remember to comment this out
            # @TODO remove this don't cleanup locally - for now add this back in
            sys_call("rm -r "+args.path)


if __name__ == '__main__':
    main()
