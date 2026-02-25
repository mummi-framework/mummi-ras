"""
Python module for coordinating the analysis and running of MD simulations.

Specifically, this module handles the argument parsing for the analysis as well
as handles the IO and maintains the datastructures accumulated during the analysis run.

It also handles the file management, moving data from local node SSD and lustre filesystems.
"""


import argparse
import mummi_ras

class Analyzer:
    """
    The analyzer class coordinates the analysis and online feedback of a running MD simulaiton.
    """

    def __init__(self):
        """
        Constructor for analysis run.

        Reads arguments, set IO, sets options datastructure, sets up output buffer.
        """
        parser = self.setup_argparser()
        args = parser.parse_args()

        LOGGER.info("Running analysis, options are:")
        LOGGER.info(args)
        LOGGER.info("path        = {}".format(args.path))
        LOGGER.info("remote - gpfs - path = {}".format(args.pathremote))

        # Initilize iointerface 
        self.iointerface = mummi_ras.get_io(args.fstype)
        if args.nofeedback:
            self.fbinterface = None
        else:
            self.fbinterface = mummi_ras.get_io(args.fbio)
        LOGGER.info('Using IO_Interface [%s] for data output', self.iointerface)
        LOGGER.info('Using IO_Interface [%s] for feedback output', self.fbinterface)

        # Read in all arguments
        self.simname = args.simname
        self.path = args.path
        self.test = args.test
        self.noninteractive = args.noninteractive
        self.fcounter = read_arg(args.fcount)  # Warning this has to be correct first frame
        self.pathremote = args.pathremote
        self.fstype = args.fstype
        self.fbio = args.fbio
        self.nt = args.nprocs

        self.currentCheckpointIndex = 0

        # @TODO this should be in resorce yaml files - and read in here.
        # @WARNING this needs to match what ddcMD is doing - use 100000 for
        self.counterIncr = 25000
        # @WARNING this needs to match what ddcMD is doing - e.g. !!!
        # don't !!! remove restart files
        self.counterSaveIncr = 1000
        #self.counterSaveIncr = 100000
        # @WARNING this needs to match what ddcMD is doing - set in
        # object.data - checkpointrate
        self.counterCheckpointIncr = 500000
        # Stop Simulations and analysis after this max frame is reached
        # for C3 we expect to use 250000000 or 5us - about 5 days of running
        self.counterSimStop = args.maxsimtime
        # Min frame needs to be reached before stopping simulations and analysis with the flag
        self.counterSimStopMinLimit = args.minsimtime
        # Sets how many frames (snapshots, analysis and feedback) are processed together
        # This affects the GPFS load default is 5 for large runs something larger like 20 is suggested
        self.frameProcessBatchSize = args.frameProcessBatchSize

        self.dynamicsEngine = args.dynamicsEngine

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

        parser.add_argument("--minsimtime", type=int, default=100000000,
                            help="Don't stop simulation and analysis before this time step is reached")

        parser.add_argument("--frameProcessBatchSize", type=int, default=5,
                            help="Show many frames are processed together")

        # Simulation specific arguments
        parser.add_argument("--siminputs", type=str, required=True,
                            help="Path to collection of simulation inputs.")

        parser.add_argument("--simbin", type=str, required=True,
                            help="Path to simulation binary to use for simulation. Empty string is don't run simulation")

        parser.add_argument("--nofeedback", action="store_true",
                            help="Turn off the feedback calculation.")

        parser.add_argument("--useGromacs", action="store_true",
                            help="Turn off online ddcMD and use only gromacs files")

        parser.add_argument("--dynamicsEngine", type=str, required=True,
                            help="Specify which MD engine to use for this run")

        return parser

    def processOutput(options, cOutputList, iointerface):
        """Parse the cOutoutList dictionary and send feedback messages and write to analysis files."""

        # Sort output list by frame number - so they are processed in order
        cOutputList.sort(key=lambda x: x["frameNum"])

        # loop over frames 
        # this loops only processes data for the feedback mechanism. Rest of cOutput gets written later
        for cOutput in cOutputList:

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

            #if depth == 'n/a':
            #    depth = np.nan
            LOGGER.info(f'Add feedback for frameNum {cOutput["frameNum"]} - state {cOutput["states"]}')
            # we are removing this! as no atomistic in PC. Normally put in ROOT/ml/cg
            #options.cgs.add_snapshot(options.simname, cOutput["frameNum"], tilt, rot, depth)
            if options.fb != None:
                options.fb.aggregate(options.simname, cOutput["frameNum"], cOutput["prot_positions"], cOutput["states"])



        LOGGER.info("out of the loop")
        ## this is outside the "for cOutput in cOutputList" loop
        #if len(cOutputList[0]["states"]) > 0:
        # In C4 Only do any feedback if one protein in the patch - not for lipids only and not for dimers
        LOGGER.info('report now')
        if options.fb != None:
            options.fb.report()
        LOGGER.info('reported')

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

        return

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

        return


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

    return 


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

    return
