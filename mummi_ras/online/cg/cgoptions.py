#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------
"""Python module for setting up options for CGanalysis"""

import yaml
import os
from mummi_core.utils import Naming
from mummi_ras.feedback.manager.feedback_manager_cg2ml import FeedbackManager_cg2ml

from mummi_core.workflow.feedback_manager import FeedbackManagerType
from mummi_core.utils.utilities import read_arg, sys_call
from mummi_ras.ml.selectors import CGSelector, CGSelectorType
from mummi_ras.online.cg.cgselections import CgSelections
from mummi_ras.online import MDEngine

# set up logging
from logging import getLogger
LOGGER = getLogger(__name__)

# Class for neede CgAnalyser options and helper functions
# @WARNING this one should be light as it's passed to the worker threads
class CgAnalyzerOptions:
    """A class for carrying around cganalysis parameters as well as some utility functions."""

    def __init__(self, args, iointerface, fbinterface):

        self.simname = args.simname
        self.path = args.path
        self.test = args.test
        self.noninteractive = args.noninteractive
        self.fcounter = read_arg(args.fcount)  # Warning this has to be correct first frame
        self.pathremote = args.pathremote
        self.fstype = args.fstype
        self.fbio = args.fbio
        self.nt = args.nprocs # Number of cores for analysis
        self.nsims_cores = args.simcores # Number of cores for simulation
        # hardcoded lookup for yaml at canonical file location
        self.parameter_file = os.path.join(Naming.dir_res("analysis"), "cg", "cg_analysis_parameters.yaml")

        LOGGER.info(f"Loading parameters from {self.parameter_file}")

        self.iointerface = iointerface
        self.fbinterface = fbinterface

        self.currentCheckpointIndex = 0

        # The backend used (GROMACS, DDCMD or NOSIM)
        self.MDEngine = args.backend

        if self.MDEngine == MDEngine.DDCMD:
            ### SPECIFIC TO DDCMD - Should be moved into specialized subclass
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
            ### END to SPECIFIC TO DDCMD
        else:
            # Specific to GROMACS
            self.xtcTrajFile = os.path.join(self.path, "traj_comp.xtc")
            self.counterSaveIncr = 1000
            self.counterIncr = 25000
            self.bundlingsize = 10                          # How many frames are in each .xtc bundel
            self.counterCheckpointIncr = 25000

        # Stop Simulations and analysis after this max frame is reached
        # for C3 we expect to use 250000000 or 5us - about 5 days of running
        self.counterSimStop = args.maxsimtime
        # Min frame needs to be reached before stopping simulations and analysis with the flag
        self.counterSimStopMinLimit = args.minsimtime
        # Sets how many frames (snapshots, analysis and feedback) are processed together
        # This affects the GPFS load default is 5 for large runs something larger like 20 is suggested
        self.frameProcessBatchSize = args.frameProcessBatchSize

        # If using GROMACS or testing make this False if ddcMD real analysis
        # use True
        if self.test:
            # Be very verbose about testing
            LOGGER.info("WARNING! THIS IS A TEST RUN ONLY")
            self.MDEngine = MDEngine.NOSIM
        

        if self.MDEngine == MDEngine.DDCMD:
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
        if self.fbinterface == None:
            self.fb = None
        else:
            self.fb = FeedbackManager_cg2ml(FeedbackManagerType.Worker, self.simname,
                                           cg_dir, cg_dir, cg_dir,
                                           self.iointerface, self.fbinterface, nbins = 100)

        self.cgs = CGSelector(CGSelectorType.Worker, self.simname,
                              os.path.join(Naming.dir_root('ml'), 'cg'), {})
        # --------------------------------------------------------------------------

        # read interface file 
        with open(self.parameter_file, "r") as fh:
            self.params = yaml.load(fh, Loader=yaml.Loader)

        # set up interfaces  
        self._build_interface_selections()

        return

    def get_interface_selections(self):
        """Return a list of tuples contaning the (selection_i, selection_j) 
        pairs for each interface we are tracking."""
        
        return self.interface_atom_ids 

    def _build_interface_selections(self):
        """Set up selections from interface specifications.""" 

        interface_sel = []
        for key in self.params['interfaces']:
            interface = self.params['interfaces'][key] 
            atom_id_1 = interface['domain_1']['atom_ids'] 
            atom_id_2 = interface['domain_2']['atom_ids'] 
            interface_sel.append((atom_id_1, atom_id_2))

        self.interface_atom_ids = interface_sel
        return 

    def set_cgSelections(self, cg):
        """Set cgSelections"""

        self.cgSel = CgSelections(cg)

        return

    def get_frame_folder_name(self, fcounter):
        """Return the folder name on disk for the given frame."""
        return 'snapshot.{:012d}'.format(fcounter)

    def get_sim_file_name(self, fcounter, extension):
        """Return the simulation name corresponding to the frame."""

        if self.MDEngine == MDEngine.DDCMD:
            # fname = 'snapshot.{:08d}/atoms#000000'.format(fcounter)
            fname = "{}/subset#000000" \
                .format(self.get_frame_folder_name(fcounter))
            if len(extension) > 0:
                fname = fname + '.{}'.format(extension)
        else:
            # GROMACS case
            fname = self.get_sim_analysis_name(fcounter, extension)
        fname = os.path.join(self.path, fname)
        return fname

    def get_sim_analysis_name(self, fcounter: int, extension = "xtc"):
        """Return the name of the simulaiton we are analyzing."""
        if self.MDEngine == MDEngine.GROMACS:
            return f"traj_comp.xtc"
        elif self.MDEngine == MDEngine.GROMACS_PARTS:
            intCount = int(fcounter * self.counterIncr)
            return f"traj_comp.step{intCount:012d}.{extension}"