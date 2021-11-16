# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
###############################################################################
# @todo add Pilot2-splash-app disclaimer
# @authors  Harsh Bhatia  <bhatia4@llnl.gov>
#           Helgi I. Ingolfsson <ingolfsson1@llnl.gov>
###############################################################################

"""Module to read output of CG simulations."""

import numpy as np
#import datetime
import logging
import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
import mummi_ras.transformations.createsim.lipiddef as lipiddef

# Set mda defaults
#mda.core.flags['use_periodic_selections'] = True
#mda.core.flags['use_KDTree_routines'] = False

# This flag has to be false as lipids are -pbc mol
# (e.g. all lipid atoms are together untill COM moves across the PBC)
#mda.core.flags['use_pbc'] = False
#np.set_printoptions(precision=3)

LOGGER = logging.getLogger(__name__)


class CGParser:

    def __init__(self, fnametraj, fnametopo, empty=False):
        self.fnametopo = fnametopo
        self.fnametraj = fnametraj
        self.empty = empty

        if not empty:
            self.setup(fnametraj, fnametopo)

    def setup(self, fnametraj, fnametopo):
        LOGGER.info("Initializing CGParser")
        LOGGER.info("\tLoading: ({}) ({})".format(self.fnametraj, self.fnametopo))
        # ----------------------------------------------------------------------
        # Load files
        self.syst = mda.Universe(fnametopo, fnametraj)
        '''
        x = self.syst.atoms.select_atoms("resname CYF and name SC1")
        print x.atoms.positions
        '''

        LOGGER.info("\tNumber of frames in sim: {:d}".format(len(self.syst.trajectory)))
        LOGGER.info("\tNumber of atoms in sim: {:d}".format(len(self.syst.atoms)))
        LOGGER.info("\tSim dimensions: " + str(self.syst.dimensions[0:3]))
        #print "\tNumber of frames in sim: {:d}".format(len(self.syst.trajectory))
        #print "\tNumber of atoms in sim: {:d}".format(len(self.syst.atoms))
        #print "\tSim dimensions: " + str(self.syst.dimensions[0:3])

        # ----------------------------------------------------------------------
        # Get list of all supported lipid types
        lipidTypeNames = np.array([l[0] for l in lipiddef.lipidTypeList])

        # Get all resnames in this simulation
        self.resnames = np.unique(self.syst.atoms.resnames)

        # Make a lipid type dictionary
        LOGGER.warn(' Resnames ' +
            str(self.resnames[np.logical_not(np.in1d(self.resnames, lipidTypeNames))])
            + ' either not lipids or not supported!')

        self.lipidTypeDic = {}
        for i in np.where(np.in1d(lipidTypeNames, self.resnames))[0]:
            self.lipidTypeDic[lipidTypeNames[i]] = lipiddef.LipidType(lipiddef.lipidTypeList[i])

        # ----------------------------------------------------------------------
        # Define top/bottom leaflets

        # Get beads that define the bilayer leaflets
        # (not in bilayer middle, but closely packed, use linker + first tail bead of lipids that do not flip-flop)
        self.defHeadgroups = self.syst.atoms[[]]
        for i in list(self.lipidTypeDic.values()):
            self.defHeadgroups += i.getLeafletSelection(self.syst)

        # get leaflets
        # LeafletFinder .optimize give errors with very asymetrical bilayers - use ourown function
        #print(self.defHeadgroups)
        #rcutoff, n = optimize_cutoff(self.syst, self.defHeadgroups, dmin=15.0, dmax=20.0, step=0.5)
        #lfls = LeafletFinder(self.syst, self.defHeadgroups, cutoff=rcutoff, pbc=True)
        lfls = []
        rcutoff = 15
        while 1:
            lfls = LeafletFinder(self.syst, self.defHeadgroups, cutoff=rcutoff, pbc=True)
            if len(lfls.groups()) <= 2:
                break
            else:
                rcutoff += 0.5

        if len(lfls.groups()) != 2:
            error = "Simulation doesn't have x2 leaflets."
            LOGGER.error(error)
            raise ValueError(error)

        # check if they're even
        self.top_head, self.bot_head = lfls.groups()
        rt = float(len(self.top_head)) / float(len(self.bot_head))
        if rt > 1.5 or rt < 0.5:
            raise ValueError("Uneven leaflets.")

        LOGGER.warn('Found leaflets of size {} {}'.format(len(self.top_head), len(self.bot_head)))

        # Make top_head as the actual top/upper/outer leaflet
        if self.top_head.center_of_geometry()[2] < self.bot_head.center_of_geometry()[2]:
            self.top_head, self.bot_head = (self.bot_head, self.top_head)

        # Get beads that define the heads of all flip-flopping none leaflet defining lipids
        self.defFlipFlopHeadgroups = self.syst.atoms[[]]
        for i in list(self.lipidTypeDic.values()):
            self.defFlipFlopHeadgroups += i.getNoneLeafletSelection(self.syst)

        # Setup done now not empty
        self.empty = False

        # ----------------------------------------------------------------------
        # Done init!
        LOGGER.debug('__init__ done!')


    def update_frame(self, frameName):
        # self.syst.load_new(pdbName)  # This replaces this old cordinates with the new
        self.syst.load_new(frameName, format='DDCMD')  # This replaces this old cordinates with the new

    '''
    def get_layers(self, frame_idx):
        # this sets the syst frame to frame_idx for everything that references that syst
        ts = self.syst.trajectory[frame_idx]
        LOGGER.info("\tFrame: {:5d}, Time: {:8.3f}ps".format(ts.frame, self.syst.trajectory.time))
        return self.get_layers()
    '''

    def get_layers(self):
        # Get all lipids in top/bot leaflets (including flip-flop lipids - therefore has to be done for each frame)
        tp = self.top_head + self.defFlipFlopHeadgroups.select_atoms("around 12 global group topsel", topsel=self.top_head)
        bt = self.bot_head + self.defFlipFlopHeadgroups.select_atoms("around 12 global group botsel", botsel=self.bot_head)
        ce = (self.defHeadgroups + self.defFlipFlopHeadgroups) - bt - tp

        # @TODO How long does this take? Haven't ever seen it fail
        # This one should never fail if "normal" lipids don't flip-flop and the bilayers are defect free
        if len(tp.select_atoms("group bt", bt=bt)):
            LOGGER.warn("\tFrame %d: #%d common atoms between leaflets.\n" %
                (self.syst.trajectory.frame, len(tp.select_atoms("group bt", bt=bt))))

        #print 'returning from get_layers!'
        return tp, ce, bt

    '''
    def get_positions(self, frame_idx):
        tp, bt = self.get_layers(frame_idx)
        frame = []
        for layer, layer_name in zip((tp, bt), ('outer', 'inner')):
            for resname in np.unique(layer.resnames):
                atoms = layer.select_atoms('resname ' + resname)
                positions = atoms.positions[:,0:2]
                frame.append((positions, resname, layer_name))

        return (frame, frame_idx)
    '''
