#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2024, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
#
# @authors  Harsh Bhatia  <bhatia4@llnl.gov>
#           Helgi I. Ingolfsson <ingolfsson1@llnl.gov>
#           Jeremy O. B. Tempkin <tempkin1@llnl.gov>
#           Loic Pottier <pottier1@llnl.gov>
# -----------------------------------------------------------------------------

import numpy as np
import os
import logging

import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder

import mummi_ras.transformations.createsim.lipiddef as lipiddef
import mummi_ras.parsers.cg_aligner as cga 

LOGGER = logging.getLogger(__name__)


class CGParser:
    """Class that manages interface to an underlying MDAnalysis Universe with the purpose to read output of CG simulaitons. 
    
    This class implements io, setup / teardown, alignment and returns MDA selections required to compute 
    observables for other modules in MuMMI. 
    """

    def __init__(self, fnametraj:str , fnametopo: str, align: bool = True, empty: bool = False, gmx_trajfile: str = None):
        """Constructor for a CGParser object.
        """

        self.fnametopo = fnametopo
        self.fnametraj = fnametraj
        self.empty = empty
        self.align = align
        self.syst = None
        self.gmx_trajfile = gmx_trajfile
        # Record the last time the XTC (gmx_trajfile) was modified
        self.fnametraj_ts = -1

        if self.align:
            self.cga = cga.CGFrameAligner("membrane", True)

        if not os.path.isfile(self.fnametraj):
            LOGGER.info(f"Frame file {self.fnametraj} does not exist yet.")
            self.empty = True

        if not self.empty:
            self.setup(fnametraj, fnametopo)

        return

    def isSetup(self):
        return self.empty == False

    def setup(self, fnametraj, fnametopo):
        """Main setup routine for initializing CGParser components.
        """

        LOGGER.info(f"Initializing CGParser (MDAnalysis version {mda.__version__})")
        LOGGER.info("\tLoading: ({}) ({})".format(self.fnametraj, self.fnametopo))
        # ----------------------------------------------------------------------
        # Load files
        self.syst = mda.Universe(fnametopo, fnametraj, refresh_offsets=True)

        LOGGER.info("\tNumber of frames in sim: {:d}".format(len(self.syst.trajectory)))
        LOGGER.info("\tNumber of atoms in sim: {:d}".format(len(self.syst.atoms)))
        LOGGER.info("\tSim dimensions: " + str(self.syst.dimensions[0:3]))

        # Get all resnames in this simulation
        self.resnames = np.unique(self.syst.atoms.resnames)

        # setup lipids
        self._setup_lipids()

        # aligner setup
        self._setup_cg_aligner()

        # Setup done now not empty
        self.empty = False

        LOGGER.debug('__init__ done!')

        return

    def _setup_lipids(self):
        """Setup routine for lipid definitions and leaflets"""

        # ----------------------------------------------------------------------
        # Get list of all supported lipid types
        lipidTypeNames = np.array([l[0] for l in lipiddef.lipidTypeList])


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

        return

    def remove_duplicate_frame(self, fname: str, keep_original: bool = True):
        """
        Remove duplicate frames in a MDA universe. It does it by rewriting a file
        if keep_original is True, the original XTC is kept as fname.orig
        """
        if not self.syst:
            return

        temp = f"new.{fname}"
        uniq_idx = {}
        for x in self.syst.trajectory:
            if x.time not in uniq_idx:
                uniq_idx[x.time] = x.frame

        if len(uniq_idx) == self.syst.trajectory.n_frames:
            LOGGER.info(f"No duplicate frame found in {fname} ({self.syst.trajectory.n_frames})")
        else:
            LOGGER.warning(f"{self.syst.trajectory.n_frames-len(uniq_idx)} duplicate frames found in {fname}")

        slice_idx = np.array([x for x in uniq_idx.values()], dtype=np.int64)

        t = self.syst.atoms
        t.write(temp, frames = self.syst.trajectory[slice_idx])

        if keep_original:
            os.rename(fname, f"{fname}.orig")
            # LOGGER.info(f"listdir after rename of {fname}: {os.listdir()}")
        os.rename(temp, fname)
        # LOGGER.info(f"listdir after rename of {temp} into {fname}: {os.listdir()}")

    def _setup_cg_aligner(self):
        """Setup the aligner"""

        # setup the aligner.
        self.cga.set_mdu(self.syst)

        # set the protein configuration in aligner
        self.anchor_resids = cga.alignment_anchors['RAS-BRAF']

        self.cga.set_anchors(self.anchor_resids["origin"],
                             self.anchor_resids["xaxis"],
                             self.anchor_resids["xyaxis"])

        self.cga.init()

        return

    def get_layers(self):
        """Return the MDA selections correspinding to the top and bottom leaflets."""
        # Get all lipids in top/bot leaflets (including flip-flop lipids - therefore has to be done for each frame)
        tp = self.top_head + self.defFlipFlopHeadgroups.select_atoms("around 12 global group topsel", topsel=self.top_head)
        bt = self.bot_head + self.defFlipFlopHeadgroups.select_atoms("around 12 global group botsel", botsel=self.bot_head)
        ce = (self.defHeadgroups + self.defFlipFlopHeadgroups) - bt - tp

        # @TODO How long does this take? Haven't ever seen it fail
        # This one should never fail if "normal" lipids don't flip-flop and the bilayers are defect free
        if len(tp.select_atoms("group bt", bt=bt)):
            LOGGER.warning("\tFrame %d: #%d common atoms between leaflets.\n" %
                (self.syst.trajectory.frame, len(tp.select_atoms("group bt", bt=bt))))

        #print 'returning from get_layers!'
        return tp, ce, bt

    def get_selection(self, selection):
        """Return the MDAnalysis selection for the passed string.
        """
        return self.syst.select_atoms(selection)

    def align_frame(self):
        """Align the current frame.
        """

        self.cga.align()

        return

