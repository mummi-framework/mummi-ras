#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------

import os
import enum
import logging
import time
import logging
from multiprocessing import  Pool
import numpy as np
import MDAnalysis as mda
from typing import List

from mummi_core.utils.naming import MuMMI_NamingUtils as  Naming
from mummi_ras.datastructures.cgSnapshot import CGSnapshot
from mummi_ras.interfaces import get_io
from mummi_ras.online.cg.get_tiltrot_z_state_multiple import get_protein_info, getKRASstates
from mummi_core.utils.utilities import sig_ign_and_rename_proc

# dynim 0.1 in Spack does not have dynim.SamplerBinned so we import local dynim
# try:
#     import dynim
# except ImportError:
from ... import dynim

LOGGER = logging.getLogger(__name__)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
class UCGSelectorType(enum.Enum):
    Unknown = 0
    Worker = 1
    Manager = 2

# ------------------------------------------------------------------------------
class UCGSelector(object):

    # name of the selector (will have simname and hostname)
    namefmt = '{}'  # removed hostname and cg-selector

    # --------------------------------------------------------------------------
    def __init__(self, type: UCGSelectorType, simname: str, workspace: str, config: dict):

        assert isinstance(workspace, str) and isinstance(config, dict)
        assert isinstance(type, UCGSelectorType)
        assert type in [UCGSelectorType.Manager, UCGSelectorType.Worker]

        self.type = type
        self.simname = simname
        self.name = UCGSelector.namefmt.format(simname)
        self.workspace = workspace
        self.io = get_io('simple')

        if type == UCGSelectorType.Worker:
            self.last_frame = -1

        else:
            min_cands_b4_sel = int(config['min_cands_b4_sel'])
            novelty_factor = float(config['novelty_factor'])

            assert min_cands_b4_sel >= 0
            assert 0. <= novelty_factor <= 1.
            assert simname == 'manager'

            self.sampler = dynim.SamplerBinned(self.name, workspace,
                                               novelty_factor=novelty_factor,
                                               binning_type='3d-spherical_and_cartesian',
                                               min_cands_b4_sel=min_cands_b4_sel)
            self.recent_selections = []

            self.proc_path = os.path.join(self.workspace, "processed_frames")
            if not os.path.isdir(self.proc_path):
                os.mkdir(self.proc_path)

    # --------------------------------------------------------------------------
    def __str__(self):
        if self.type == UCGSelectorType.Manager:
            return f'<{self.name}: #candidates = {self.sampler.num_candidates()}>'
        else:
            return f'<{self.name}: last_frame = {self.last_frame}>'

    def __repr__(self):
        return self.__str__()

    # --------------------------------------------------------------------------
    def num_candidates(self):
        return self.sampler.num_candidates()

    def num_selected(self):
        return self.sampler.num_selected()

    def update(self):
        self.sampler.update()

    # --------------------------------------------------------------------------
    def add_snapshot(self, structure_file: str, topology_file: str, cleanup: bool = False):
        """
        Add a valid frame (snapshot) to the selection process.

        :param structure_file: Path to the structure file (.gro)
        :type structure_file: str
        :param topology_file: Path to the corresponding topology file (.tpr)
        :type topology_file: str
        :param cleanup: If True, we delete the topology file after computing the tilt, rotation and depth.
        :type cleanup: bool
        """

        simname = os.path.splitext(structure_file)[0]
        frame_id = int(self._remove_but_digits(simname.split('_')[-1]))
        LOGGER.debug(f"fnamestrct = {structure_file}, fnametopo = {topology_file}, simname = {simname}")

        syst = mda.Universe(topology_file, structure_file)

        # Obtain trd info from backmapped structures 
        kras_indices = get_protein_info(syst, 'resname ACE1 and name BB')
        # Here for UCG lipids_exist MUST be True
        states_list = getKRASstates(syst, kras_indices, lipids_exist = False)

        tilt = np.deg2rad(states_list[0][0][2])
        rot = np.deg2rad(states_list[0][0][3])
        depth = states_list[0][0][4]

        self._add_snapshot(simname, frame_id, tilt, rot, depth)
        # We delete the topology files
        if cleanup:
            try:
                os.remove(topology_file)
            except OSError as e:
                LOGGER.error(f"Could delete file {e.filename} => {e.strerror}")
                pass

    # --------------------------------------------------------------------------
    def _add_snapshot(self, simname: str, frame_id: int, tilt: float, rot: float, depth: float):
        """
        XXX

        :param simname: XXX
        :type simname: str
        :param frame_id: XXX
        :type frame_id: str
        :param tilt: XXX
        :type tilt: float
        :param rot: XXX
        :type rot: float
        :param depth: XXX
        :type depth: float
        """
        LOGGER.debug(f'adding [sim {simname}, frame({frame_id})]: t={tilt} r={rot} d={depth}')

        cg = CGSnapshot(simname, frame_id, tilt, rot, depth)
        if cg.data is None:
            LOGGER.warning(f'skipping [sim {simname}, frame({frame_id})]: t={tilt} r={rot} d={depth}')
            return

        # for a worker, write this file
        if self.type == UCGSelectorType.Worker:
            assert simname == self.simname

            fname = Naming.ucgframe(self.simname, frame_id)
            self.io.save_npz(self.workspace, fname, cg, CGSnapshot.save_npz)
            self.last_frame = max(self.last_frame, frame_id)

        # for a manager, add to datastructure
        else:
            hdpoints = [dynim.HDPoint(Naming.ucgframe(cg.simname, cg.frame_id), cg.data)]
            self.sampler.add_candidates(hdpoints)

    # --------------------------------------------------------------------------
    def add_snapshots(self, snapshots):
        n = len(snapshots)
        if n == 0 or self.type == UCGSelectorType.Worker:
            return

        LOGGER.info(f'adding {n} snapshots')
        snapshots = list(filter(lambda x: CGSnapshot.is_valid(x[2],x[3],x[4]),
                                snapshots))
        if len(snapshots) < n:
            LOGGER.warning(f'filtered {n-len(snapshots)} invalid snapshots')

        snapshots = [CGSnapshot(d[0], d[1], d[2], d[3], d[4]) for d in snapshots]

        # for a manager, add to datastructure
        hdpoints = [dynim.HDPoint(Naming.ucgframe(cg.simname, cg.frame_id), cg.data) for cg in snapshots]
        self.sampler.add_candidates(hdpoints)

    # --------------------------------------------------------------------------
    @staticmethod
    def _load_cg_snapshot(filename):
        try:
            _cg = CGSnapshot.load_npz(filename)
            _id = Naming.ucgframe(_cg.simname, _cg.frame_id)
            assert os.path.basename(filename) == _id+'.npz', f'Invalid id (= {_id}) read for key {_cg.simname}'
            return dynim.HDPoint(_id, _cg.data)

        except Exception:
            return None

    def update_manager(self, max_keys = 300000):

        # TODO: also apply the same fix
        #keypattern = Naming.cgframe('*', '*')
        keypattern = '*_0*_f0*'

        keys = np.array(self.io.list_keys(self.workspace, keypattern))
        keys = keys[:max_keys]
        nkeys = len(keys)

        LOGGER.info(f'Updating manager. found {nkeys} keys in ({self.workspace})')
        if nkeys == 0:
            return

        # ----------------------------------------------------------------------
        ts = time.strftime("%Y%m%d-%H%M%S")
        outpath = os.path.join(self.proc_path, f"processed.{ts}")
        os.makedirs(outpath, exist_ok=True)

        load_pool = Pool(processes=10, initializer=sig_ign_and_rename_proc, initargs=('pool_load_cgSelector',))
        hdpoints = []
        chunk_size = 10000

        try:
            for chunk_start in range(0, nkeys, chunk_size):
                keys_chunk = keys[chunk_start:chunk_start + chunk_size]
                LOGGER.debug(f'processing {len(keys_chunk)} frames starting with {chunk_start} (total={nkeys})')

                # parallel: load frames
                fnames = [os.path.join(self.workspace, k) for k in keys_chunk]
                frames = load_pool.map(UCGSelector._load_cg_snapshot, fnames, chunksize=1000)
                frames = np.array(frames)

                # figure out which ones were successful
                success_idxs = np.where(frames != None)[0]

                # remove the ones that failed to load
                if len(success_idxs) < len(keys_chunk):
                    LOGGER.warning(f'found {len(keys_chunk) - len(success_idxs)} frames = None')
                    frames = frames[success_idxs]
                    keys_chunk = keys_chunk[success_idxs]

                hdpoints.append(frames)

                LOGGER.debug(f'moving {len(keys_chunk)} keys from ({self.workspace}) to ({outpath})')
                for key in keys_chunk:
                    try:
                        os.rename(os.path.join(self.workspace, key), os.path.join(outpath, key))
                    except Exception as _expt:
                        LOGGER.error(f'Failed to move ({key}): {_expt}')

        finally:
            load_pool.close()
            load_pool.join()

        # now add these to the sampler
        hdpoints = np.concatenate(hdpoints)
        self.sampler.add_candidates(hdpoints) #, do_filter=False)
        # self.sampler.update()
        LOGGER.profile(f'Updated manager with {len(hdpoints)} candidate hdpoints')
        self.checkpoint()

    # --------------------------------------------------------------------------
    def select(self, k: int, do_confirm: bool = True) -> List[dynim.HDPoint]:
        """
        Select k points from the pool of available structures.

        :param k: Number of structure to select
        :type k: int
        :param do_confirm: Confirm the selection (see Dynim documentation)
        :type do_confirm: bool

        :return: The selected high-dimensional points
        :rtype: List[dynim.HDPoint]
        """
        assert isinstance(k, int)
        assert k >= 0

        # Sep 05, 2021
        # HB moved the sampler update from add_candidates() to select()
        self.sampler.update()

        LOGGER.info(f'Selecting {k} UCG frames from ({self.name})')

        self.recent_selections = self.sampler.select(k, do_confirm)
        selections = self.recent_selections

        # no need to maintain this if i am confirming
        if do_confirm:
            self.recent_selections = []

        LOGGER.debug(f'Selected {len(selections)} UCG frames from ({self.name})')
        self.checkpoint()
        return selections

    def _remove_but_digits(self, s: str) -> str:
        """ Remove everything but digits from a string"""
        return ''.join([x for x in s if x.isdigit()])

    # --------------------------------------------------------------------------
    def test(self):
        return self.sampler.test()

    def checkpoint(self):
        self.sampler.checkpoint()

    def restore(self):
        return self.sampler.restore()

# ------------------------------------------------------------------------------
