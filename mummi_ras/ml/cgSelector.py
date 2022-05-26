# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
import os
import enum
import logging
import time
import logging
from multiprocessing import  Pool
import numpy as np

from mummi_core.interfaces.utils import move_to, tar_and_remove
from mummi_core.utils.naming import MuMMI_NamingUtils as  Naming
from mummi_ras.datastructures.cgSnapshot import CGSnapshot
from mummi_ras.interfaces import get_io
from mummi_core.utils.utilities import sig_ign_and_rename_proc

if __package__ or "." in __name__:
    from . import dynim
else:
    import dynim


LOGGER = logging.getLogger(__name__)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
class CGSelectorType(enum.Enum):
    Unknown = 0
    Worker = 1
    Manager = 2


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
class CGSelector(object):

    # name of the selector (will have simname and hostname)
    namefmt = '{}'  # removed hostname and cg-selector

    # --------------------------------------------------------------------------
    def __init__(self, type: CGSelectorType, simname: str, workspace: str, config: dict):

        assert isinstance(workspace, str) and isinstance(config, dict)
        assert isinstance(type, CGSelectorType)
        assert type in [CGSelectorType.Manager, CGSelectorType.Worker]

        self.type = type
        self.simname = simname
        self.name = CGSelector.namefmt.format(simname)
        self.workspace = workspace
        self.io = get_io('simple')

        if type == CGSelectorType.Worker:
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
        if self.type == CGSelectorType.Manager:
            return f'<{self.name}: #candidates = {self.sampler.num_candidates()}>'
        else:
            return f'<{self.name}: last_frame = {self.last_frame}>'

    def __repr__(self):
        return self.__str__()

    # --------------------------------------------------------------------------
    def num_candidates(self):
        return self.sampler.num_candidates()

    def update(self):
        self.sampler.update()

    # --------------------------------------------------------------------------
    def add_snapshot(self, simname, frame_id, tilt, rot, depth):

        LOGGER.debug(f'adding [sim {simname}, frame({frame_id})]: t={tilt} r={rot} d={depth}')

        cg = CGSnapshot(simname, frame_id, tilt, rot, depth)
        if cg.data is None:
            LOGGER.warning(f'skipping [sim {simname}, frame({frame_id})]: t={tilt} r={rot} d={depth}')
            return

        # for a worker, write this file
        if self.type == CGSelectorType.Worker:
            assert simname == self.simname

            fname = Naming.cgframe(self.simname, frame_id)
            self.io.save_npz(self.workspace, fname, cg, CGSnapshot.save_npz)
            self.last_frame = max(self.last_frame, frame_id)

        # for a manager, add to datastructure
        else:
            hdpoints = [dynim.HDPoint(Naming.cgframe(cg.simname, cg.frame_id), cg.data)]
            self.sampler.add_candidates(hdpoints)

    # --------------------------------------------------------------------------
    def add_snapshots(self, snapshots):
        assert 0
        n = len(snapshots)
        if n == 0:
            return

        LOGGER.info(f'adding {n} snapshots')
        snapshots = list(filter(lambda x: CGSnapshot.is_valid(x[2],x[3],x[4]),
                                snapshots))
        if len(snapshots) < n:
            LOGGER.warning(f'filtered {n-len(snapshots)} invalid snapshots')

        snapshots = [CGSnapshot(d[0], d[1], d[2], d[3], d[4]) for d in snapshots]

        # for a worker, write this file
        if self.type == CGSelectorType.Worker:
            assert False

        # for a manager, add to datastructure
        else:
            hdpoints = [dynim.HDPoint(Naming.cgframe(cg.simname, cg.frame_id), cg.data) for cg in snapshots]
            self.sampler.add_candidates(hdpoints)

    # --------------------------------------------------------------------------
    @staticmethod
    def _load_cg_snapshot(filename):
        try:
            _cg = CGSnapshot.load_npz(filename)
            _id = Naming.cgframe(_cg.simname, _cg.frame_id)
            assert os.path.basename(filename) == _id+'.npz', f'Invalid id (= {_id}) read for key {key}'
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
        #write_pool = Pool(processes=10, initializer=sig_ign_and_rename_proc, initargs=('pool_write_cgSelector',))

        hdpoints = []
        chunk_size = 10000

        try:
            for chunk_start in range(0, nkeys, chunk_size):
                keys_chunk = keys[chunk_start:chunk_start + chunk_size]
                LOGGER.debug(f'processing {len(keys_chunk)} frames starting with {chunk_start} (total={nkeys})')

                # parallel: load frames
                fnames = [os.path.join(self.workspace, k) for k in keys_chunk]
                frames = load_pool.map(CGSelector._load_cg_snapshot, fnames, chunksize=1000)
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

                #write_pool.apply_async(tar_and_remove, (self.proc_path, proc_name, self.workspace,
                #                                        keys_chunk, frames,
                #                                        CGSnapshot.save_npz))

        finally:
            load_pool.close()
            #write_pool.close()

            load_pool.join()
            #write_pool.join()

        # now add these to the sampler
        hdpoints = np.concatenate(hdpoints)
        self.sampler.add_candidates(hdpoints) #, do_filter=False)
        # self.sampler.update()
        LOGGER.profile(f'Updated manager with {len(hdpoints)} candidate hdpoints')
        self.checkpoint()

    # --------------------------------------------------------------------------
    def select(self, k, do_confirm=True):

        assert isinstance(k, int)
        assert k >= 0

        # Sep 05, 2021
        # HB moved the sampler update from add_candidates() to select()
        self.sampler.update()

        LOGGER.info(f'Selecting {k} CG frames from ({self.name})')

        self.recent_selections = self.sampler.select(k, do_confirm)
        selections = self.recent_selections

        # no need to maintain this if i am confirming
        if do_confirm:
            self.recent_selections = []

        LOGGER.debug(f'Selected {len(selections)} CG frames from ({self.name})')
        self.checkpoint()
        return selections

    # --------------------------------------------------------------------------
    def test(self):
        return self.sampler.test()

    def checkpoint(self):
        self.sampler.checkpoint()

    def restore(self):
        return self.sampler.restore()

    #def restore_lspace(self, filepath):
    #    self.sampler.hdspace.restore(filepath)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
