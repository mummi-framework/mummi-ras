"""
Feedback manager for handling converting CG simulation data for ML.
"""


import numpy as np
import os
import collections

from mummi_core.workflow.feedback_manager import FeedbackManager, FeedbackManagerType

from logging import getLogger
LOGGER = getLogger(__name__)


class FeedbackManager_cg2ml(FeedbackManager):  

    def __init__(self, type, name, workspace, outpath, fbpath, iointerface, fbinterface, nbins = 100):
        """
        Constructor for a cg to ml manager.
        """
        assert type in [FeedbackManagerType.Manager, FeedbackManagerType.Worker]
        super().__init__(type, name)

        # defines header path
        self.workspace = workspace

        # outpath for the feedback (where final files are written in npz)  
        self.iointerface = iointerface
        self.outpath = outpath

        # fbpath where communication between worker and manager will happen
        self.fbinterface = fbinterface
        self.fbpath = fbpath

        # give me a full path as the fb path
        # and i will figure out if i need to trim down to the dbr namespace
        assert len(os.path.dirname(self.fbpath)) > 0
        assert self.fbpath == self.outpath

        self.fname_chkpt = f'{self.name}.chk'
        self.fname_chkpt_data = f'{self.name}.chk.npz'

        self.fname_master = os.path.join(self.fbpath, 'xyz_master_history.tar')

        self.structure_buffer = collections.deque()

        # set to zero at construction
        self.nstructures = 0

        # Number of bins (i.e, folders created in feedback-cg2ml/)
        # Each bins/folder will get all the structure which id is: structure_ID % nbins
        self.nbins = nbins

        # ----------------------------------------------------------------------
        LOGGER.info(f'Initialized {super().__str__()}')

        return

    # report the analysis
    def report(self):

        LOGGER.info(f'Writing out new {len(self.structure_buffer)} new structures to disk')

        for p in self.structure_buffer:
            _simname = p[0]
            _fid = p[1]
            _pos = p[2]
            _state = p[3]

            #import ipdb; ipdb.set_trace()
            # we spread feedback into nbins folder to speedup filesystem access
            sim_id = int(_simname.split('_')[2]) % self.nbins
            sim_id_str = f"{sim_id:03}"
            outpath = os.path.join(self.outpath, sim_id_str)

            self.iointerface.save_npz(outpath,
                                      f"{_simname}_f{_fid:012}.npz",
                                      {"positions" : _pos,
                                       "state" : _state},
                                      interfaces=['simple'])

            # increment record of structures
            self.nstructures += 1

        # we wrote all structures from th ebuffer to disk, now clean it out  
        self.structure_buffer.clear()
        return

    # report the analysis
    def report(self):
        _simlist = []
        _data = []
        simnames = set()
        LOGGER.debug(f'Writing out new {len(self.structure_buffer)} new structures to disk')
        for p in self.structure_buffer:
            _simname = p[0]
            simnames.add(_simname)
            _fid = p[1]
            _pos = p[2]
            _state = p[3]
            _simlist.append(f"{_simname}_f{_fid:012}.npz")
            _data.append({"positions" : _pos, "state" : _state})
            # increment record of structures
            self.nstructures += 1

        if len(simnames) != 1:
            LOGGER.error(f"This function cannot handle more than one structure at the time: {simnames}")
            raise RuntimeError(f"This function cannot handle more than one structure at the time: {simnames}")

        # we spread feedback into nbins folder to speedup filesystem access
        sim_id = next(iter(simnames)) # Hack to retrieve an element from a set without removing it
        sim_id = int(sim_id.split('_')[2]) % self.nbins
        sim_id_str = f"{sim_id:03}"
        outpath = os.path.join(self.outpath, sim_id_str)

        self.iointerface.save_npz(outpath, _simlist, _data, interfaces=['simple'])
        LOGGER.info(f"Added {len(_simlist)} frames (from {', '.join([x for x in simnames])}) to feedback ({self.outpath})")
        # we wrote all structures from the ebuffer to disk, now clean it out  
        self.structure_buffer.clear()
        return

    def aggregate(self, sn, fnum, pos, st):
        """
        Append structure to feedback write buffer.
        """

        assert isinstance(pos, np.ndarray)

        LOGGER.debug(f'Aggregating a new position: {self}.')

        self.structure_buffer.append([sn, fnum, pos, st])

        LOGGER.debug(f'Aggregated data: {self}')

        return

    def checkpoint(self):
        return

    def load(self):
        return

    def restore(self):
        return

    def test(self):
        return
