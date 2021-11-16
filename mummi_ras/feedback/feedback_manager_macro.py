# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
import os
import numpy as np
from time import strftime, time
import time, datetime
import yaml
import multiprocessing
import random

from logging import getLogger
LOGGER = getLogger(__name__)

import mummi_core
from mummi_core.utils import Naming
from mummi_ras.datastructures.rdf import RDF
from mummi_core.workflow.feedback_manager import FeedbackManager, FeedbackManagerType


# ------------------------------------------------------------------------------
class MacroFeedbackManager(FeedbackManager):

    # --------------------------------------------------------------------------
    def __init__(self, type, name, workspace,
                 outpath, fbpath, iointerface, fbinterface, pselector):

        assert type in [FeedbackManagerType.Manager, FeedbackManagerType.Worker]
        super().__init__(type, name)

        if pselector is not None:
            pselector.test()
        self.pselector = pselector

        self.naggregates = 0            # number of aggregations
        self.nreports = 0               # number of reported aggregations
        self.naggs_at_last_report = 0   # number of aggregations in last report

        self.workspace = workspace

        # outpath for the feedback (where final files will be saved for macro)
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
        self.fname_chkpt_data = f'{self.name}.rdf.chk.npz'

        self.fname_master = os.path.join(self.fbpath, 'rdf_master_history.tar')

        self.rdfs = RDF()
        self.rdfs.set_zero()

        # ----------------------------------------------------------------------
        # 2020.10.22 add npz file from helgi
        # 2021.07.12 now loads all 18 x 8 x 599 (both RAS and RAS-RAF for both RAS 4B and 4A)
        if self.type == FeedbackManagerType.Manager:
            _fname = os.path.join(Naming.MUMMI_RESOURCES,
                                  'feedback_macro_initial_rdfs',
                                  'RDFs_counts_all_sims.npz')
            LOGGER.info(f'Initializing master RDF from ({_fname})')
            init_rdfs = np.load(_fname, allow_pickle=True)['rdfs_total']
            self.rdfs.set(init_rdfs)
            LOGGER.info(f'Initialized rdf: {self.rdfs.data.shape}')

        # ----------------------------------------------------------------------
        LOGGER.info(f'Initialized {super().__str__()}')

    def __str__(self):
        return f'({self.name}): ' \
               f'#reports = {self.nreports}; #aggregates = {self.naggregates}; ' \
               f'#aggregates_at_last_report = {self.naggs_at_last_report}; ' \
               f'sum_rdf = {self.rdfs.sum()}'

    def __repr__(self):
        return self.__str__()

    # --------------------------------------------------------------------------
    # load the data (do we need this as a separate function?)
    def load(self):

        # ----------------------------------------------------------------------
        if self.type == FeedbackManagerType.Manager:
            LOGGER.error('Not implemented!')
            return

        # ----------------------------------------------------------------------
        elif self.type == FeedbackManagerType.Worker:
            LOGGER.error('Not implemented!')
            return
        # ----------------------------------------------------------------------

    # --------------------------------------------------------------------------
    def reset(self):
        if self.type == FeedbackManagerType.Manager:
            LOGGER.warning('Should not be called for a Manager!')

        elif self.type == FeedbackManagerType.Worker:
            self.rdfs.set_zero()

    # --------------------------------------------------------------------------
    # aggregate the data
    def aggregate(self, rdfs=None, lock_patch_selector=None):

        # ----------------------------------------------------------------------
        if self.type == FeedbackManagerType.Manager:

            IO_READ = 'redis'  # read from this interface!

            servers_to_keys = \
                self.fbinterface.get_keys_for_macro_feedback(self.fbpath,
                                                             Naming.fb_macro_key('worker'),
                                                             interfaces=[IO_READ])
            keys = [k for d in servers_to_keys.values() for k in d]

            # ------------------------------------------------------------------
            nkeys = len(keys)
            LOGGER.debug(f'found {nkeys} keys in ({self.fbpath}) '
                         f'with {self.fbinterface.get_type()} interface')
            if nkeys == 0:
                return

            # ------------------------------------------------------------------
            # compute weight for each frame using patch selector
            frame_weights = np.ones(nkeys, dtype=np.float32)
            if self.pselector is not None:
                # procedure: (1) group different keys (frames) into sims/patches
                #            (2) get (normalized) sim weights from ml
                #            (3) convert sim weights to frame weights

                # step (1): extract the simnames from these keys!
                sim_names = [Naming.fb_macro_key_extract_simname(k) for k in keys]
                sim_names, key2sim = np.unique(np.array(sim_names), return_inverse=True)
                LOGGER.debug(f'sim_names = {sim_names.shape}')

                # step (2): use patch selector to compute sim weights
                # 2021.02.11: moved the locking here
                p = multiprocessing.current_process()
                LOGGER.debug(f'Acquiring lock on patch selector ({p.name}, {p.pid})')
                with lock_patch_selector:
                    LOGGER.debug(f'Acquired lock on patch selector ({p.name}, {p.pid})')
                    _, sim_weights = self.pselector.get_weights(sim_names)
                    LOGGER.debug(f'Released lock on patch selector ({p.name}, {p.pid})')

                # step (3): assign same sim weight to each frame of that sim
                frame_weights = sim_weights[key2sim]

            # ------------------------------------------------------------------
            rdfs = self.fbinterface.load_rdfs_for_macro_feedback(self.fbpath, servers_to_keys, interfaces=[IO_READ])
            LOGGER.debug(f'found {len(rdfs)} rdfs in ({self.fbpath}) '
                         f'with {self.fbinterface.get_type()} interface')

            # ------------------------------------------------------------------
            success_idxs = []
            failed_idxs = []
            for i in range(nkeys):
                if rdfs[i] is not None:
                    success_idxs.append(i)
                else:
                    failed_idxs.append(i)
            success_idxs = np.array(success_idxs)
            if len(failed_idxs) > 0:
                nfailed = len(keys) - len(success_idxs)
                LOGGER.warning(f'Found {nfailed} rdfs=None')
                for i in failed_idxs:
                    failed_key = keys[i]
                    for server in servers_to_keys:
                        while failed_key in servers_to_keys[server]:
                            servers_to_keys[server].remove(failed_key)

                rdfs = [rdfs[i] for i in success_idxs]
                keys = [keys[i] for i in success_idxs]
                frame_weights = frame_weights[success_idxs]

            rdfs = np.array(rdfs)

            # ------------------------------------------------------------------
            # serial: scale and aggregate
            t0 = time.time()
            self.rdfs.scale_and_aggregate_many(rdfs, frame_weights)
            LOGGER.debug(f'aggregated {len(rdfs)} took {time.time() - t0} sec.')

            # ------------------------------------------------------------------
            # TODO: this is incorrect.
            #  we should only delete the keys that were successfully loaded
            self.fbinterface.remove_keys_for_macro_feedback(self.fbpath, servers_to_keys, interfaces=[IO_READ])

        # ----------------------------------------------------------------------
        elif self.type == FeedbackManagerType.Worker:
            assert isinstance(rdfs, np.ndarray)
            assert rdfs.shape == RDF.SHAPE and rdfs.dtype == RDF.DTYPE

            LOGGER.debug(f'Aggregating a new rdf: {self}. sum(rdf) = {np.sum(rdfs)}')
            self.rdfs.aggregate(rdfs)

        # ----------------------------------------------------------------------
        self.naggregates += 1
        LOGGER.debug(f'Aggregated data: {self}')

    # --------------------------------------------------------------------------
    # report the analysis
    def report(self):

        # ----------------------------------------------------------------------
        if self.type == FeedbackManagerType.Worker:

            # Reporting analysis
            key = Naming.fb_macro_key('worker', simname=self.name, nreports=self.nreports)
            LOGGER.info(f'Reporting analysis ({self.name}) fbpath = ({self.fbpath}), '
                        f'key ({key}) using interface ({self.fbinterface.get_type()})')
            self.fbinterface.save_rdfs(self.fbpath, key, self.rdfs, interfaces=['simple', 'redis'])

        # ----------------------------------------------------------------------
        elif self.type == FeedbackManagerType.Manager:

            if self.naggregates == self.naggs_at_last_report:
                LOGGER.info(f'Skip reporting stale analysis: {self}')
                return

            key = Naming.fb_macro_key('manager', nreports=self.nreports)
            LOGGER.info(f'Reporting analysis {self.name} '
                        f'fbpath = ({self.fbpath}), key ({key})')

            # collect the history in tar
            key = key + '-' + time.strftime("%Y%m%d-%H%M%S")
            self.iointerface.save_rdfs(self.fname_master, key, self.rdfs, interfaces=['taridx'])

            # write the latest one as npz
            key = 'rdf_main.npz'
            self.fbinterface.save_rdfs(self.fbpath, key, self.rdfs)

            # write the latest ones in text+npz format for macro
            RDF.write_rdfs_individually(self.rdfs.data, self.outpath,
                                        mummi_core.get_io('simple'))

        # ----------------------------------------------------------------------
        self.naggs_at_last_report = self.naggregates
        self.nreports += 1
        LOGGER.info('Reported analysis: {}'.format(self.__str__()))

    # --------------------------------------------------------------------------
    def checkpoint(self):

        _type = 'm' if FeedbackManagerType.Manager else 'w'

        LOGGER.info('Checkpointing: {}'.format(self.__str__()))
        state = {'naggregates': self.naggregates, 'nreports': self.nreports,
                 'naggs_at_last_report': self.naggs_at_last_report,
                 'type': _type, 'name': self.name,
                 'rdf': {'naggregate': self.rdfs.naggregate}
                 }

        self.iointerface.save_checkpoint(os.path.join(self.workspace, self.fname_chkpt), state)

        self.iointerface.take_backup(os.path.join(self.workspace, self.fname_chkpt_data))
        self.iointerface.save_npz(self.workspace, self.fname_chkpt_data, self.rdfs,
                                  writer_func=RDF.save_npz, interfaces=['simple'])

    def restore(self):

        #state = self.iointerface.load_checkpoint(os.path.join(self.workspace, self.fname_chkpt))

        filename = os.path.join(self.workspace, self.fname_chkpt)
        if not os.path.isfile(filename):
            LOGGER.info('Checkpoint file {} does not exist!'.format(filename))
            return

        with open(filename, 'r') as infile:
            data = yaml.load(infile, Loader=yaml.Loader)

        LOGGER.info('Restored checkpoint file {} from {}'.format(filename, data['ts']))
        state = data
        # -------

        if len(state) == 0:
            return

        self.nreports = int(state['nreports'])
        self.naggregates = int(state['naggregates'])
        self.naggs_at_last_report = int(state['naggs_at_last_report'])

        if self.type == FeedbackManagerType.Manager:
            assert 'm' == state['type']
        else:
            assert 'w' == state['type']
        assert self.nreports >= 0
        assert self.naggregates >= 0
        assert self.naggs_at_last_report >= 0
        assert self.naggs_at_last_report <= self.naggregates

        self.rdfs = self.iointerface.load_npz(self.workspace,
                                              self.fname_chkpt_data,
                                              reader_func=RDF.load_npz,
                                              interfaces=['simple'])
        assert self.rdfs.naggregate == state['rdf']['naggregate']
        LOGGER.info('Restored: {}'.format(self.__str__()))

    def test(self):
        raise NotImplementedError('todo')

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
