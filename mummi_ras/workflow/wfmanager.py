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
import signal
import sys
import yaml
import timeit
import traceback
import multiprocessing
from multiprocessing.managers import SyncManager
import pika

import mummi_core
from mummi_core.workflow.job import JOB_TYPES, JOB_NEXT_QUEUE
from mummi_core.workflow.flux_env import flux_uri
from mummi_core.workflow.jobTracker import JobTracker
from mummi_core.utils.timer import Timer
from mummi_core.brokers import RPCClient

import mummi_ras
from mummi_ras import Naming
from mummi_ras.transformations.patch_creator import MacroPatchCreator
from mummi_ras.ml.selectors import CGSelector, CGSelectorType
from mummi_ras.feedback import MacroFeedbackManager, FeedbackManagerType, FeedbackManager_AA2CG

from logging import getLogger
LOGGER = getLogger(__name__)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
class ProxyManager(SyncManager): pass


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
class WFManager():

    # --------------------------------------------------------------------------
    def __init__(self, job_name):#, exec_bucket):

        self.job_name = job_name
        # self.exception_queue = exec_bucket
        self.processes = {}

        self._exit = multiprocessing.Event()     # Trigger the daemon to exit.
        self._error = multiprocessing.Event()    # Trigger the daemon to exit.
        self._setup = multiprocessing.Event()    # Wait until the daemon is setup.

        # use this event to manage when wf is ready for
        # the other processes to take over
        # currently, we set it after the first wf loop has finished
        self._wf_ready = multiprocessing.Event()

    # --------------------------------------------------------------------------
    def is_exit(self):
        return self._exit.is_set()

    def is_error(self):
        return self._error.is_set()

    def exit(self):
        self._exit.set()

    def error(self):
        self._error.set()

    # --------------------------------------------------------------------------
    def signal_wrapper(self, name, pid):
        def handler(signum, frame):
            LOGGER.warning(f'Received (SIGNUM={signum}) for "{name}[pid={pid}]"'
                           f'frame "{frame.f_back}"')
            self.checkpoint()
            traceback.print_exc()
            #self.exception_queue.put(sys.exc_info())
            #self.error()  # copied here but not sure if we want this
            #self.exit()
            self.stop()

        return handler

    # --------------------------------------------------------------------------
    def cleanup(self):
        pass

    # --------------------------------------------------------------------------
    def setup(self, env, config):

        if self._setup.is_set():
            return

        mummi_core.init()
        try:
            # ------------------------------------------------------------------
            # read the config
            self.config = config
            self.wconfig = config['wfmanager']['config']

            # ------------------------------------------------------------------
            # the different run modes
            self.do_gc = bool(self.wconfig['is_gc'])
            assert self.do_gc == bool(self.config['macro_patch_creator']['params']['is_gc'])
            self.config['createsim']['config']['is_gc'] = self.do_gc

            self.do_patchcreation = bool(self.wconfig['do_patchcreation'])
            self.do_workflow = bool(self.wconfig['do_workflow'])
            self.do_schedulejobs = bool(self.wconfig['do_schedulejobs'])
            self.do_mlserver = bool(self.wconfig['do_mlserver'])
            self.do_ucgserver = bool(self.wconfig['do_ucgserver'])
            assert self.do_mlserver ^ self.do_ucgserver, "Currently, we do not support ML and UCG server mode at the same time, please chose one"

            self.do_patchselection = bool(self.wconfig['do_patchselection'])
            self.do_feedback_cg2mc = bool(self.wconfig['do_feedback_cg2mc'])

            self.do_cgselection = bool(self.wconfig['do_cgselection'])
            self.do_feedback_aa2cg = bool(self.wconfig['do_feedback_aa2cg'])

            # ------------------------------------------------------------------
            # control the work to be done per iteration
            self.nMaxJobsPerIter = int(self.wconfig['nMaxJobsPerIter'])
            self.nReadPatchesPerIter = int(self.wconfig['nReadPatchesPerIter'])
            self.nReadCGFramesPerIter = int(self.wconfig['nReadCGFramesPerIter'])
            self.nMaxPatchesSelectionsPerIter = int(self.wconfig['nMaxPatchesSelectionsPerIter'])
            self.nMaxCGFramesSelectionsPerIter = int(self.wconfig['nMaxCGFramesSelectionsPerIter'])

            self.nMaxSelectedPatchBuffer = int(self.wconfig['nMaxSelectedPatchBuffer'])
            self.nMaxSelectedCGFrameBuffer = int(self.wconfig['nMaxSelectedCGFrameBuffer'])
            LOGGER.info(f'buffers = {self.nMaxSelectedPatchBuffer}, {self.nMaxSelectedCGFrameBuffer}')

            self.nSecPerIterWF = int(self.wconfig['nSecPerIterWF'])
            self.nSecPerIterPC = int(self.wconfig['nSecPerIterPC'])

            self.nSecPerCgSelUpdate = int(self.wconfig['nSecPerCgSelUpdate'])

            self.nSecPerMacroFeedback = int(self.wconfig['nSecPerMacroFeedback'])
            self.nSecPerCGFeedback = int(self.wconfig['nSecPerCGFeedback'])
            self.iterCounterMLServer = int(self.wconfig['mlserver_round_id'])
            self.iterCounterUCGServer = int(self.wconfig['ucgserver_round_id'])
            self.generated_patches = [] # List of patches submitted to createsim

            # ------------------------------------------------------------------
            # identify scheduler interface
            self.hostname = mummi_core.get_hostname()
            LOGGER.info(f'Workflow launched on host: {self.hostname}')

            self.flux = flux_uri(override_uri_from_file=False)
            LOGGER.info(f'  Flux  URI: [{self.flux}]')
            LOGGER.info(f'  Local URI: [{os.environ.get("FLUX_URI", None)}]')

            if self.flux is None:
                LOGGER.error('Failed to identify flux URI. Disabling job scheduling')
                self.do_schedulejobs = False
                assert False, 'Failed to identify flux URI'

            # ------------------------------------------------------------------
            # if no workspace has been specified, use the default in hierarchy
            self.workspace = self.wconfig.get('workspace', '')
            if self.workspace == '':
                self.workspace = Naming.dir_root('workspace')

            LOGGER.info('workspace = {}'.format(self.workspace))

            # checkpoint will be in workspace
            self.chkpt = self.wconfig['jobname'] + '.chk'
            self.chkpt = os.path.join(self.workspace, self.chkpt)

            # ------------------------------------------------------------------
            # iointerfaces
            self.iointerface = mummi_ras.get_io(self.wconfig['iotype'])
            self.fbinterface = mummi_ras.get_io(self.wconfig['fbtype'])
            LOGGER.info('Using [%s] for io', self.iointerface)
            LOGGER.info('Using [%s] for feedback', self.fbinterface)

            self.ns_pfpatches = Naming.dir_root('patches')

            # ------------------------------------------------------------------
            # workflow state
            # ------------------------------------------------------------------
            self.iterCounterWF = 0
            self.iterCounterPC = 0
            self.iterCounterFB_cg2macro = 0
            self.iterCounterFB_aa2cg = 0

            self.patchCounter = 0
            self.cgPatchCounter = 0
            self.cgCurrentFrameIdx = 0  # Add CG time frames incrementally

            self.timeOfLastCgSelUpdate = 0
            # ------------------------------------------------------------------

        except Exception as e:
            # self.exception_queue.put(sys.exc_info())
            traceback.print_exc()
            LOGGER.error('WF failed during setup!')
            self.error()
            raise e

        # Set that the deamon is now setup.
        self._setup.set()

    # ------------------------------------------------------------------------
    def _init_scheduling(self):

        # Maestro Adapter Settings
        self.adapter_batch = self.config['wfmanager']['batch']

        # ----------------------------------------------------------------------
        # get the nodes and split them for different jobs
        nnodes_tot = int(os.environ['MUMMI_NNODES'])
        nnodes_macro = int(os.environ.get('MUMMI_MACRO_NNODES', 0))
        nnodes_redis = int(os.environ.get('MUMMI_REDIS_NNODES', 0))
        nnodes_mlserver = int(os.environ.get('MUMMI_MLSERVER_NNODES', 0))
        nnodes_ucgserver = int(os.environ.get('MUMMI_UCGSERVER_NNODES', 0))
        nnodes_feedbacks = int(os.environ.get('MUMMI_FB_NNODES', 0))

        if self.do_mlserver and nnodes_mlserver == 0:
            LOGGER.error('ML Server is supposed to be running but MUMMI_MLSERVER_NNODES is set to 0.')
        if self.do_ucgserver and nnodes_ucgserver == 0:
            LOGGER.error('UCG Server is supposed to be running but MUMMI_UCGSERVER_NNODES is set to 0.')

        portion_csim = float(self.wconfig['portion_csim'])
        portion_cg = float(self.wconfig['portion_cg'])
        assert 0 <= portion_csim <= 1
        assert 0 <= portion_cg <= 1

        # Fixed cost Infrastructure nodes
        nnodes_flux = 0   # Dedicated full node for Flux
        nnodes_wf = 1     # Dedicated CPU-only for wfmanager

        # Node computations
        nnodes_tot = nnodes_tot - nnodes_flux
        nnodes_mlserver = 0 # ML server and wfmanger are co-located but ML server needs one GPU

        # For GROMACS we move to vertical scheduling (cganalysis alone on the node and createsims alone, no sharing anymore)
        nnodes_tot = nnodes_tot - (nnodes_wf + nnodes_feedbacks + nnodes_macro + nnodes_redis + nnodes_mlserver + nnodes_ucgserver)
        nnodes_gpu = nnodes_tot
        nnodes_cpu = nnodes_tot # createsims uses GPU in that version
        # For GROMACS we cannot have more than 1 as createsims and cganalysis use GPUs (4 per nodes)
        self.nnodes_csim : int = max(0, int(nnodes_gpu * portion_csim))
        self.nnodes_cg : int = max(0, nnodes_gpu - self.nnodes_csim)
        # END For GROMACS section

        # Classic ddcMD section 
        # nnodes_gpu = nnodes_tot
        # nnodes_cpu = nnodes_tot - (nnodes_wf + nnodes_macro + nnodes_redis + nnodes_mlserver + nnodes_ucgserver)

        # self.nnodes_csim : int = int(nnodes_cpu * portion_csim)
        # self.nnodes_cg : int = int(nnodes_gpu * portion_cg)

        # self.nnodes_bmap : int = max(0, nnodes_cpu - self.nnodes_csim)
        # self.nnodes_aa : int = max(0, nnodes_gpu - self.nnodes_cg)

        # WARNING: we are not using these jobs at the moment
        self.nnodes_bmap : int = 0
        self.nnodes_aa : int = 0

        LOGGER.info(f'NNodes: total = {nnodes_tot}: '
                    f'(wf = {nnodes_wf}, redis = {nnodes_redis},'
                    f' macro = {nnodes_macro}, mlserver = {nnodes_mlserver}'
                    f' ucgserver = {nnodes_ucgserver}, feeedback = {nnodes_feedbacks})')

        LOGGER.info(f'NNodes: total = {nnodes_tot}: '
                    f'(csim = {self.nnodes_csim}, cg = {self.nnodes_cg},'
                    f' bmap = {self.nnodes_bmap}, aa = {self.nnodes_aa})')

        # ----------------------------------------------------------------------
        self.job_trackers = {
            'createsim': JobTracker(self.config['createsim'],
                            self.nnodes_csim, self.iointerface,
                            self.adapter_batch,
                            self.workspace,
                            self.do_schedulejobs),
                            # self.exception_queue, self.do_schedulejobs),

            'cg': JobTracker(self.config['cg'],
                            self.nnodes_cg, self.iointerface,
                            self.adapter_batch,
                            self.workspace,
                            self.do_schedulejobs),
                            # self.exception_queue, self.do_schedulejobs),

            'backmapping': JobTracker(self.config['backmapping'],
                            self.nnodes_bmap, self.iointerface,
                            self.adapter_batch,
                            self.workspace,
                            self.do_schedulejobs),
                            # self.exception_queue, self.do_schedulejobs),

            'aa': JobTracker(self.config['aa'],
                            self.nnodes_aa, self.iointerface,
                            self.adapter_batch,
                            self.workspace,
                            self.do_schedulejobs)
                            # self.exception_queue, self.do_schedulejobs)}
        }

    def _init_recovery(self):

        return
        LOGGER.error('recovery is a deprecated option!\n')
        exit(1)

        if self.do_recovery == '':
            self.do_recovery = False
            return

        # read the data from the recovery file
        LOGGER.info('Starting Recovery mode using {}'.format(self.do_recovery))

        with open(self.do_recovery) as f:
            content = f.readlines()
            content = [x.strip() for x in content]

        recovery_patches = [p.split()[0] for p in content]

        # queue these patches for ddcmd
        self.job_trackers['cg'].add_to_queue(recovery_patches)
        LOGGER.info('Queued {} patches for CG'.format(len(recovery_patches)))

        # and change the name of the chkpoint file
        self.chkpt = self.wconfig['jobname'] + '.recovery.chk'
        self.chkpt = os.path.join(self.workspace, self.chkpt)
        LOGGER.info('Changed checkpointing to ({})'.format(self.chkpt))

        # recovery mode overrides all settings
        self.nMaxJobsPerIter = 24
        self.do_recovery = True
        self.do_workflow = True
        self.do_schedulejobs = True
        #self.do_cold_restart = True
        self.do_patchcreation = False
        self.do_patchselection = False
        self.do_cgselection = False

    def _init_patch_creation(self):

        self.patch_creator = None
        if not self.do_patchcreation:
            return

        self.patch_creator = MacroPatchCreator(self.config["macro_patch_creator"],
                                               self.iointerface, self.iointerface)
        self.patch_creator.restore()

    # def _init_patch_selection(self):

    #     self.pselector = None
    #     if self.do_patchselection or self.do_feedback_cg2mc:
    #         ml_config = Naming.ml('macro')
    #         stype = ml_config['selection']
    #         assert stype in ['importance', 'random']

    #         wspace = os.path.join(Naming.dir_root('ml'), 'macro')
    #         self.pselector = self.proxy_manager.PatchSelector(stype, wspace, ml_config)
    #         self.pselector.restore()

    def _init_cg_selection(self):

        self.cgselector = None

        if self.do_cgselection:
            ml_config = Naming.ml('cg')
            stype = ml_config['selection']
            assert stype == 'importance'

            wspace = os.path.join(Naming.dir_root('ml'), 'cg')
            os.makedirs(wspace, exist_ok=True)
            self.cgselector = CGSelector(CGSelectorType.Manager, 'manager',
                                         wspace, ml_config)
            self.cgselector.restore()

    def _init_feedback_cg2macro(self):

        self.fb_cg2macro = None

        # ----------------------------------------------------------------------
        if self.do_feedback_cg2mc:
            name = 'CG2MacroFeedback_Manager'
            wspace = Naming.dir_root('feedback-cg')

            # should ideally use io interface (here and below)
            #outpath = self.iointerface.namespace('feedback-cg')
            #fbpath = self.fbinterface.namespace('feedback-cg')
            outpath = wspace
            fbpath = wspace

            do_weights = bool(self.wconfig['fbcg_do_wts'])
            pselector = self.pselector if do_weights else None

            # fdinatal -- removal of dbr due to hang
            # self.fbinterface.create_namespace(fbpath)
            self.fb_cg2macro = MacroFeedbackManager(FeedbackManagerType.Manager,
                                                    name, wspace, outpath, fbpath,
                                                    self.iointerface, self.fbinterface,
                                                    pselector)
            self.fb_cg2macro.restore()

        # ----------------------------------------------------------------------

    def _init_feedback_aa2cg(self):

        self.fb_aa2cg = None

        # ----------------------------------------------------------------------
        if self.do_feedback_aa2cg:
            wspace = Naming.dir_root('feedback-aa')

            #outpath = self.iointerface.namespace('feedback-aa')
            #fbpath = self.fbinterface.namespace('feedback-aa')
            outpath = wspace
            fbpath = wspace

            #hvr_th = 0.2
            #crd_th = 0.15
            hvr_th = float(self.wconfig['fbaa_hvr_th'])
            crd_th = float(self.wconfig['fbaa_crd_th'])
            frm_incrmnt = int(self.wconfig['fbaa_frame_increment'])

            # fdinatal -- removal of dbr due to hang
            # self.fbinterface.create_namespace(fbpath)
            self.fb_aa2cg = FeedbackManager_AA2CG(FeedbackManagerType.Manager,
                                                  wspace, outpath, fbpath,
                                                  self.iointerface, self.fbinterface,
                                                  hvr_th, crd_th, frm_incrmnt)
            self.fb_aa2cg.restore()

        # ----------------------------------------------------------------------

    def _init_ml_server(self):
        """
        Initialize the client on the workflow side to talk to the ML Server.
        March 2023: 
            Current implementation supports only RabbitMQ, currently the ML server connect to
            a RabbitMQ server running on a Kubernetes cluster (defined in the credentials).
            The ML server is responsible of sampling, generating and validating new structures
            for the wfmanager.
        """
        if self.do_mlserver:
            if Naming.MUMMI_MINI:
                spath = os.path.join(mummi_ras.Naming.MUMMI_SPECS, "mini-mummi")
            else:
                spath = os.path.join(mummi_ras.Naming.MUMMI_SPECS, "workflow")

            config = mummi_core.read_specs(os.path.join(spath, "mlserver.yaml"))

            wspace = config["workspace"]["path"]
            credentials_path = os.path.join(wspace,
                config["workspace"]["credentials"]
            )   
            certificate_path = os.path.join(wspace,
                config["workspace"]["certificate"]
            )

            usr = os.environ.get("USER", "mummiusr")
            routing_key = config["broker"]["queue"] + '_' + usr
            if usr == "mummiusr":
                LOGGER.warning(f"Did not find current user: defaulted to {usr}")
            broker_interface = config["broker"]["interface"]
            try:
                # interface = get_broker(broker_interface)
                self.rpc_client = RPCClient(
                    credentials = credentials_path,
                    cacert = certificate_path,
                    routing_key = routing_key,
                )
            except Exception as e:
                traceback.print_exc()
                LOGGER.error("ML Server failed during setup!")
                self.error()
                raise e

            LOGGER.info(f"> Initializing RPC client to interact with ML Server")
            LOGGER.info(f"  > Client")
            LOGGER.info(f"    > Interface   {broker_interface}")
            LOGGER.info(f"    > Credentials {credentials_path}")
            LOGGER.info(f"    > Certificate {certificate_path}")
            LOGGER.info(f"    > Routing key {routing_key}")

    def _init_ucg_server(self):
        """
        Initialize the client on the workflow side to talk to the UCG Server.
        June 2023: 
            Current implementation supports only RabbitMQ, currently the UCG server connect to
            a RabbitMQ server running on a Kubernetes cluster (defined in the credentials).
            The UCG server is responsible of reading trajectories, selecting, generating 
            and validating new structures for the wfmanager.
        """
        if self.do_ucgserver:
            if Naming.MUMMI_MINI:
                spath = os.path.join(mummi_ras.Naming.MUMMI_SPECS, "mini-mummi")
            else:
                spath = os.path.join(mummi_ras.Naming.MUMMI_SPECS, "workflow")

            config = mummi_core.read_specs(os.path.join(spath, "ucgserver.yaml"))


            wspace = config["workspace"]["path"]
            credentials_path = os.path.join(wspace,
                config["workspace"]["credentials"]
            )
            certificate_path = os.path.join(wspace,
                config["workspace"]["certificate"]
            )

            usr = os.environ.get("USER", "mummiusr")
            routing_key = config["broker"]["queue"] + '_' + usr
            if usr == "mummiusr":
                LOGGER.warning(f"Did not find current user: defaulted to {usr}")
            broker_interface = config["broker"]["interface"]
            try:
                self.rpc_client = RPCClient(
                    credentials = credentials_path,
                    cacert = certificate_path,
                    routing_key = routing_key,
                )
            except Exception as e:
                traceback.print_exc()
                LOGGER.error("UCG Server failed during setup!")
                self.error()
                raise e

            LOGGER.info(f"> Initializing RPC client to interact with UCG Server")
            LOGGER.info(f"  > Client")
            LOGGER.info(f"    > Interface   {broker_interface}")
            LOGGER.info(f"    > Credentials {credentials_path}")
            LOGGER.info(f"    > Certificate {certificate_path}")
            LOGGER.info(f"    > Routing key {routing_key}")

    # --------------------------------------------------------------------------
    def checkpoint(self):

        state = dict(
            flux = self.flux,
            iterCounterWF = self.iterCounterWF,
            iterCounterPC = self.iterCounterPC,
            patchCounter = self.patchCounter,
            jobs_createsim = self.job_trackers['createsim'].status(),
            jobs_backmapping = self.job_trackers['backmapping'].status(),
            jobs_cg = self.job_trackers['cg'].status(),
            jobs_aa = self.job_trackers['aa'].status())

        # We keep the two last checkpoints
        self.iointerface.save_checkpoint(self.chkpt,
            state,
            use_tstamp=True,
            cleanup=True,
            keep_checkpoint=2
        )

    def restore(self):

        state = self.iointerface.load_checkpoint(self.chkpt, loader=yaml.UnsafeLoader)
        if len(state) == 0:
            return

        LOGGER.info(f"Restoring workflow as of {state['ts']}")
        sys.stdout.flush()

        # restore the data
        self.iterCounterWF = state['iterCounterWF']
        self.iterCounterPC = state['iterCounterPC']
        self.patchCounter = state['patchCounter']

        prev_flux = state['flux']

        warm_restart = self.flux == prev_flux
        wstring = 'Warm' if warm_restart else 'Cold'
        LOGGER.info(f'{wstring} restart of the Workflow! flx={self.flux}, prev_flx={prev_flux}')

        sims_success = {}
        sims_failed = {}
        for j in JOB_TYPES:
            job_state = state['jobs_' + j]
            sims_success[j], sims_failed[j] = self.job_trackers[j].restore(job_state, warm_restart)
            LOGGER.info(self.job_trackers[j].__str__())

        # if any jobs were found successful finished, need to queue for the next step!
        for j in JOB_NEXT_QUEUE.keys():
            fjobs = sims_success[j]
            if len(fjobs) > 0:
                LOGGER.info(f'Found {len(fjobs)} successful {j} sims!')
                self.job_trackers[JOB_NEXT_QUEUE[j]].add_to_queue(fjobs)

        LOGGER.info(f"Restored WFManager from {state['ts']}. "
                f"iterCounterWF = {self.iterCounterWF}, "
                f"iterCounterPC = {self.iterCounterPC}, "
                f"patchCounter = {self.patchCounter}")

    # ------------------------------------------------------------------------
    # Main tasks to be done by the wf manager
    # ------------------------------------------------------------------------
    def _task_add_new_patches_to_ml(self, lock_patch_io, lock_patch_select):

        if not self.do_patchselection:
            LOGGER.debug(f"No patchselection as do_patchselection = {self.do_patchselection}")
            return 0
        LOGGER.debug(f"Patch selection started")

        p = multiprocessing.current_process()

        # ----------------------------------------------------------------------
        patches = []
        if self.do_gc:
            patches = self.patch_creator.run_for_gc()

        else:
            patch_ids = [Naming.pfpatch(self.patchCounter+i)
                         for i in range(self.nReadPatchesPerIter)]

            if not lock_patch_io.acquire(False):  # non-blocking acquire
                LOGGER.debug("Failed to acquire lock on patches ({}, {})".format(p.name, p.pid))
                return 0

            LOGGER.debug("Acquired lock on patches ({}, {})".format(p.name, p.pid))
            patches = self.iointerface.load_patches(self.ns_pfpatches, patch_ids)

            lock_patch_io.release()
            LOGGER.debug("Released lock on patches ({}, {})".format(p.name, p.pid))

        # ----------------------------------------------------------------------
        npatches = len(patches)
        if npatches == 0:
            return npatches

        LOGGER.debug("Acquiring lock on patch selector ({}, {})".format(p.name, p.pid))
        with lock_patch_select:
            LOGGER.debug("Acquired lock on patch selector ({}, {})".format(p.name, p.pid))
            self.pselector.add_candidates(patches)
        LOGGER.debug("Released lock on patch selector ({}, {})".format(p.name, p.pid))

        self.patch_creator.checkpoint()

        self.patchCounter += npatches
        LOGGER.info('Added {} patches to the selector'.format(npatches))
        return npatches

    def _task_add_cgframes_to_ml(self):

        if not self.do_cgselection:
            LOGGER.debug(f"No cgselection as do_cgselection = {self.do_patchselection}")
            return 0
        LOGGER.debug(f"CG selection started")

        ctime = timeit.default_timer() - self.timeOfLastCgSelUpdate
        if ctime < self.nSecPerCgSelUpdate:
            return

        try:
            self.cgselector.update_manager()
        except Exception as e:
            LOGGER.error(f'Got error from cgselector.update_manager: {e}')
        self.timeOfLastCgSelUpdate = timeit.default_timer()

    # --------------------------------------------------------------------------
    def _task_select_pfpatches(self, lock_patch_select):

        if not self.do_patchselection:
            return 0

        nPendingPatches = self.job_trackers['cg'].nqueued_sims() + \
                          self.job_trackers['createsim'].nqueued_sims() + \
                          self.job_trackers['createsim'].nrunning_sims()

        LOGGER.debug(f"pending patches = {nPendingPatches} = {self.job_trackers['cg'].nqueued_sims()} "
                     f"+ {self.job_trackers['createsim'].nqueued_sims()} "
                     f"+ {self.job_trackers['createsim'].nrunning_sims()}")
        nPatches = min(self.nMaxSelectedPatchBuffer - nPendingPatches,
                       self.nMaxPatchesSelectionsPerIter)

        LOGGER.debug(f"npatches = {nPatches} = min({self.nMaxSelectedPatchBuffer} "
                     f"- {nPendingPatches}, {self.nMaxPatchesSelectionsPerIter}")

        if nPatches <= 0:
            return 0

        # select new patches
        p = multiprocessing.current_process()
        LOGGER.debug("Acquiring lock on patch selector ({}, {})".format(p.name, p.pid))
        with lock_patch_select:
            LOGGER.debug("Acquired lock on patch selector ({}, {})".format(p.name, p.pid))
            LOGGER.info('Select {} new patches'.format(nPatches))
            selections = self.pselector.select(nPatches)
            LOGGER.debug("Released lock on patch selector ({}, {})".format(p.name, p.pid))

        n = len(selections)
        if n == 0:
            return n

        # test these selections
        # Split a list of sims based on their status.
        # Returns: sims_success, sims_failed, sims_unknown
        _s, _f, _u, _stop = self.job_trackers['createsim'].split_sims_on_status(selections)

        # this is the correct scenario
        if len(_s) == 0 and len(_f) == 0:
            self.job_trackers['createsim'].add_to_queue(selections)

        # looks like some of these selections were already simulated
        else:
            LOGGER.error('Found (%d %d) selected patches for which createsim was already run. '
                         'Were you using old sims with new ml?', len(_s), len(_f))
            self.job_trackers['createsim'].add_to_queue(_u)
            self.job_trackers['cg'].add_to_queue(_s)

        return n

    def _request_ml_structures(self, nPatches: int):
        selections = self.rpc_client.call({'k_samples': nPatches, 'iteration_id': self.iterCounterMLServer, 'strict': False})
        if isinstance(selections, str):
            LOGGER.error(f"[.] RPC server sent back an exception: {selections}")
            self.checkpoint()
            return []
        selections = [x.split('/')[-1] for x in selections]
        LOGGER.info(f" [.] Got generate_new_samples(iter={self.iterCounterMLServer}, k={nPatches}) = {selections}")
        for sample in selections:
            if sample in self.generated_patches:
                LOGGER.error(f"Sample {sample} has already been generated. This will cause trouble. {sample} should be removed.")
        self.generated_patches += selections
        return selections

    # --------------------------------------------------------------------------
    # Instead of using PatchSelector we query the ML Server for structures to
    # to give to createsim
    # --------------------------------------------------------------------------
    def _task_add_new_patches_mlserver(self):
        if not self.do_mlserver:
            return 0

        # Running and Queued CG represent a certain number of createsims that have been executed in the past
        # Running Createsims represent the current createsim running
        # Queued Createsims represent the future createsim that will be run at some point.
        nPendingPatches = self.job_trackers['cg'].nqueued_sims() + \
                          self.job_trackers['createsim'].nqueued_sims() + \
                          self.job_trackers['createsim'].nrunning_sims()

        LOGGER.info(f"pending patches = {nPendingPatches} = {self.job_trackers['cg'].nqueued_sims()} "
                    f"+ {self.job_trackers['createsim'].nqueued_sims()} "
                    f"+ {self.job_trackers['createsim'].nrunning_sims()}")
        nPatches = min(self.nMaxSelectedPatchBuffer - nPendingPatches, self.nMaxPatchesSelectionsPerIter)

        LOGGER.info(f"npatches => min({self.nMaxSelectedPatchBuffer} - {nPendingPatches}, "
                    f"{self.nMaxPatchesSelectionsPerIter}) = {nPatches}")

        # We have enough patches
        if nPatches <= 0:
            return 0

        # if nPatches is too small and the ML model success rate is low (<40%)
        # Then requesting 5 or 10 will never return valid patches statistically
        if nPatches < self.nMaxPatchesSelectionsPerIter:
            nPatches = self.nMaxPatchesSelectionsPerIter

        # select new patches
        p = multiprocessing.current_process()
        LOGGER.debug(f"> Running RPC client {p.name} on [{self.hostname}/{p.pid}]")
        LOGGER.info(f" [x] Requesting generate_new_samples(iter={self.iterCounterMLServer}, k={nPatches})")
        # If strict equal True then call() will ALWAYS return the requested number of VALID patches (but it will take longer to run)
        selections = []
        try:
            # selections = self.rpc_client.call({'k_samples': nPatches, 'iteration_id': self.iterCounterMLServer, 'strict': False})
            # if isinstance(selections, str):
            #     LOGGER.error(f"[.] RPC server sent back an exception: {selections}")
            #     self.checkpoint()
            #     # self.stop()
            #     return -1
            # selections = [x.split('/')[-1] for x in selections]
            # LOGGER.info(f" [.] Got generate_new_samples(iter={self.iterCounterMLServer}, k={nPatches}) = {selections}")
            # for sample in selections:
            #     if sample in self.generated_patches:
            #         LOGGER.error(f"Sample {sample} has already been generated. This will cause trouble. {sample} should be removed.")
            # self.generated_patches += selections
            selections = self._request_ml_structures(nPatches)
        except pika.exceptions.AMQPError as e:
            # we try to restart the connection
            LOGGER.warning(f"We got {e}")
            LOGGER.info(f"We try to reinitialize the connection with ML server")
            self._init_ml_server()
            LOGGER.info(f"Connection to ML server reinitialized.")
            # Retrying
            LOGGER.info(f" [x] Requesting generate_new_samples(iter={self.iterCounterMLServer}, k={nPatches})")
            selections = self._request_ml_structures(nPatches)

        n = len(selections)
        if n == 0:
            return n

        # test these selections
        # Split a list of sims based on their status.
        # Returns: sims_success, sims_failed, sims_unknown, sims_stop
        sims_success, sims_failed, sims_unknown, sims_stop = self.job_trackers['createsim'].split_sims_on_status(selections)

        LOGGER.debug(f"Found sims_success = {sims_success}, sims_failed = {sims_failed}, sims_unknown= {sims_unknown}, sims_unknown= {sims_stop}")

        # this is the correct scenario (the new patches have not been explored yet)
        if len(sims_success) == 0 and len(sims_failed) == 0:
            self.job_trackers['createsim'].add_to_queue(selections)

        # looks like some of these selections were already simulated
        else:
            LOGGER.error('Found (%d %d) selected patches for which createsim was already run. '
                         'Were you using old sims with new ML?', len(sims_success), len(sims_failed))
            self.job_trackers['createsim'].add_to_queue(sims_unknown)
            self.job_trackers['cg'].add_to_queue(sims_success)

        return n

    def _task_add_new_patches_ucgserver(self):
        if not self.do_ucgserver:
            return 0

        # Running and Queued CG represent a certain number of createsims that have been executed in the past
        # Running Createsims represent the current createsim running
        # Queued Createsims represent the future createsim that will be run at some point.
        nPendingPatches = self.job_trackers['cg'].nqueued_sims() + \
                          self.job_trackers['createsim'].nqueued_sims() + \
                          self.job_trackers['createsim'].nrunning_sims()

        LOGGER.info(f"pending patches = {nPendingPatches} = {self.job_trackers['cg'].nqueued_sims()} "
                    f"+ {self.job_trackers['createsim'].nqueued_sims()} "
                    f"+ {self.job_trackers['createsim'].nrunning_sims()}")
        nPatches = min(self.nMaxSelectedPatchBuffer - nPendingPatches, self.nMaxPatchesSelectionsPerIter)

        LOGGER.info(f"npatches => min({self.nMaxSelectedPatchBuffer} - {nPendingPatches}, "
                    f"{self.nMaxPatchesSelectionsPerIter} = {nPatches})")

        if nPatches <= 0:
            return 0

        # select new patches
        p = multiprocessing.current_process()
        LOGGER.debug(f"> Running RPC client {p.name} on [{self.hostname}/{p.pid}]")
        LOGGER.info(f" [x] Requesting generate_new_samples(iter={self.iterCounterUCGServer}, k={nPatches})")
        # If strict equal True then call() will always return the requested number of VALID patches (but it will take longer to run)
        selections = []
        try:
            selections = self.rpc_client.call({'k_samples': nPatches, 'iteration_id': self.iterCounterUCGServer, 'strict': False})
            if isinstance(selections, str):
                LOGGER.error(f"[.] RPC server sent back an exception: {selections}")
                self.checkpoint()
                # self.stop()
                return -1
            LOGGER.info(f" [.] Got generate_new_samples(iter={self.iterCounterUCGServer}, k={nPatches}) = {selections}")
            selections = [x.split('/')[-1] for x in selections]

            for sample in selections:
                if sample in self.generated_patches:
                    LOGGER.warning(f"Sample {sample} has already been generated. This will cause trouble. {sample} should be removed.")
            self.generated_patches += selections
        except pika.exceptions.AMQPError as e:
            # we try to restart the connection
            LOGGER.warning(f"We got {e}")
            LOGGER.info(f"We try to reinitialize the connection with UCG server")
            self._init_ucg_server()
            LOGGER.info(f"Connection to UCG server reinitialized.")
            return 0

        n = len(selections)
        if n == 0:
            return n

        # test these selections
        # Split a list of sims based on their status.
        # Returns: sims_success, sims_failed, sims_unknown
        sims_success, sims_failed, sims_unknown, _ = self.job_trackers['createsim'].split_sims_on_status(selections)

        LOGGER.debug(f"Found sims_success = {sims_success}, sims_failed = {sims_failed}, sims_unknown= {sims_unknown}")

        # this is the correct scenario (the new patches have not been explored yet)
        if len(sims_success) == 0 and len(sims_failed) == 0:
            self.job_trackers['createsim'].add_to_queue(selections)

        # looks like some of these selections were already simulated
        else:
            LOGGER.error('Found (%d %d) selected patches for which createsim was already run. '
                         'Were you using old sims with new ML?', len(sims_success), len(sims_failed))
            self.job_trackers['createsim'].add_to_queue(sims_unknown)
            self.job_trackers['cg'].add_to_queue(sims_success)

        return n

    # --------------------------------------------------------------------------
    def _task_select_cgframes(self):

        if not self.do_cgselection:
            return 0

        nPendingFrames = self.job_trackers['aa'].nqueued_sims() + \
                          self.job_trackers['backmapping'].nqueued_sims() + \
                          self.job_trackers['backmapping'].nrunning_sims()

        LOGGER.debug(f"pending frames = {nPendingFrames} = {self.job_trackers['aa'].nqueued_sims()} "
                     f"+ {self.job_trackers['backmapping'].nqueued_sims()} "
                     f"+ {self.job_trackers['backmapping'].nrunning_sims()}")

        nFrames = min(self.nMaxSelectedCGFrameBuffer - nPendingFrames,
                       self.nMaxCGFramesSelectionsPerIter)

        LOGGER.debug(f"nframes = {nFrames} = min({self.nMaxSelectedCGFrameBuffer} "
                     f"- {nPendingFrames}, {self.nMaxCGFramesSelectionsPerIter}")

        if nFrames <= 0:
            return 0

        # select new frames
        LOGGER.info('Select {} new cg frames'.format(nFrames))
        selections = self.cgselector.select(nFrames)

        n = len(selections)
        LOGGER.info('Selected {} cg frames'.format(n))
        if n == 0:
            return n

        selections = [s.id for s in selections]

        # test these selections
        # Split a list of sims based on their status.
        # Returns: sims_success, sims_failed, sims_unknown, sims_stop
        _s, _f, _u, _ = self.job_trackers['backmapping'].split_sims_on_status(selections)

        # this is the correct scenario
        if len(_s) == 0 and len(_f) == 0:
            self.job_trackers['backmapping'].add_to_queue(selections)

        # looks like some of these selections were already simulated
        else:
            LOGGER.error('Found (%d %d) selected patches for which createsim was already run. '
                         'Were you using old sims with new ml?', len(_s), len(_f))
            self.job_trackers['backmapping'].add_to_queue(_u)
            self.job_trackers['aa'].add_to_queue(_s)

        return n

    # --------------------------------------------------------------------------
    def _task_update_jobs(self):

        succeeded = {}
        failed = {}
        started = {}

        # ----------------------------------------------------------------------
        # first, we will update job tracker
        for j in JOB_TYPES:
            succeeded[j], failed[j] = self.job_trackers[j].update()
            LOGGER.debug(f"   JOB_TYPES={j} and self.job_trackers[{j}] = {self.job_trackers[j]}")

        LOGGER.debug('succeeded = {} and failed = {}'.format(succeeded, failed))

        # ----------------------------------------------------------------------
        # add the finished jobs to the next queue!
        for j in JOB_NEXT_QUEUE.keys():
            fjobs = succeeded[j]
            if len(fjobs) > 0:
                LOGGER.debug('Found {} successful {} sims. next, queue them to {}!'.format(len(fjobs), j, JOB_NEXT_QUEUE[j]))
                self.job_trackers[JOB_NEXT_QUEUE[j]].add_to_queue(fjobs)

        LOGGER.debug('JOB_NEXT_QUEUE = {}'.format(JOB_NEXT_QUEUE.keys()))
        # ----------------------------------------------------------------------
        # finally, we will start any new jobs!
        LOGGER.debug(f"JOB_TYPES = {JOB_TYPES} started = {started}")
        for j in JOB_TYPES:
            LOGGER.debug(f"   JOB_TYPES={j} and self.job_trackers[{j}] = {self.job_trackers[j]}")
            nJobs, started[j] = self.job_trackers[j].start_jobs(self.nMaxJobsPerIter)
            LOGGER.debug(f"   nJobs={nJobs}, started[{j}]={started[j]}")

        LOGGER.debug(f"started = {started}, succeeded = {succeeded}, failed = {failed}")
        # ----------------------------------------------------------------------
        # return the dictionaries?
        return started, succeeded, failed


    # --------------------------------------------------------------------------
    # the main workflow task
    # --------------------------------------------------------------------------
    def run_workflow(self, n, lock_patch_io, lock_patch_select):

        p = multiprocessing.current_process()
        signal.signal(signal.SIGTERM, self.signal_wrapper("workflow", p.pid))
        signal.signal(signal.SIGINT, self.signal_wrapper("workflow", p.pid))
        # ----------------------------------------------------------------------
        try:
            self._init_scheduling()
            self._init_cg_selection()
            if self.do_gc:
               self._init_patch_creation()
            if self.do_mlserver:
                self._init_ml_server()
            if self.do_ucgserver:
                self._init_ucg_server()
            self.restore()

            previously_ckpt_jobs = 0
            for j in JOB_TYPES:
                previously_ckpt_jobs += self.job_trackers[j].nrunning_jobs() + self.job_trackers[j].nqueued_sims()

            if previously_ckpt_jobs == 0:
                LOGGER.info(f"No previously checkpointed jobs found (={previously_ckpt_jobs}) in checkpoints, starting workflow")
                self._wf_ready.set()

            slp_time = 0
            loop_timer = Timer()
            while not self._exit.wait(slp_time):

                LOGGER.info('Starting {} iteration {}'.format(p.name, self.iterCounterWF))
                LOGGER.profile('Starting {} iteration {}'.format(p.name, self.iterCounterWF))

                loop_timer.start()

                # --------------------------------------------------------------
                # in the restore phase, let's focus only on starting the jobs
                if not self._wf_ready.is_set():
                    njobs_remaining = 0
                    for j in JOB_TYPES:
                        self.job_trackers[j].start_jobs(self.nMaxJobsPerIter)
                        njobs_remaining += self.job_trackers[j].njobs_2start()

                    if njobs_remaining <= 0:
                        LOGGER.info("Concluding restore phase")
                        self._wf_ready.set()   # kick off other processes

                # --------------------------------------------------------------
                # otherwise, we need to do the regular tasks
                else:
                    # ----------------------------------------------------------
                    # task 1: read new patches/cg frames and add to ML selectors
                    
                    # Option 1. using ML (Latent Space) server to generate and select patches
                    self._task_add_new_patches_mlserver()
                    # Option 2. Using UCG server to generate and select patches
                    self._task_add_new_patches_ucgserver()
                    # Option 3. Using PatchCreator with macro model to create patches (these patches will be selected below)
                    self._task_add_new_patches_to_ml(lock_patch_io, lock_patch_select)

                    LOGGER.profile('    > Iteration {}: added new patches'.format(self.iterCounterWF))
                    self._task_add_cgframes_to_ml()
                    if self.do_cgselection:
                        LOGGER.profile('    > Iteration {}: added CG frame patches'.format(self.iterCounterWF))

                    # ----------------------------------------------------------
                    # 2021.02.07: HB moved this task up
                    # task 2: start the jobs
                    self._task_update_jobs()
                    LOGGER.profile('    > Iteration {}: updated jobs'.format(self.iterCounterWF))

                    # ----------------------------------------------------------
                    # task 3: select new candidates for simulations
                    self._task_select_pfpatches(lock_patch_select)
                    if self.do_patchselection:
                        LOGGER.profile('    > Iteration {}: selected patches'.format(self.iterCounterWF))

                    self._task_select_cgframes()
                    if self.do_cgselection:
                        LOGGER.profile('    > Iteration {}: select CG frames'.format(self.iterCounterWF))

                    # ----------------------------------------------------------
                    # validate state and checkpoint
                    if self.do_patchselection:
                        self.pselector.test()
                        LOGGER.profile('    > Iteration {}: patch selector test'.format(self.iterCounterWF))
                    if self.do_cgselection:
                        self.cgselector.test()
                        LOGGER.profile('    > Iteration {}: CG selector test'.format(self.iterCounterWF))

                # --------------------------------------------------------------
                for j in JOB_TYPES:
                    self.job_trackers[j].test()
                LOGGER.profile('    > Iteration {}: testing all jobs'.format(self.iterCounterWF))
                self.checkpoint()

                # ------------------------------------------------------------------
                LOGGER.info("{} iteration {} finished: {} {}"
                            .format(p.name, self.iterCounterWF, self.is_exit(), self.is_error()))

                # --------------------------------------------------------------
                loop_time = loop_timer.elapsed()
                self.iterCounterWF += 1

                # slp_time = 0 # max(0, self.nSecPerIterWF - loop_time)
                slp_time = max(0, self.nSecPerIterWF - loop_time)
                LOGGER.debug('{} waiting for {} seconds'.format(p.name, slp_time))

            LOGGER.debug("AFTER LOOP: Workflow loop ending.")

        # ----------------------------------------------------------------------
        except Exception as e:
            # self.exception_queue.put(sys.exc_info())
            traceback.print_exc()
            self.error()
            self.exit()
            raise e

        # ----------------------------------------------------------------------
        if self.is_error():
            LOGGER.info('{} process is exiting due to error flag'.format(p.name))
        elif self.is_exit():
            LOGGER.info('{} process is exiting due to exit flag'.format(p.name))

    # --------------------------------------------------------------------------
    # patch creation task
    # --------------------------------------------------------------------------
    def run_patch_creation(self, n, lock_patch_io):

        p = multiprocessing.current_process()
        signal.signal(signal.SIGTERM, self.signal_wrapper('patch_creator', p.pid))

        # ----------------------------------------------------------------------
        try:
            self._wf_ready.wait()
            self._init_patch_creation()

            slp_time = 0
            loop_timer = Timer()
            while not self._exit.wait(slp_time):
                LOGGER.info('Starting {} iteration {}'.format(p.name, self.iterCounterPC))
                LOGGER.profile('Starting {} iteration {}'.format(p.name, self.iterCounterPC))

                loop_timer.start()

                # --------------------------------------------------------------

                patches = self.patch_creator.run_for_macro()

                # if there were some patches created!
                if len(patches) > 0:
                    # write them (use blocking acquire of lock)
                    with lock_patch_io:
                        LOGGER.debug("Acquired lock on patches ({}, {})".format(p.name, p.pid))
                        self.iointerface.save_patches(Naming.dir_root('patches'), patches)
                        self.patch_creator.checkpoint()
                    LOGGER.debug("Released lock on patches ({}, {})".format(p.name, p.pid))

                # --------------------------------------------------------------
                LOGGER.info("{} iteration {} finished: {} {}"
                            .format(p.name, self.iterCounterPC, self.is_exit(), self.is_error()))

                # --------------------------------------------------------------
                loop_time = loop_timer.elapsed()
                self.iterCounterPC += 1


                slp_time = max(0, self.nSecPerIterPC - loop_time)
                LOGGER.info('{} waiting for {} seconds'.format(p.name, slp_time))

            LOGGER.debug("AFTER LOOP: Patch creator loop ending.")

        # ----------------------------------------------------------------------
        except Exception as e:
            # self.exception_queue.put(sys.exc_info())
            traceback.print_exc()
            self.error()
            self.exit()
            raise e

        # ----------------------------------------------------------------------
        if self.is_error():
            LOGGER.info('{} process is exiting due to error flag'.format(p.name))
        elif self.is_exit():
            LOGGER.info('{} process is exiting due to exit flag'.format(p.name))

    # --------------------------------------------------------------------------
    # feedback tasks
    # --------------------------------------------------------------------------
    def run_feedback_cg2macro(self, n, lock_patch_select):

        if not self.do_feedback_cg2mc:
            return

        p = multiprocessing.current_process()
        signal.signal(signal.SIGTERM, self.signal_wrapper('feedback_cg2macro', p.pid))

        # ----------------------------------------------------------------------
        try:
            self._wf_ready.wait()
            self._init_feedback_cg2macro()

            slp_time = 0
            loopTimer = Timer()
            while not self._exit.wait(slp_time):
                LOGGER.info('Starting {} iteration {}'.format(p.name, self.iterCounterFB_cg2macro))
                LOGGER.profile('Starting {} iteration {}'.format(p.name, self.iterCounterFB_cg2macro))
                loopTimer.start()

                # --------------------------------------------------------------
                self.fb_cg2macro.aggregate(lock_patch_selector=lock_patch_select)
                self.fb_cg2macro.report()
                self.fb_cg2macro.checkpoint()

                # --------------------------------------------------------------
                LOGGER.info("{} iteration {} finished: {} {}"
                            .format(p.name, self.iterCounterFB_cg2macro,
                                    self.is_exit(), self.is_error()))

                # --------------------------------------------------------------
                loop_time = loopTimer.elapsed()
                self.iterCounterFB_cg2macro += 1

                slp_time = max(0, self.nSecPerMacroFeedback - loop_time)
                LOGGER.info('{} waiting for {} seconds'.format(p.name, slp_time))

            LOGGER.debug("AFTER LOOP: Feedback[CG2Macro] loop ending.")

        # ----------------------------------------------------------------------
        except Exception as e:
            # self.exception_queue.put(sys.exc_info())
            traceback.print_exc()
            self.error()
            self.exit()
            raise e

        # ----------------------------------------------------------------------
        if self.is_error():
            LOGGER.info('{} process is exiting due to error flag'.format(p.name))
        elif self.is_exit():
            LOGGER.info('{} process is exiting due to exit flag'.format(p.name))

    def run_feedback_aa2cg(self, n, _fake_argument):

        if not self.do_feedback_aa2cg:
            return

        p = multiprocessing.current_process()
        signal.signal(signal.SIGTERM, self.signal_wrapper('feedback_aa2cg', p.pid))

        # ----------------------------------------------------------------------
        try:
            self._wf_ready.wait()
            self._init_feedback_aa2cg()

            slp_time = 0
            loopTimer = Timer()
            while not self._exit.wait(slp_time):
                LOGGER.info('Starting {} iteration {}'.format(p.name, self.iterCounterFB_aa2cg))
                loopTimer.start()

                self.fb_aa2cg.aggregate()
                self.fb_aa2cg.report()
                self.fb_aa2cg.checkpoint()

                # --------------------------------------------------------------
                LOGGER.info("{} iteration {} finished: {} {}"
                            .format(p.name, self.iterCounterFB_aa2cg,
                                    self.is_exit(), self.is_error()))

                # --------------------------------------------------------------
                loop_time = loopTimer.elapsed()
                self.iterCounterFB_aa2cg += 1

                slp_time = max(0, self.nSecPerCGFeedback - loop_time)
                LOGGER.info('{} waiting for {} seconds'.format(p.name, slp_time))

            LOGGER.debug("AFTER LOOP: Feedback[AA2CG] loop ending.")

        # ----------------------------------------------------------------------
        except Exception as e:
            # self.exception_queue.put(sys.exc_info())
            traceback.print_exc()
            self.error()
            self.exit()
            raise e

        # ----------------------------------------------------------------------
        if self.is_error():
            LOGGER.info('{} process is exiting due to error flag'.format(p.name))
        elif self.is_exit():
            LOGGER.info('{} process is exiting due to exit flag'.format(p.name))

    # --------------------------------------------------------------------------
    # Stop main processes
    # --------------------------------------------------------------------------
    def stop(self):
        LOGGER.info('marking the wfmanager to stop!')
        self.exit()
        sys.exit(1)

    def join(self):
        for key, process in self.processes.items():
            LOGGER.info("JOINING: Attempting to join process '%s'...", key)
            process.join()

    # --------------------------------------------------------------------------
    # the main run task
    # --------------------------------------------------------------------------
    def start(self):

        self._setup.wait()

        s = "gc = {}, patchcreation = {}, workflow = {}, mlserver = {}, ucgserver = {}, feedback = {}/{}, schedulejobs = {}"
        s = s.format(self.do_gc, self.do_patchcreation, self.do_workflow, self.do_mlserver, self.do_ucgserver,
                     self.do_feedback_cg2mc, self.do_feedback_aa2cg, self.do_schedulejobs)

        LOGGER.info("Starting workflow_manager ({})".format(s))
        lock_patch_io = multiprocessing.Lock()
        lock_patch_select = multiprocessing.Lock()

        # ----------------------------------------------------------------------
        # May 14, 2024: PatchSelector is deprecated for now
        # ProxyManager.register("PatchSelector", PatchSelector)
        # self.proxy_manager = ProxyManager()
        # self.proxy_manager.start()
        # self._init_patch_selection()

        # ----------------------------------------------------------------------
        if self.do_patchcreation and not self.do_gc:
            _process = multiprocessing.Process(
                name='patch_creator',
                args=(1, lock_patch_io),
                target=self.run_patch_creation)
            _process.daemon = False
            self.processes['patch_creator'] = _process

        if self.do_workflow:
            _process = multiprocessing.Process(
                name='wf_manager',
                args=(0, lock_patch_io, lock_patch_select),
                target=self.run_workflow)
            _process.daemon = False
            self.processes['wf_manager'] = _process
        else:
            self._wf_ready.set()

        if self.do_feedback_cg2mc:
            _process = multiprocessing.Process(
                name='fb_cg2macro_manager',
                args=(2, lock_patch_select),
                target=self.run_feedback_cg2macro)
            _process.daemon = False
            self.processes['fb_cg2macro_manager'] = _process

        if self.do_feedback_aa2cg:
            _process = multiprocessing.Process(
                name='fb_aa2cg_manager',
                args=(3, None),
                target=self.run_feedback_aa2cg)
            _process.daemon = False
            self.processes['fb_aa2cg_manager'] = _process

        # ----------------------------------------------------------------------
        for key, process in self.processes.items():
            LOGGER.info("STARTING: Starting process '%s'...", key)
            process.start()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
