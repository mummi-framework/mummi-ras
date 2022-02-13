# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
import os
import signal
import sys
import yaml
import timeit
import traceback
import multiprocessing
from multiprocessing.managers import SyncManager

from mummi_core.workflow.job import JOB_TYPES, JOB_NEXT_QUEUE
from mummi_core.workflow.flux_env import flux_uri
from mummi_core.workflow.jobTracker import JobTracker

import mummi_core
import mummi_ras
from mummi_core.utils.timer import Timer
from mummi_ras import Naming
from mummi_ras.transformations.patch_creator import MacroPatchCreator
from mummi_ras.ml import PatchSelector, CGSelector, CGSelectorType
from mummi_ras.feedback.feedback_manager_macro import MacroFeedbackManager, FeedbackManagerType
from mummi_ras.feedback.feedback_manager_aatocg import FeedbackManager_AA2CG

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

            traceback.print_exc()
            #self.exception_queue.put(sys.exc_info())
            #self.error()  # copied here but not sure if we want this
            #self.exit()

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

            # ------------------------------------------------------------------
            # identify scheduler interface
            self.flux = flux_uri()
            LOGGER.info(f'Identified flux_uri = [{self.flux}]')
            if self.flux == '':
                LOGGER.error('Failed to identify flux uri. Disabling job scheduling')
                self.do_schedulejobs = False
                assert False, 'Failed to identify flux uri'

            # ------------------------------------------------------------------
            # if no workspace has been specified, use the default in hierarchy
            self.workspace = self.wconfig.get('workspace', '')
            if self.workspace == '':
                self.workspace = Naming.dir_root('workspace')

            os.makedirs(self.workspace, exist_ok=True)
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
        self.adapter_type = self.adapter_batch['type']

        # ----------------------------------------------------------------------
        # get the nodes and split them for different jobs
        nnodes_tot = int(os.environ['MUMMI_NNODES'])

        # Fixed cost Infrastructure nodes
        nnodes_flux: int = 0   # Dedicated full node for Flux
        nnodes_wf: int = 1     # Dedicated CPU-only for wfmanager

        # Variable cost Infrastructure nodes
        # CPU-only Services
        #nnodes_dbr: int = int(os.environ.get('MUMMI_DBR_TOTAL_NNODES', 0))
        nnodes_mac: int = int(os.environ.get('MUMMI_MACRO_NNODES', 0))
        nnodes_redis: int = int(os.environ.get('MUMMI_REDIS_NNODES', 0))
        nnodes_infra: int = nnodes_mac + nnodes_redis

        # Node computations
        nnodes_tot:int = nnodes_tot - nnodes_flux
        nnodes_gpu: int = nnodes_tot
        nnodes_cpu: int = nnodes_tot - nnodes_wf - nnodes_infra

        portion_csim = float(self.wconfig['portion_csim'])
        portion_cg = float(self.wconfig['portion_cg'])
        assert 0 <= portion_csim <= 1
        assert 0 <= portion_cg <= 1

        nnodes_csim : int = int(nnodes_cpu * portion_csim)
        nnodes_cg : int = int(nnodes_gpu * portion_cg)
        nnodes_bmap : int = nnodes_cpu - nnodes_csim
        nnodes_aa : int = nnodes_gpu - nnodes_cg

        LOGGER.info(f'NNodes: total = {nnodes_tot}: '
                    f'(csim = {nnodes_csim}, cg = {nnodes_cg},'
                    f' bmap = {nnodes_bmap}, aa = {nnodes_aa})')


        # ----------------------------------------------------------------------
        
        self.job_trackers = {
            'createsim': JobTracker(self.config['createsim'],
                            nnodes_csim, self.iointerface,
                            self.adapter_batch,
                            self.do_schedulejobs),
                            # self.exception_queue, self.do_schedulejobs),

            'cg': JobTracker(self.config['cg'],
                            nnodes_csim, self.iointerface,
                            self.adapter_batch,
                            self.do_schedulejobs),
                            # self.exception_queue, self.do_schedulejobs),

            'backmapping': JobTracker(self.config['backmapping'],
                            nnodes_bmap, self.iointerface,
                            self.adapter_batch,
                            self.do_schedulejobs),
                            # self.exception_queue, self.do_schedulejobs),

            'aa': JobTracker(self.config['aa'],
                            nnodes_aa, self.iointerface,
                            self.adapter_batch,
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

    def _init_patch_selection(self):

        self.pselector = None
        if self.do_patchselection or self.do_feedback_cg2mc:
            ml_config = Naming.ml('macro')
            stype = ml_config['selection']
            assert stype in ['importance', 'random']

            wspace = os.path.join(Naming.dir_root('ml'), 'macro')
            self.pselector = self.proxy_manager.PatchSelector(stype, wspace, ml_config)
            self.pselector.restore()

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

    # --------------------------------------------------------------------------
    def checkpoint(self):

        state = dict(
            flux=self.flux,
            iterCounterWF=self.iterCounterWF,
            iterCounterPC=self.iterCounterPC,
            patchCounter=self.patchCounter,
            jobs_createsim=self.job_trackers['createsim'].status(),
            jobs_backmapping=self.job_trackers['backmapping'].status(),
            jobs_cg=self.job_trackers['cg'].status(),
            jobs_aa=self.job_trackers['aa'].status())

        self.iointerface.save_checkpoint(self.chkpt, state, use_tstamp=True)

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
            return 0

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
            return 0

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

        LOGGER.debug(f"pending patches = {nPendingPatches} = {self.job_trackers['cg'].nqueued_sims()} + {self.job_trackers['createsim'].nqueued_sims()} + {self.job_trackers['createsim'].nrunning_sims()}")
        nPatches = min(self.nMaxSelectedPatchBuffer - nPendingPatches,
                       self.nMaxPatchesSelectionsPerIter)

        LOGGER.debug(f"npatches = {nPatches} = min({self.nMaxSelectedPatchBuffer} - {nPendingPatches}, {self.nMaxPatchesSelectionsPerIter}")

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
        _s, _f, _u = self.job_trackers['createsim'].split_sims_on_status(selections)

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

    # --------------------------------------------------------------------------
    def _task_select_cgframes(self):

        if not self.do_cgselection:
            return 0

        nPendingFrames = self.job_trackers['aa'].nqueued_sims() + \
                          self.job_trackers['backmapping'].nqueued_sims() + \
                          self.job_trackers['backmapping'].nrunning_sims()

        LOGGER.debug(f"pending frames = {nPendingFrames} = {self.job_trackers['aa'].nqueued_sims()} + {self.job_trackers['backmapping'].nqueued_sims()} + {self.job_trackers['backmapping'].nrunning_sims()}")

        nFrames = min(self.nMaxSelectedCGFrameBuffer - nPendingFrames,
                       self.nMaxCGFramesSelectionsPerIter)

        LOGGER.debug(f"nframes = {nFrames} = min({self.nMaxSelectedCGFrameBuffer} - {nPendingFrames}, {self.nMaxCGFramesSelectionsPerIter}")

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
        _s, _f, _u = self.job_trackers['backmapping'].split_sims_on_status(selections)

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

        # ----------------------------------------------------------------------
        # add the finished jobs to the next queue!
        for j in JOB_NEXT_QUEUE.keys():
            fjobs = succeeded[j]
            if len(fjobs) > 0:
                LOGGER.info('Found {} successful {} sims. next, queue them to {}!'.format(len(fjobs), j, JOB_NEXT_QUEUE[j]))
                self.job_trackers[JOB_NEXT_QUEUE[j]].add_to_queue(fjobs)

        # ----------------------------------------------------------------------
        # finally, we will start any new jobs!
        for j in JOB_TYPES:
            nJobs, started[j] = self.job_trackers[j].start_jobs(self.nMaxJobsPerIter)

        # ----------------------------------------------------------------------
        # return the dictionaries?


    # --------------------------------------------------------------------------
    # the main workflow task
    # --------------------------------------------------------------------------
    def run_workflow(self, n, lock_patch_io, lock_patch_select):

        p = multiprocessing.current_process()
        signal.signal(signal.SIGTERM, self.signal_wrapper('workflow', p.pid))

        # ----------------------------------------------------------------------
        try:
            self._init_scheduling()
            self._init_cg_selection()
            if self.do_gc:
               self._init_patch_creation()

            self.restore()

            slp_time = 0
            loop_timer = Timer()
            while not self._exit.wait(slp_time):

                LOGGER.info('Starting {} iteration {}'.format(p.name, self.iterCounterWF))
                LOGGER.profile('Starting {} iteration {}'.format(p.name, self.iterCounterWF))

                loop_timer.start()

                # --------------------------------------------------------------
                # in the restore phase, let's focus only on starting the jobs
                if not self._wf_ready.is_set():

                    LOGGER.info('Doing restore phase')
                    njobs_remaining = 0
                    for j in JOB_TYPES:
                        self.job_trackers[j].start_jobs(self.nMaxJobsPerIter)
                        njobs_remaining += self.job_trackers[j].njobs_2start()

                    if njobs_remaining <= 0:
                        LOGGER.info('Concluding restore phase')
                        self._wf_ready.set()   # kick off other processes

                # --------------------------------------------------------------
                # otherwise, we need to do the regular tasks
                else:
                    # ----------------------------------------------------------
                    # task 1: read new patches/cg frames and add to ML selectors
                    self._task_add_new_patches_to_ml(lock_patch_io, lock_patch_select)
                    self._task_add_cgframes_to_ml()

                    # ----------------------------------------------------------
                    # 2021.02.07: HB moved this task up
                    # task 2: start the jobs
                    self._task_update_jobs()

                    # ----------------------------------------------------------
                    # task 3: select new candidates for simulations
                    self._task_select_pfpatches(lock_patch_select)
                    self._task_select_cgframes()

                    # ----------------------------------------------------------
                    # validate state and checkpoint
                    if self.do_patchselection:
                        self.pselector.test()
                    if self.do_cgselection:
                        self.cgselector.test()

                # --------------------------------------------------------------
                for j in JOB_TYPES:
                    self.job_trackers[j].test()
                self.checkpoint()

                # ------------------------------------------------------------------
                LOGGER.info("{} iteration {} finished: {} {}"
                            .format(p.name, self.iterCounterWF, self.is_exit(), self.is_error()))

                # --------------------------------------------------------------
                loop_time = loop_timer.elapsed()
                self.iterCounterWF += 1

                slp_time = 0 # max(0, self.nSecPerIterWF - loop_time)
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

    def join(self):
        for key, process in self.processes.items():
            LOGGER.info("JOINING: Attempting to join process '%s'...", key)
            process.join()

    # --------------------------------------------------------------------------
    # the main run task
    # --------------------------------------------------------------------------
    def start(self):

        self._setup.wait()

        s = 'gc = {}, patchcreation = {}, workflow = {}, feedback = {}/{}, schedulejobs = {}'
        s = s.format(self.do_gc, self.do_patchcreation, self.do_workflow,
                     self.do_feedback_cg2mc, self.do_feedback_aa2cg,
                     self.do_schedulejobs)

        LOGGER.info('Starting workflow_manager ({})'.format(s))
        lock_patch_io = multiprocessing.Lock()
        lock_patch_select = multiprocessing.Lock()

        # ----------------------------------------------------------------------
        ProxyManager.register('PatchSelector', PatchSelector)
        self.proxy_manager = ProxyManager()
        self.proxy_manager.start()
        self._init_patch_selection()

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
