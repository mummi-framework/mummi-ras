#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2023, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------
# @authors
#           Loic Pottier <pottier1@llnl.gov>
# ------------------------------------------------------------------------------

import os
import signal
import traceback
import numpy as np
from typing import List, Union, Any
import time
import datetime
import tempfile
import sys
from timeit import default_timer as timer
import multiprocessing
from multiprocessing.pool import ThreadPool as Pool
import threading

import mummi_core

from mummi_ras.feedback.sampling import LSSamplingDB

from logging import getLogger, Logger
LOGGER = getLogger(__name__)

class FeedbackCreatesims_C3:
    """
    Class that constantly check the status (failure, success or running) of createsims
    """
    def __init__(self,
                 dir_to_watch: str,
                 sampling_db: str,
                 model_name: str,
                 ls_encoder,
                 waittime: int = 60,
                 nthreads: int = 1,
                 lockfile: str = None,
                 logger: Logger = LOGGER
        ) -> None:
        self.sampling_db = sampling_db
        self.model_name = model_name
        self.nthreads = nthreads
        self.logger = logger
        self.hostname = mummi_core.get_hostname(contract_hostname = False)

        self.ls_encoder = ls_encoder

        # Files to look for
        self.file_patterns = ["createsims_success", "createsims_failure"]
        self.lockfile = lockfile or tempfile.mkstemp(prefix="mummi-sampling-", suffix=".lock", dir=None)[1]
        self.waittime = waittime # in seconds

        self.dir_to_watch = dir_to_watch

        p = multiprocessing.current_process()
        signal.signal(signal.SIGTERM, self.signal_wrapper("Feedback_Createsims", p.pid))
        signal.signal(signal.SIGINT, self.signal_wrapper("Feedback_Createsims", p.pid))

        self.logger.info(f"Starting FeedbackCreatesims")
        self.logger.info(f" > Directory               {self.dir_to_watch}")
        self.logger.info(f" > ML model                {self.model_name}")
        self.logger.info(f" > Feeback Sampling DB     {self.sampling_db}")
        self.logger.info(f" > Every (seconds)         {self.waittime}")
        self.logger.info(f" > Threads                 {self.nthreads}")

    def signal_wrapper(self, name, pid):
        def handler(signum, frame):
            print("")
            self.logger.warning(f"Received (SIGNUM={signum}) for {name}[pid={pid}], stopping the process...")
            try:
                self.logger.debug(f'Removing {self.lockfile}')
                os.remove(self.lockfile)
            except FileNotFoundError as _:
                self.logger.warning(f'Lock file not found {self.lockfile}')
                pass
            sys.exit(1)
        return handler

    def launch(self):
        # We declare the LSSampling DB
        self.database = LSSamplingDB(
            database_path = self.sampling_db,
            model_name = self.model_name,
            ls_encoder = self.ls_encoder,
            lockfile = self.lockfile,
            logger = self.logger
        )

        self.logger.info(f"Start monitoring for feedback: {self.dir_to_watch}")

        # We call it only once to sstart
        while True:
            self.watch_createsims(nthreads = self.nthreads)
            time.sleep(self.waittime)
            # self.logger.info(self.database.print_db())

    def process_createsims(self, folder: os.DirEntry, file_patterns: List[str]) -> dict:
        """
        Check all directories in createsims_dir to see if they contains createsims_success or createsims_failure.
        If status =
            - 0 => success
            - 1 => Failure
            - 2 => unknown (default status, neither success or failure)
        """

        struct_name = os.path.basename(folder.path)
        success = os.path.exists(os.path.join(folder.path, file_patterns[0]))
        failure = os.path.exists(os.path.join(folder.path, file_patterns[1]))
        status = 2 # default status
        tid = threading.current_thread().name

        if success:
            status = 1 # True
        if failure:
            status = 0 # False
        if success and failure:
            self.logger.warning(f"[{tid}] Found success and failure for {struct_name}")

        self.logger.debug(f"[{tid}] Checking {folder.path} => found {status}")
        return status

    def watch_createsims(self, nthreads: int = 1):
        """
        Check all directories in createsims_dir to see if they contains createsims_success or createsims_failure.
        then update sampling_db file. This function will wait for "wait" seconds between each check.
        """
        start = timer()

        def print_failed(caught):
            traceback.print_exc(file = sys.stdout)

        createsims_status = {}
        futures = {}
        with Pool(nthreads) as pool:
            list_folders = list(os.scandir(self.dir_to_watch))
            for folder in list_folders:
                struct_name = os.path.basename(folder.path)
                futures[struct_name] = pool.apply_async(self.process_createsims, (folder, self.file_patterns), error_callback = print_failed)

            for k,v in futures.items():
                createsims_status[k] = v.get()

        end = timer() - start
        if len(list_folders) > 0:
            self.logger.info(f"Scanned {len(list_folders)} structures in {end:.3f} seconds ({len(list_folders)/end} struct/sec) with {nthreads} threads")

        data = self.database.read_data()
        structure_names = data["structure_names"]
        previous_createsims_status = data["createsims_status"]
        # v_coords = data["v_coords"]
        # v2ucg_coords = data["v2ucg_coords"]
        # ucg2cg_coords = data["ucg2cg_coords"]

        nb_updated_structures = 0
        for struct, status in zip(structure_names, previous_createsims_status):
            # if struct not in createsims_status:
            #     nb_updated_structures += 1
            if struct in createsims_status and status != createsims_status[struct]:
                nb_updated_structures += 1

        new_structure_names = []
        # validation_status = []
        # new_v_coords = []
        # new_v2ucg_coords = []
        # new_ucg2cg_coords = []

        updated_createsims = list(previous_createsims_status)

        # Update existing createsims status
        for updated_struct, status in createsims_status.items():
            self.logger.debug(f"Updating {updated_struct} status={status}")
            index = np.where(structure_names == updated_struct)[0]
            if len(index) == 0:
                self.logger.warning(
                    f"{updated_struct} is running as a createsims but "
                    f"is not present in that DB. Index = {index}")

                # # We add a new structure that do not have any ls coords
                # # this could happen if the previous DB got corrupted
                # # We fill the v's, v2ucg and ucg2cg's with zeroes
                # new_structure_names.append(updated_struct)
                # new_v_coords.append(np.empty(v_coords.shape[1:]))
                # new_v2ucg_coords.append(np.empty(v2ucg_coords.shape[1:]))
                # new_ucg2cg_coords.append(np.empty(ucg2cg_coords.shape[1:]))
                # validation_status.append(False)
                # updated_createsims.append(status)
            elif len(index) == 1:
                k = index[0]
                try:
                    if status != updated_createsims[k]:
                        self.logger.info(f"Updated {updated_struct}: {updated_createsims[k]} => {status}")
                    updated_createsims[k] = status
                except IndexError as e:
                    self.logger.warning(f"Could not update {updated_struct} as {index} is out of bound ({len(updated_createsims)})")
            else:
                for k in index:
                    try:
                        if status != updated_createsims[k]:
                            self.logger.info(f"Updated {updated_struct}: {updated_createsims[k]} => {status}")
                        updated_createsims[k] = status
                    except IndexError as e:
                        self.logger.warning(f"Could not update {updated_struct} as {index} is out of bound ({len(updated_createsims)})")
                    self.logger.error(f"{updated_struct} has more than one entry in {self.sampling_db}. The DB is likely corrupted. (Index={index})")

        updated_field = self.database.update_db(
            structure_names = new_structure_names,
            # v_coords = new_v_coords,
            # v2ucg_coords = new_v2ucg_coords,
            # ucg2cg_coords = new_ucg2cg_coords,
            # validation_status = validation_status,
            createsims_status = updated_createsims
        )
        end = timer() - start

        if nb_updated_structures > 0:
            self.logger.info(f"Updated {nb_updated_structures} structures in {end:.3f} seconds ({nb_updated_structures/end} struct/sec) with {nthreads} threads")

