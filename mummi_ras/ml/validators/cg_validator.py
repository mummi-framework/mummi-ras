#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------

import os
import logging
import numpy as np
import subprocess as subp
import multiprocessing
from typing import List, Any
from timeit import default_timer as timer
from shutil import rmtree, which
import signal

# import MDAnalysis as mda
# from MDAnalysis.analysis import distances
# from MDAnalysis.analysis.rms import rmsd

from .validator import Validator

LOGGER = logging.getLogger(__name__)

class CGBaseValidator(Validator):
    """
    Class to validate the generated CG protein structures.
    """
    def __init__(self, iteration_path: str, resource_path: str, complex_name: str, nthreads: int = None, healing: bool = True, cleanup: bool = True, mini_mummi: bool = False):
        super().__init__()
        self.cleanup = cleanup
        self.healing = healing
        self.iteration_path = iteration_path
        self.resource_dir: str = resource_path
        self.complex_name: str = complex_name # Name of the complex .gro
        self.max_mf = float("1.0e+06")
        self.max_pe = float("1.0e+07")
        self.mini_mummi = mini_mummi

        try:
            self.num_cores_prev = nthreads or len(os.sched_getaffinity(0))
        except AttributeError as e:
            self.num_cores_prev = nthreads or multiprocessing.cpu_count()

        # self.num_cores = max(1, (self.num_cores_prev // self.get_logical_count()) - 2) # we want physical cores count
        # No need to divide by get_logical_count because when run inside of Flux, Flux set the right thing and the number of cores = number of threads
        self.num_cores = max(1, int(self.num_cores_prev) - 2) # We substract 2 for the ML server sampler/generator/feedback

        self.mpcontext = "fork" # fork or spawn

        self.current_rpath = iteration_path
        LOGGER.info(f"> Initialization of validator with")
        LOGGER.info(f"  > iteration_path        = {self.iteration_path}")
        LOGGER.info(f"  > resource_path         = {self.resource_dir}")
        LOGGER.info(f"  > complex_name          = {self.complex_name}")
        LOGGER.info(f"  > maximum force         = [-{self.max_mf}, {self.max_mf}]")
        LOGGER.info(f"  > potential energy      = [-{self.max_pe}, {self.max_pe}]")
        LOGGER.info(f"  > cleanup               = {self.cleanup}")
        LOGGER.info(f"  > healing               = {self.healing}")
        LOGGER.info(f"  > logical cores         = {self.num_cores_prev}")
        LOGGER.info(f"  > available cores       = {self.num_cores}")

        gmx_path = which("gmx")
        if gmx_path is None:
            LOGGER.error(f"GROMACS 'gmx' not found. Please verify that GROMACS is installed and reachable. Aborting")
            raise RuntimeError("GROMACS executable 'gmx' is missing in PATH.")

        if not os.path.isdir(self.resource_dir):
            LOGGER.error(f"The validator directory {self.resource_dir} is not valid.")
            raise NotADirectoryError(f"Invalid validator directory: {self.resource_dir}")

        complex_file = os.path.isfile(os.path.join(self.resource_dir, self.complex_name))
        if not complex_file:
            LOGGER.error(f"The complex file {complex_file} does not exist.")
            raise FileNotFoundError(f"{complex_file} does not exist.")

    def get_logical_count(self):
        """
        Return the number of logical cores per physical cores.
        On POWERPC, psutil is not working properly
        """
        p = subp.run(["lscpu"], capture_output=True, text=True)
        for l in p.stdout.split('\n'):
            if l.startswith("Thread(s) per core:"):
                try:
                    nthreads_per_core = int(l.split()[-1])
                except ValueError as e:
                    LOGGER.error(f"Could not find the number of threads/core in lscpu output")
                    nthreads_per_core = None
        return nthreads_per_core

    def num_validated(self, iteration_id: int) -> int:
        """
        Return the number of validated structures for a given structure.

        :param iteration_id: Iteration ID (positive integer starting at 0)
        :type iteration_id: int
        
        :return: Number of validated structures
        :rtype: int
        """

        iteration_id = "{:02d}".format(iteration_id)
        RPATH = os.path.join(self.iteration_path, f"iter{iteration_id}")
        file = os.path.join(RPATH, f"table_all_predictions.npz")
        nb_valid = 0
        if os.path.isfile(file):
            data = np.load(file, allow_pickle=True)["status"]
            nb_valid = sum(data == True)
        return nb_valid

    def validate(self, iteration_id: int, name: str, positions: np.array) -> List:
        iteration_id = "{:02d}".format(iteration_id)
        RPATH = os.path.join(self.iteration_path, f"iter{iteration_id}")
        return self._internal_validate([name, positions, RPATH, self.resource_dir, self.complex_name, self.max_mf, self.max_pe, self.healing, self.cleanup, self.mini_mummi])

    def validateArray(self, iteration_id: int, nameArray: np.array, positionsArray: np.array) -> np.array:
        results = []
        start = timer()
        for i in range(len(positionsArray)):
            results.append(self.validate(iteration_id, nameArray[i], positionsArray[i]))
        LOGGER.info(f"Validator returned with {len(results)} elements and took: {timer()-start} sec ({len(results) / (timer()-start)} sample/sec)")
        return np.array(results, dtype=object)

    def init_worker(self):
        """ Needed to avoid that worker catch signal intended for the server. """
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGTERM, signal.SIG_IGN)

    def _cleanup_dir(self, cwd: str, path: str, change_dir: bool = True) -> None:
        if not self.cleanup:
            return
        try:
            LOGGER.debug(f"Cleaning up {path}")
            rmtree(path)
            if change_dir:
                os.chdir(cwd)
        except Exception as e:
            LOGGER.debug(f"Cleanup failed {path} ({e.__class__.__name__}: {e})")

    def validateArrayThreaded(self, iteration_id: int, nameArray: np.array, positionsArray: np.array, num_proc: int = None) -> np.array:
        assert self.mpcontext in ["fork", "spawn"]
        dataArray = []
        iteration_id = "{:02d}".format(iteration_id)
        RPATH = os.path.join(self.iteration_path, f"iter{iteration_id}")
        if not num_proc:
            num_proc = self.num_cores
        num_proc = min(num_proc, len(positionsArray))

        for i in range(len(positionsArray)):
            dataArray.append([nameArray[i], positionsArray[i], RPATH, self.resource_dir, self.complex_name, self.max_mf, self.max_pe, self.healing, self.cleanup, self.mini_mummi])
            # We make dir for each structure not in parallel
            # Sometimes Lustre fails when to many conccurent mkdir
            fullPathDir = os.path.join(RPATH, nameArray[i])
            os.makedirs(fullPathDir, exist_ok=True)

        start = timer()

        with multiprocessing.get_context(self.mpcontext).Pool(num_proc, self.init_worker) as pool:
            try:
                results = pool.map(self._internal_validate, dataArray, chunksize = 1)
            except Exception as e:
                LOGGER.error(f"validateArrayThreaded: {e}")
                pool.terminate()
                results = []

        # with ProcessPoolExecutor(max_workers = num_proc) as pool:
        #     results = list(pool.map(CGValidator._internal_validate, dataArray, chunksize = 1))

        LOGGER.info(f"Validator returned with {len(results)} elements and took: {timer()-start} sec ({len(results) / (timer()-start)} sample/sec)")
        return np.array(results, dtype=object)
