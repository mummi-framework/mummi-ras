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
import MDAnalysis as mda
from typing import List, Any
from timeit import default_timer as timer
from shutil import rmtree, which

from .validator import Validator

LOGGER = logging.getLogger(__name__)

class UCGValidator(Validator):
    """
    Class to validate the generated (U)CG protein structures. This class adds a healing phase
    to the validation step.
    """
    def __init__(self, iteration_path: str, resource_path: str, cleanup: bool = True):
        super().__init__()
        self.cleanup = cleanup
        self.iteration_path = iteration_path
        self.resource_dir: str = resource_path
        self.max_mf = float("1.0e+06")
        self.max_pe = float("1.0e+06")
        self.num_cores = len(os.sched_getaffinity(0))
        self.mpcontext = "fork" # spawn
        # self.cwd = os.path.join(os.getcwd(), iteration_path)
        self.current_rpath = iteration_path
        LOGGER.info(f"> Initialization of validator with")
        LOGGER.info(f"  > iteration_path        = {iteration_path}")
        LOGGER.info(f"  > resource_path         = {resource_path}")
        LOGGER.info(f"  > maximum force         = [-{self.max_mf}, {self.max_mf}]")
        LOGGER.info(f"  > potential energy      = [-{self.max_pe}, {self.max_pe}]")
        LOGGER.info(f"  > num_cores             = {self.num_cores}")
        LOGGER.info(f"  > mpcontext             = {self.mpcontext}")

        gmx_path = which("gmx")
        if gmx_path is None:
            LOGGER.error(f"GROMACS 'gmx' not found. Please verify that GROMACS is installed and reachable. Aborting")
            raise ValueError("GROMACS executable 'gmx' is missing.")

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
            nb_valid = sum(data == "True")
        return nb_valid

    def validate(self, iteration_id: int, name: str, positions: np.array) -> List:
        iteration_id = "{:02d}".format(iteration_id)
        RPATH = os.path.join(self.iteration_path, f"iter{iteration_id}")
        return UCGValidator._internal_validate([name, positions, RPATH, self.resource_dir, self.max_mf, self.max_pe, self.cleanup])

    def validateArray(self, iteration_id: int, nameArray: np.array, positionsArray: np.array) -> np.array:
        retStatus = []
        for i in range(len(positionsArray)):
            retStatus.append(self.validate(iteration_id, nameArray[i], positionsArray[i]))
        return np.array(retStatus)

    def validateArrayThreaded(self, iteration_id: int, nameArray: np.array, positionsArray: np.array, num_proc: int = None) -> np.array:
        assert self.mpcontext in ["fork", "spawn"]
        dataArray = []
        iteration_id = "{:02d}".format(iteration_id)
        RPATH = os.path.join(self.iteration_path, f"iter{iteration_id}")
        if not num_proc:
            num_proc = self.num_cores
        num_proc = min(num_proc, len(positionsArray))

        for i in range(len(positionsArray)):
            dataArray.append([nameArray[i], positionsArray[i], RPATH, self.resource_dir, self.max_mf, self.max_pe, self.cleanup])
        start = timer()

        with multiprocessing.get_context(self.mpcontext).Pool(num_proc) as pool:
            results = pool.map(UCGValidator._internal_validate, dataArray, chunksize = 1)

        LOGGER.info(f"Validator returned with {len(results)} elements and took {timer()-start}s")
        return results

    def write_validation_info(self, iteration_id: int, return_array: np.array, all_structure_names: np.array) -> tuple[np.array, np.array]:

        iteration_id = "{:02d}".format(iteration_id)
        RPATH = os.path.join(self.iteration_path, f"iter{iteration_id}")
        self.current_rpath = RPATH
        if not os.path.isdir(RPATH):
            os.makedirs(RPATH, exist_ok=True)

        npatches = return_array.shape[0]
        files = np.empty(npatches, dtype='object')
        master_file = os.path.join(RPATH, f"master_iter{iteration_id}.npz")

        # Find valid patches
        valid_ids = np.where(return_array[:,0] == 'True')[0]
        valid_files = return_array[valid_ids, 1]

        # files = np.append(files, valid_files)
        files = np.array(valid_files)

        # we do not erase the previous master file
        if os.path.isfile(master_file):
            previous_data = np.load(master_file, allow_pickle=True)["files"]
            files = np.append(previous_data, files)

        # # Save the list of the files into the master file
        np.savez_compressed(master_file, files = files)
        
        # # Save the details for all the predictions
        name_table_all_predictions = os.path.join(os.path.join(RPATH, "table_all_predictions.npz"))
        status = return_array[:,0]
        gromacs_files = return_array[:,1]
        maximum_force = return_array[:,2]
        potential_energy = return_array[:,3]

        if os.path.isfile(name_table_all_predictions):
            previous_data = np.load(name_table_all_predictions, allow_pickle=True)
            all_structure_names = np.append(previous_data['structure_names'], all_structure_names)
            status = np.append(previous_data['status'], status)
            gromacs_files = np.append(previous_data['gromacs_files'], gromacs_files)
            if "maximum_force" in previous_data.files:
                maximum_force = np.append(previous_data['maximum_force'], maximum_force)
            if "potential_energy" in previous_data.files:
                potential_energy = np.append(previous_data['potential_energy'], potential_energy)

        np.savez_compressed(
            name_table_all_predictions,
            structure_names = all_structure_names,
            files = files,
            status = status,
            gromacs_files = gromacs_files,
            maximum_force = maximum_force,
            potential_energy = potential_energy
        )

        LOGGER.info(f'Found {valid_ids.shape[0]} valid patches out of {npatches} patches')
        return valid_ids, valid_files

    @staticmethod
    def _internal_validate(data: List[Any]) -> List:
        name, positions, cwd, resource_dir, max_mf, max_pe, cleanup = data
        gromacs_file = None
        valid_status = False
        start = timer()

        # Make temp dir to run minimization in
        fullPathDir = os.path.join(cwd, name)
        old_cwd = os.getcwd()
        if os.path.isdir(fullPathDir) != True:
            os.mkdir(fullPathDir)
        os.chdir(fullPathDir)

        # Save .gro file - @CHECK faster if MDAnalysis universe is in class and here only replace and save BUT maybe prolematic for later threaded version ???
        gromacs_ref_file = os.path.join(resource_dir, 'RAS-RAF_state_ma_0.gro')
        u = mda.Universe(gromacs_ref_file)
        u.atoms.positions = positions
        u.atoms.write(f"{name}.gro")
        pe = np.inf
        mf = np.inf

        # Run minimization
        try:
            LOGGER.info(f"Run GROMACS minimization in dir {name}")

            # Run energy minimization (steepest descent, with protein constreins)
            grompp = subp.run(f'gmx grompp -c {name}.gro -f {resource_dir}/martini_v2.x_new-rf-em-sc-test.mdp -p {resource_dir}/system-sc.top -o topol-sc.tpr -maxwarn 5', shell=True, capture_output=True)
            # @NOTE: some gmx versions are compiled without thread-MPI and do not support setting the number of threads "-nt" so fix OMP_NUM....
            mdrun = subp.run(f'export OMP_NUM_THREADS=1; gmx mdrun -nt 1 -v -deffnm topol-sc -c {name}-em-sc.gro  > md-sc.out 2>&1', shell=True, capture_output=True)
            grompp2 = subp.run(f'gmx grompp -c {name}-em-sc.gro -f {resource_dir}/martini_v2.x_new-rf-em-test.mdp -p {resource_dir}/system.top -o topol.tpr -maxwarn 5', shell=True, capture_output=True)
            mdrun2 = subp.run(f'export OMP_NUM_THREADS=1; gmx mdrun -nt 1 -v -deffnm topol -c {name}-em.gro  > md.out 2>&1', shell=True, capture_output=True)

            if os.path.isfile(f"{name}-em.gro") == True:
                if os.path.isfile("md.out") == True:
                    # open md.out and get energy values
                    with open('md.out') as f:
                        for line in f:
                            if "Potential Energy" in line:
                                pe = float(line.rstrip().split()[3])
                            elif "Maximum force" in line:
                                mf = float(line.rstrip().split()[3])
                    LOGGER.debug(f"minimization results for {name}: potential energy of {pe} and maximum force of {mf}")
                    gromacs_file = f"{name}.gro"
                    if mf > -max_mf and mf < max_mf and pe > -max_pe and pe < max_pe:
                        valid_status = True
                        ## This is a healing section, that's the new addition compared to classic CGValidator class
                        gmx_output = subp.run(f'echo 0 | gmx trjconv -f {name}-em.gro -s topol.tpr -o {name}-em_whole.gro -pbc whole', shell=True, capture_output=True)
                        whole_gro_file = f"{name}-em_whole.gro"
                        u = mda.Universe(whole_gro_file)
                        positions_new = u.atoms.positions
                        positions_new[:,2] = positions_new[:,2]  - positions_new[413,2]
                        positions_new[:,1] = positions_new[:,1]  - positions_new[413,1]

                        u.atoms.positions = positions_new
                        u.atoms.write(whole_gro_file)
                        subp.call(f'mv {name}-em_whole.gro {cwd}/{name}.gro', shell=True)
                        # We also save the topology for the selection process
                        subp.call(f'mv topol.tpr {cwd}/{name}.tpr', shell=True)
                        end = timer() - start
                        LOGGER.info(f"Structure {name} VALID (potential energy of {pe} and maximum force of {mf}) in {end:.3f} seconds")
                    else:
                        end = timer() - start
                        LOGGER.warning(f"Structure {name} INVALID in {end:.3f} seconds")
                        if mf >= max_mf:
                            LOGGER.warning(f"  > maximum force of {mf} higher than threshold {max_mf}")
                        elif pe >= max_pe:
                            LOGGER.warning(f"  > potential energy  of {pe} higher than threshold {max_pe}")
                        else:
                            LOGGER.warning(f"  > Unknown problem {gmx_output}")
                else:
                    LOGGER.error(f"ERROR in GROMACS minimization of {name} - file md.out not found => code={grompp2.returncode} stderr={grompp2.stderr}")
            else:
                LOGGER.error(f"ERROR in GROMACS minimization of {name} - file {name}-em.gro not found =>  code={grompp.returncode} stderr={grompp.stderr}")

        except Exception as error:
            gromacs_file = None
            valid_status = False
            LOGGER.error(f"ERROR in GROMACS minimization of {name} - error was: {error}")

        finally:
            os.chdir(old_cwd)
            if cleanup:
                rmtree(fullPathDir)
        return [valid_status, gromacs_file, mf, pe]