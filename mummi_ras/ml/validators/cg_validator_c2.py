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

import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.rms import rmsd

from .cg_validator import CGBaseValidator

LOGGER = logging.getLogger(__name__)

class CGValidator_C2(CGBaseValidator):
    """
    Class to validate the generated CG protein structures for C2.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def write_validation_info(self, iteration_id: int, return_array: np.array, all_structure_names: np.array):
        iteration_id = "{:02d}".format(iteration_id)
        RPATH = os.path.join(self.iteration_path, f"iter{iteration_id}")
        self.current_rpath = RPATH
        if not os.path.isdir(RPATH):
            os.makedirs(RPATH, exist_ok=True)

        if len(return_array) == [] or all_structure_names == []:
            return [], []

        npatches = return_array.shape[0]
        files = np.empty(npatches, dtype='object')
        master_file = os.path.join(RPATH, f"master_iter{iteration_id}.npz")

        # Find valid patches
        valid_ids = np.where(return_array[:,0])[0]
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
            try:
                maximum_force = np.append(previous_data['maximum_force'], maximum_force)
                potential_energy = np.append(previous_data['potential_energy'], potential_energy)
            except KeyError as e:
                maximum_force = return_array[:,2]
                potential_energy = return_array[:,3]
                pass

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

    # Two functions from Tim C. to check structure and pick what FF to use
    # Note cutoff values should be in Angstroms
    # @TODO this analysis will also go into cgAnalysis and this test will also be used in createsims consolidate code
    @staticmethod
    def zn_cutoff(u: mda.Universe, cutoff: float, name: str) -> bool:
        ZN740 = u.select_atoms('name ZN and index 1805') #and resid 740
        CYS280 = u.select_atoms('name SC1 and index 750') #and resid 280
        ZN_distances = mda.analysis.distances.distance_array(ZN740.atoms.positions,CYS280.atoms.positions)
        if ZN_distances < cutoff:
            LOGGER.info("Protein structure {} has ZN740 with {} of CYS280 less than cutoff {} use -DACYM280 parameters".format(name, ZN_distances, cutoff))
            return True
        else:
            LOGGER.info("Protein structure {} has ZN740 with {} of CYS280 more than cutoff {} dont use -DACYM280 parameters".format(name, ZN_distances, cutoff))
            return False

    @staticmethod
    def SEP365_cutoff(u: mda.Universe, cutoff: float, name: str) -> bool:
        SEP365 = u.select_atoms('resname SEP and index 942 and name SC1')
        SEP365_ARG56 = u.select_atoms('resname ARG and index 1943 and name SC2') # and resid 56
        SEP365_ARG127 = u.select_atoms('resname ARG and index 2116 and name SC2') #and resid 127
        SEP_distances1 = mda.analysis.distances.distance_array(SEP365.atoms.positions,SEP365_ARG56.atoms.positions)
        SEP_distances2 = mda.analysis.distances.distance_array(SEP365.atoms.positions,SEP365_ARG127.atoms.positions)
        if SEP_distances1 > cutoff or SEP_distances2 > cutoff:
            LOGGER.info("Protein structure {} has SEP365 {} and {} from ARG56 and 127, more than cutoff {} use -DASEP365 parameters".format(name, SEP_distances1, SEP_distances2, cutoff))
            return True
        else:
            LOGGER.info("Protein structure {} has SEP365 {} and {} from ARG56 and 127, within cutoff {} dont use -DASEP365 parameters".format(name, SEP_distances1, SEP_distances2, cutoff))
            return False

    # Four check functions from Fikret to check structure for violation related to knot formation and
    # problematic loop configuration in CRD-Loop1, RAS peptide-14-3-3 and potential RBD-CRD overlap
    @staticmethod
    def check_knot_crd_loop1(u: mda.Universe, min_dist: float, name: str) -> bool:
        # Define three parts (p1, p2 and p3) in CRD and Loop1 that have potential to form a knot
        crdp1 = u.select_atoms('index 696 700 703 712') # resid 259 260 261 266 and name BB
        crdp2 = u.select_atoms('index 728 730 733 735 737 739 741') # resid 270:276 and name BB
        crdp3 = u.select_atoms('index 685 687') # resid 254 255 and name BB
        linkp1 = u.select_atoms('index 753 755 760 762') # resid 282:285 and name BB
        linkp2 = u.select_atoms('index 780 783 787 791') # resid 293:296 and name BB
        linkp3 = u.select_atoms('index 739') # resid 275 and name BB

        # Get center of masses for each part of CRD and Linker 1
        crdp1_com = crdp1.center_of_mass(compound='group')
        crdp2_com = crdp2.center_of_mass(compound='group')
        crdp3_com = crdp3.center_of_mass(compound='group')
        linkp1_com = linkp1.center_of_mass(compound='group')
        linkp2_com = linkp2.center_of_mass(compound='group')
        linkp3_com = linkp3.center_of_mass(compound='group')

        # Calculate the distances between COMs of corresponding parts
        linkp1_crdp1_dist = distances.distance_array(linkp1_com, crdp1_com,
                                        box=u.dimensions)
        linkp2_crdp2_dist = distances.distance_array(linkp2_com, crdp2_com,
                                        box=u.dimensions)
        linkp3_crdp3_dist = distances.distance_array(linkp3_com, crdp3_com,
                                        box=u.dimensions)

        # Check if either of the distances is below minimum allowed distance to detect a knot
        if linkp1_crdp1_dist < min_dist or linkp2_crdp2_dist < min_dist or linkp3_crdp3_dist < min_dist:
            LOGGER.info(f"Protein structure {name} has a potential knot formation as at least "
                f"one distance between parts of CRD and Loop 1 (d1: {linkp1_crdp1_dist} "
                f"d2: {linkp2_crdp2_dist} d3: {linkp2_crdp2_dist}) is below minimum "
                f"allowed distance {min_dist}")
            return False
        else:
            LOGGER.info(f"Protein structure {name} doesn't have a knot in CRD-Loop 1 region "
                f"Distances between parts of CRD and Loop 1 are d1:  {linkp1_crdp1_dist} "
                f"d2: {linkp2_crdp2_dist}, and d3: {linkp3_crdp3_dist} which are above "
                f"minimum allowed distance {min_dist}")
            return True

    @staticmethod
    def check_knot_hvr_rbd(u: mda.Universe, min_dist: float, name: str) -> bool:
        # Define a part (p1) in HVR and RBD that have potential to form a knot
        hvrp1 = u.select_atoms('index 417 420 422 425') # resid 180 181 182 183 and name BB
        rbdp1 = u.select_atoms('index 454 459 483 499 549 580 599 603') # resid 157 159 169 177 199 212 220 222 and name BB

        # Get center of masses for each part of HVR and RBD
        hvrp1_com = hvrp1.center_of_mass(compound='group')
        rbdp1_com = rbdp1.center_of_mass(compound='group')

        # Calculate the distances between COMs of corresponding parts
        hvrp1_rbdp1_dist = distances.distance_array(hvrp1_com, rbdp1_com,
                                                    box=u.dimensions)

        # Check if the distance is below minimum allowed distance to detect a knot
        if hvrp1_rbdp1_dist < min_dist:
            LOGGER.info(f"Protein structure {name} has a potential knot formation as "
                f"the distance between parts of HVR and RBD ({hvrp1_rbdp1_dist}) "
                f"allowed distance {min_dist}")
            return False
        else:
            LOGGER.info(f"Protein structure {name} doesn't have a knot in HVR-RBD region "
                f"Distance between parts of HVR and RBD is {hvrp1_rbdp1_dist} "
                f"which is above minimum allowed distance {min_dist}")
            return True

    @staticmethod
    def check_rbd_crd_distance(u: mda.Universe, min_dist: float, name: str) -> bool:
        # Define RBD and CRD
        crd = u.select_atoms('name BB and index 637 to 740')
        rbd = u.select_atoms('name BB and index 452 to 618')

        # Get center of masses for RBD and CRD
        crd_com = crd.center_of_mass(compound='group')
        rbd_com = rbd.center_of_mass(compound='group')

        # Calculate the distances between COMs of RBD and CRD
        rbd_crd_dist = distances.distance_array(rbd_com, crd_com,
                                            box=u.dimensions)

        # Check if the COM distances is below minimum allowed distance to detect a bad configuration
        if rbd_crd_dist < min_dist:
            LOGGER.info(f"Protein structure {name} has a potentially bad RBD-CRD configuration "
                f"the COM distance {rbd_crd_dist} is below minimum allowed distance {min_dist}")
            return False
        else:
            LOGGER.info(f"Protein structure {name} has a proper RBD-CRD distance "
                f"the COM distance {rbd_crd_dist} is above minimum allowed distance {min_dist}")
            return True

    @staticmethod
    def check_peptide_1433_distance(u: mda.Universe, min_dist: float, name: str) -> bool:
        # Define RAF peptide and 14-3-3 helix regions
        peptide_rafp1 = u.select_atoms('index 1778 1780 1782') # resid 727 728 729 and name BB
        helices_1433p1 = u.select_atoms('index 2756 2758 2814 2816 2864 2866 2868 2870 2872 2874') # resid 157 159 169 177 199 212 220 222 and name BB

        peptide_rafp2 = u.select_atoms('index 1797 1799 1800') # resid 736 737 738 and name BB
        helices_1433p2 = u.select_atoms('index 2616 2619 2630 2640 2708') # resid 115 116 119 123 151 and name BB

        peptide_rafp1_com = peptide_rafp1.center_of_mass(compound='group')
        helices_1433p1_com = helices_1433p1.center_of_mass(compound='group')

        peptide_rafp2_com = peptide_rafp2.center_of_mass(compound='group')
        helices_1433p2_com = helices_1433p2.center_of_mass(compound='group')

        peptiderafp1_helices1433p1_dist = distances.distance_array(peptide_rafp1_com, helices_1433p1_com,
                                                            box=u.dimensions)

        peptiderafp2_helices1433p2_dist = distances.distance_array(peptide_rafp2_com, helices_1433p2_com,
                                                            box=u.dimensions)

        # Check if the COM distanes is below minimum allowed distance to detect a bad configuration
        if peptiderafp1_helices1433p1_dist < min_dist or peptiderafp2_helices1433p2_dist < min_dist:
            LOGGER.info(f"Protein structure {name} has a potential issue as at least "
                f"one distance between parts of RAF peptide and 14-3-3 helices "
                f"(d1: {peptiderafp1_helices1433p1_dist} d2: {peptiderafp2_helices1433p2_dist}) "
                f"is below minimum allowed distance {min_dist}")
            return False
        else:
            LOGGER.info(f"Protein structure {name} doesn't have an issue in RAF peptide - 14-3-3 helix region "
                f"Distances between parts of RAF peptide and 14-3-3 helices are "
                f"d1: {peptiderafp1_helices1433p1_dist} d2: {peptiderafp2_helices1433p2_dist} "
                f"is above minimum allowed distance {min_dist}")
            return True

    @staticmethod
    def RBD_rmsd(u: mda.Universe, ref: mda.Universe, rmsd_cutoff: float, dist_cutoff: float, name: str) -> bool:
        # Looks to see if the RBD rmsd is close to the reference, which would indicate that the fold is correct
        # Also checks that the domain is not mirrored
        # Also need to check that the RBD-CRD loop does not thread back inside the domain
        RBD = u.select_atoms('name BB and index 452 to 618')
        #refu = mda.Universe(str(ref))
        RBDref = ref.select_atoms('name BB and index 452 to 618')
        RBD_CRD_loop = u.select_atoms('name BB and index 619 to 636')
        RBD_com = RBD.atoms.center_of_geometry()
        loop_com = RBD_CRD_loop.atoms.center_of_geometry()
        dist = mda.analysis.distances.distance_array(RBD_com,loop_com)
        RMSD = rmsd(RBD.positions,RBDref.positions,superposition=True)

        # Check if the RMSD is below the threshold to determine that the loop is correct
        # and check if the loop COM is far enough away from the RBD that it is not entangled
        if RMSD <= rmsd_cutoff and dist > dist_cutoff:
            LOGGER.info("Protein structure {} has a proper RBD structure".format(name))
            LOGGER.info("The RMSD {} is below the threshold {}".format(RMSD, rmsd_cutoff))
            LOGGER.info("The RBD loop is {} away from the RBD COM which is above the threshold {}".format(dist[0][0], dist_cutoff))
            return True
        elif RMSD > rmsd_cutoff:
            LOGGER.info("Protein structure {} has a potentially bad RBD structure".format(name))
            LOGGER.info("The RMSD {} is above the threshold {}".format(RMSD, rmsd_cutoff))
            return False
        elif dist < dist_cutoff:
            LOGGER.info("Protein structure {} has a potentially bad RBD structure".format(name))
            LOGGER.info("The RBD loop is {} away from the RBD COM which is below the threshold {}".format(dist[0][0], dist_cutoff))
            return False

    @staticmethod
    def CRD_rmsd(u: mda.Universe, ref: mda.Universe, rmsd_cutoff: float, name: str) -> bool:
        # The RMSD checks that the domain is not mirrored
        CRD =   u.select_atoms('name BB and index 637 to 740')
        #refu = mda.Universe(str(ref))
        CRDref = ref.select_atoms('name BB and index 637 to 740')
        RMSD = rmsd(CRD.positions,CRDref.positions,superposition=True)
        if RMSD <= rmsd_cutoff:
            LOGGER.info("Protein structure {} has a proper CRD structure".format(name))
            LOGGER.info("The RMSD {} is below the threshold {}".format(RMSD, rmsd_cutoff))
            return True
        else:
            LOGGER.info("Protein structure {} has a potentially bad CRD structure".format(name))
            LOGGER.info("The RMSD {} is above the threshold {}".format(RMSD, rmsd_cutoff))
            return False

    @staticmethod
    def Gdom_rmsd(u: mda.Universe, ref: mda.Universe, rmsd_cutoff: float, dist_cutoff: float, name: str) -> bool:
        # Looks to see if the G-domain rmsd is close to the reference, which would indicate that the fold is correct
        # Also checks that the domain is not mirrored
        # Also need to check that the HVR does not thread back inside the domain
        Gdom =   u.select_atoms('name BB and index 0 to 387')
        #refu = mda.Universe(str(ref))
        Gdomref = ref.select_atoms('name BB and index 0 to 387')
        HVRstart = u.select_atoms('name BB and index 388 to 420')
        Gdom_com = Gdom.atoms.center_of_geometry()
        HVRstart_com = HVRstart.atoms.center_of_geometry()
        dist = mda.analysis.distances.distance_array(Gdom_com,HVRstart_com)
        RMSD = rmsd(Gdom.positions,Gdomref.positions,superposition=True)
        if RMSD <= rmsd_cutoff and dist > dist_cutoff:
            LOGGER.info("Protein structure {} has a proper G-domain structure".format(name))
            LOGGER.info("The RMSD {} is below the threshold {}".format(RMSD, rmsd_cutoff))
            LOGGER.info("The HVR is {} away from the G-domain COM which is above the threshold {}".format(dist[0][0], dist_cutoff))
            return True
        elif RMSD > rmsd_cutoff:
            LOGGER.info("Protein structure {} has a potentially bad G-domain structure".format(name))
            LOGGER.info("The RMSD {} is above the threshold {}".format(RMSD, rmsd_cutoff))
            return False
        elif dist < dist_cutoff:
            LOGGER.info("Protein structure {} has a potentially bad G-domain structure".format(name))
            LOGGER.info("The HVR is {} away from the G-domain COM which is below the threshold {}".format(dist[0][0], dist_cutoff))
            return False

    def _internal_validate(self, data: List[Any]) -> List:
        name, positions, cwd, resource_dir, complex_name, max_mf, max_pe, healing, cleanup, mini_mummi = data
        gromacs_file = None
        valid_status = False
        timeout_gmx = 60 # After 1 mins we kill GROMACS
        start = timer()

        # Make temp dir to run minimization in
        fullPathDir = os.path.join(cwd, name)
        gromacs_file = f"{name}.gro"
        old_cwd = os.getcwd()
        if not os.path.isdir(fullPathDir):
            os.makedirs(fullPathDir, exist_ok=True)

        try:
            os.chdir(fullPathDir)
        except Exception as e:
            LOGGER.error(f"chdir failed {fullPathDir} ({e.__class__.__name__} {e})")

        LOGGER.debug(f"Switched dir to fullPathDir {fullPathDir} ({os.getcwd()})")

        # Save .gro file - @CHECK faster if MDAnalysis universe is in class and here only replace and save BUT maybe prolematic for later threaded version ???
        gromacs_ref_file = os.path.join(resource_dir, complex_name)
        pe = np.inf # Potential energy
        mf = np.inf # Maximum force

        try:
            u = mda.Universe(gromacs_ref_file)
            u.atoms.positions = positions

            if not mini_mummi:
                ref = u.copy()
                if (CGValidator_C2.check_knot_crd_loop1(u, 8, name) and
                    CGValidator_C2.check_knot_hvr_rbd(u, 8, name) and
                    CGValidator_C2.check_peptide_1433_distance(u, 6, name) and
                    CGValidator_C2.check_rbd_crd_distance(u, 20, name) and
                    CGValidator_C2.Gdom_rmsd(u,ref,4,17,name) and
                    CGValidator_C2.CRD_rmsd(u,ref,5.5,name) and
                    CGValidator_C2.RBD_rmsd(u,ref,2,11,name)):
                    LOGGER.info(f"Structure {name} passed the CRD-Loop1 knot, HVR-RBD knot, "
                        f"RAF peptide-14-3-3 helix distance, RBD-CRD COM distance checks, "
                        f"G-domain and HVR checks, CRD-mirror checks, and RBD-mirror and mis-threading checks"
                        f"Proceed with writing the gro file")
                else:
                    end = timer() - start
                    LOGGER.info(f"Structure {name} failed at one of the checks in {end:.3f} seconds")
                    if cleanup:
                        self._cleanup_dir(old_cwd, fullPathDir)
                    return [False, gromacs_file, -1, -1]

            u.atoms.write(f"{name}.gro")
        except ValueError as e:
            end = timer() - start
            LOGGER.error(f"Writing atoms {name} pos=[min={positions.min()}, max={positions.max()}]: {e} (in {end:.3f} seconds)")
            if cleanup:
                self._cleanup_dir(old_cwd, fullPathDir)
            return [False, gromacs_file, -1, -1]

        # Check structure for needed FF changes
        added_mdp_options = ""
        if mini_mummi:
            added_mdp_options += " -DACYM280"
        else:
            if CGValidator_C2.zn_cutoff(u, 8, name):
                added_mdp_options += " -DACYM280"
            if CGValidator_C2.SEP365_cutoff(u, 8, name):
                added_mdp_options += " -DASEP365"

        subp.call(f"cp {resource_dir}/martini_v2.x_new-rf-em-sc-test.mdp .", shell=True)
        subp.call(f"cp {resource_dir}/martini_v2.x_new-rf-em-test.mdp .", shell=True)
        subp.call(f"sed -ie 's/XXX/{added_mdp_options}/g' martini_v2.x_new-rf-em-test.mdp", shell=True)
        subp.call(f"sed -ie 's/XXX/{added_mdp_options}/g' martini_v2.x_new-rf-em-sc-test.mdp", shell=True)

        # Run minimization
        try:
            LOGGER.info(f"Run GROMACS minimization in dir {fullPathDir}/{name}")
            try:
                # Run energy minimization (steepest descent, with protein constreins)
                grompp = subp.run(
                    f'gmx grompp -c {name}.gro -f martini_v2.x_new-rf-em-sc-test.mdp -p {resource_dir}/system-sc.top -o topol-sc.tpr -maxwarn 5',
                    shell=True,
                    capture_output=True,
                    timeout=timeout_gmx
                )
                # @NOTE: some gmx versions are compiled without thread-MPI and do not support setting the number of threads "-nt" so fix OMP_NUM....
                mdrun = subp.run(
                    f'export OMP_NUM_THREADS=1; gmx mdrun -nt 1 -ntmpi 1 -v -deffnm topol-sc -c {name}-em-sc.gro  > md-sc.out 2>&1',
                    shell=True,
                    capture_output=True,
                    timeout=timeout_gmx
                )
                grompp2 = subp.run(
                    f'gmx grompp -c {name}-em-sc.gro -f martini_v2.x_new-rf-em-test.mdp -p {resource_dir}/system.top -o topol.tpr -maxwarn 5',
                    shell=True,
                    capture_output=True,
                    timeout=timeout_gmx
                )
                mdrun2 = subp.run(
                    f'export OMP_NUM_THREADS=1; gmx mdrun -nt 1 -ntmpi 1 -v -deffnm topol -c {name}-em.gro  > md.out 2>&1',
                    shell=True,
                    capture_output=True,
                    timeout=timeout_gmx
                )
            except subp.TimeoutExpired as e:
                end = timer() - start
                LOGGER.error(f"GROMACS timeout for {name}: {e} (in {end:.3f} seconds)")
                if cleanup:
                    self._cleanup_dir(old_cwd, fullPathDir)
                return [False, gromacs_file, -1, -1]

            if os.path.isfile(f"{name}-em.gro"):
                if os.path.isfile("md.out"):
                    # open md.out and get energy values
                    with open('md.out') as f:
                        for line in f:
                            if "Potential Energy" in line:
                                pe = float(line.rstrip().split()[3])
                            elif "Maximum force" in line:
                                mf = float(line.rstrip().split()[3])
                    LOGGER.debug(f"minimization results for {name}: potential energy of {pe} and maximum force of {mf}")
                    if mf > -max_mf and mf < max_mf and pe > -max_pe and pe < max_pe:
                        valid_status = True
                        if healing:
                            ## This is a healing section, that's the new addition compared to classic CGValidator_C2 class
                            try:
                                gmx_output = subp.run(f'echo 0 | gmx trjconv -f {name}-em.gro -s topol.tpr -o {name}-em_whole.gro -pbc whole', shell=True, capture_output=True, timeout=timeout_gmx)
                            except subp.TimeoutExpired as e:
                                end = timer() - start
                                LOGGER.error(f"GROMACS timeout for {name}: {e} (in {end:.3f} seconds)")
                                if cleanup:
                                    self._cleanup_dir(old_cwd, fullPathDir)
                                return [False, gromacs_file, -1, -1]

                            whole_gro_file = f"{name}-em_whole.gro"
                            try:
                                u = mda.Universe(whole_gro_file)
                                positions_new = u.atoms.positions
                                positions_new[:,2] = positions_new[:,2]  - positions_new[430,2]
                                positions_new[:,1] = positions_new[:,1]  - positions_new[430,1]
                                if positions_new[430,0] < 10:
                                    positions_new[:,0] += 500 # Box size
                                u.atoms.positions = positions_new
                                u.atoms.write(whole_gro_file)
                            except ValueError as e:
                                end = timer() - start
                                LOGGER.error(f"Validating {name} (whole_gro_file) in {end:.3f} seconds pos=[min={positions_new.min()}, max={positions_new.max()}]: {e}")
                                if cleanup:
                                    self._cleanup_dir(old_cwd, fullPathDir)
                                return [False, gromacs_file, -1, -1]

                            subp.call(f'mv {name}-em_whole.gro {cwd}/{name}.gro', shell=True)
                            # # We also save the topology for the selection process (only needed for UCG)
                            # subp.call(f'mv topol.tpr {cwd}/{name}.tpr', shell=True)
                        else:
                            subp.call(f'mv {name}.gro {cwd}/', shell=True)
                        end = timer() - start
                        LOGGER.info(f"Structure {name} VALID (potential energy of {pe}, maximum force of {mf}) in {end:.3f} seconds")
                    else:
                        end = timer() - start
                        LOGGER.info(f"Structure {name} INVALID in {end:.3f} seconds")
                        if mf >= max_mf:
                            LOGGER.warning(f"  > maximum force of {mf} higher than threshold {max_mf}")
                        elif pe >= max_pe:
                            LOGGER.warning(f"  > potential energy  of {pe} higher than threshold {max_pe}")
                        else:
                            LOGGER.warning(f"  > Unknown problem {gmx_output}")
                else:
                    end = timer() - start
                    LOGGER.error(f"ERROR in GROMACS minimization of {name} - file md.out not found")
                    print(grompp2.stderr.decode())
                    print(grompp2.stdout.decode())
            else:
                LOGGER.error(f"ERROR in GROMACS minimization of {name} - file {name}-em.gro not found")
                print(grompp.stderr.decode())
                print(grompp.stdout.decode())
        except Exception as error:
            gromacs_file = None
            valid_status = False
            LOGGER.error(f"ERROR in GROMACS minimization of {name} - error was: {error.__class__.__name__}: {error}")
        finally:
            if cleanup:
                self._cleanup_dir(old_cwd, fullPathDir)
        return [valid_status, gromacs_file, mf, pe]
