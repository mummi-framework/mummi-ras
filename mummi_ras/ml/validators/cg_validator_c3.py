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
import math
import subprocess as subp
import multiprocessing
from typing import List, Any
from timeit import default_timer as timer
from shutil import rmtree, which
import signal

import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis import transformations

from .cg_validator import CGBaseValidator

import contextlib

LOGGER = logging.getLogger(__name__)

class CGValidator_C3(CGBaseValidator):
    """
    Class to validate the generated CG protein structures of C3.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_status(self, return_array: np.array) -> np.array:
        """"
        Return the status of structures based on the result from validate()
        """
        if return_array.shape[0] == 0:
            return np.array([])
        return np.array([ff for struct in return_array[:,0] for ff in struct])

    def get_structure_names(self, return_array: np.array) -> np.array:
        """"
        Return the structure names of result from validate()
        """
        if return_array.shape[0] == 0:
            return np.array([])
        return np.array([ff for struct in return_array[:,1] for ff in struct])

    def get_valid_files(self, return_array: np.array) -> np.array:
        """"
        REturn the list of valid files based on the result from validate()
        """
        if return_array.shape[0] == 0:
            return np.array([])
        status = self.get_status(return_array)
        gromacs_files = self.get_structure_names(return_array)
        valid_ids = np.where(status)[0]
        return gromacs_files[valid_ids]

    def write_validation_info(self, iteration_id: int, return_array: np.array, all_structure_names: np.array):
        iteration_id = "{:02d}".format(iteration_id)
        RPATH = os.path.join(self.iteration_path, f"iter{iteration_id}")
        self.current_rpath = RPATH
        if not os.path.isdir(RPATH):
            os.makedirs(RPATH, exist_ok=True)

        if len(return_array) == 0 or all_structure_names == []:
            return [], []

        master_file = os.path.join(RPATH, f"master_iter{iteration_id}.npz")

        # We need to flatten the input array

        status = np.array([ff for struct in return_array[:,0] for ff in struct])
        gromacs_files = np.array([ff for struct in return_array[:,1] for ff in struct])
        maximum_force = np.array([ff for struct in return_array[:,2] for ff in struct])
        potential_energy = np.array([ff for struct in return_array[:,3] for ff in struct])

        npatches = status.shape[0]

        # Find valid patches
        valid_ids = np.where(status)[0]
        valid_files = gromacs_files[valid_ids]

        files = np.array(valid_files)

        # we do not erase the previous master file
        if os.path.isfile(master_file):
            previous_data = np.load(master_file, allow_pickle=True)["files"]
            files = np.append(previous_data, files)

        # # Save the list of the files into the master file
        np.savez_compressed(master_file, files = files)

        # # Save the details for all the predictions
        name_table_all_predictions = os.path.join(os.path.join(RPATH, "table_all_predictions.npz"))

        if os.path.isfile(name_table_all_predictions):
            previous_data = np.load(name_table_all_predictions, allow_pickle=True)
            # Here structure_names and gromacs_files are identical
            # WE could remove gromacs_files. For now we keep both for legacy reasons
            all_structure_names = np.append(previous_data['structure_names'], gromacs_files)
            status = np.append(previous_data['status'], status)
            gromacs_files = np.append(previous_data['gromacs_files'], gromacs_files)
            try:
                maximum_force = np.append(previous_data['maximum_force'], maximum_force)
                potential_energy = np.append(previous_data['potential_energy'], potential_energy)
            except KeyError as e:
                maximum_force = np.array([ff for struct in return_array[:,2] for ff in struct])
                potential_energy = np.array([ff for struct in return_array[:,3] for ff in struct])
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

        success_rate = 100 * valid_ids.shape[0] / npatches
        LOGGER.info(f'Found {valid_ids.shape[0]} valid patches out of {npatches} patches (success rate = {success_rate:.2f} %)')
        return valid_ids, valid_files

    # Two functions from Tim C. to check structure and pick what FF to use
    # Note cutoff values should be in Angstroms
    # @TODO this analysis will also go into cgAnalysis and this test will also be used in createsims consolidate code
    @staticmethod
    def zn_cutoff(u: mda.Universe, cutoff: float, name: str) -> bool:
        ZN740 = u.select_atoms('name ZN and index 1805') #and resid 740
        CYS280 = u.select_atoms('name SC1 and index 750') #and resid 280

        ZN740_2 = u.select_atoms('name ZN and index 4702')#and resid 740
        CYS280_2 = u.select_atoms('name SC1 and index 3647')#and resid 280

        ZN_distances = mda.analysis.distances.distance_array(ZN740.atoms.positions,CYS280.atoms.positions)
        ZN_distances_2 = mda.analysis.distances.distance_array(ZN740_2.atoms.positions,CYS280_2.atoms.positions)
        # if ZN_distances < cutoff:
        #     LOGGER.debug("Protein structure {} has ZN740 with {} of CYS280 less than cutoff {} use -DACYM280 parameters".format(name, ZN_distances, cutoff))
        #     if ZN_distances_2 < cutoff:
        #         LOGGER.info("Protein structure {} has ZN740 with {} of CYS280 less than cutoff {} use -DACYM280 parameters (2nd monomer)".format(name, ZN_distances_2, cutoff))
        #         return True
        #     else:
        #         LOGGER.info("Protein structure {} has ZN740 with {} of CYS280 more than cutoff {} dont use -DACYM280 parameters (2nd monomer)".format(name, ZN_distances_2, cutoff))
        #         return False
        # else:
        #     LOGGER.debug("Protein structure {} has ZN740 with {} of CYS280 more than cutoff {} dont use -DACYM280 parameters".format(name, ZN_distances, cutoff))
        #     return False
        return True

    @staticmethod
    def SEP365_cutoff(u: mda.Universe, cutoff: float, name: str) -> bool:
        SEP365 = u.select_atoms('resname SEP and index 942 and name SC1')
        SEP365_ARG56 = u.select_atoms('resname ARG and index 1943 and name SC2') # and resid 56
        SEP365_ARG127 = u.select_atoms('resname ARG and index 2116 and name SC2') #and resid 127
        # SEP_distances1 = mda.analysis.distances.distance_array(SEP365.atoms.positions,SEP365_ARG56.atoms.positions)
        # SEP_distances2 = mda.analysis.distances.distance_array(SEP365.atoms.positions,SEP365_ARG127.atoms.positions)
        # if SEP_distances1 > cutoff or SEP_distances2 > cutoff:
        #     LOGGER.debug("Protein structure {} has SEP365 {} and {} from ARG56 and 127, more than cutoff {} use -DASEP365 parameters".format(name, SEP_distances1, SEP_distances2, cutoff))
        #     return True
        # else:
        #     LOGGER.debug("Protein structure {} has SEP365 {} and {} from ARG56 and 127, within cutoff {} dont use -DASEP365 parameters".format(name, SEP_distances1, SEP_distances2, cutoff))
        #     return False
        return True

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

        crdp1_2 = u.select_atoms('index 3593 4297 3600 3609') # resid 259 260 261 266 and name BB
        crdp2_2 = u.select_atoms('index 3625 3627 3630 3632 3634 3636 3638') # resid 270:276 and name BB
        crdp3_2 = u.select_atoms('index 3582 3584') # resid 254 255 and name BB
        linkp1_2 = u.select_atoms('index 3650 3652 3657 3659') # resid 282:285 and name BB
        linkp2_2 = u.select_atoms('index 3677 3680 3684 3688') # resid 293:296 and name BB
        linkp3_2 = u.select_atoms('index 3636') # resid 275 and name BB

        # Get center of masses for each part of CRD and Linker 1
        crdp1_com = crdp1.center_of_mass(compound='group')
        crdp2_com = crdp2.center_of_mass(compound='group')
        crdp3_com = crdp3.center_of_mass(compound='group')
        linkp1_com = linkp1.center_of_mass(compound='group')
        linkp2_com = linkp2.center_of_mass(compound='group')
        linkp3_com = linkp3.center_of_mass(compound='group')

        crdp1_com_2 = crdp1_2.center_of_mass(compound='group')
        crdp2_com_2 = crdp2_2.center_of_mass(compound='group')
        crdp3_com_2 = crdp3_2.center_of_mass(compound='group')
        linkp1_com_2 = linkp1_2.center_of_mass(compound='group')
        linkp2_com_2 = linkp2_2.center_of_mass(compound='group')
        linkp3_com_2 = linkp3_2.center_of_mass(compound='group')

        # Calculate the distances between COMs of corresponding parts
        linkp1_crdp1_dist = distances.distance_array(linkp1_com, crdp1_com,
                                        box=u.dimensions)
        linkp2_crdp2_dist = distances.distance_array(linkp2_com, crdp2_com,
                                        box=u.dimensions)
        linkp3_crdp3_dist = distances.distance_array(linkp3_com, crdp3_com,
                                        box=u.dimensions)

        linkp1_crdp1_dist_2 = distances.distance_array(linkp1_com_2, crdp1_com_2,
                                    box=u.dimensions)
        linkp2_crdp2_dist_2 = distances.distance_array(linkp2_com_2, crdp2_com_2,
                                        box=u.dimensions)
        linkp3_crdp3_dist_2 = distances.distance_array(linkp3_com_2, crdp3_com_2,
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
            if linkp1_crdp1_dist_2 < min_dist or linkp2_crdp2_dist_2 < min_dist or linkp3_crdp3_dist_2 < min_dist:
                LOGGER.info(f"Protein structure {name} has a potential knot formation as at least "
                    f"one distance between parts of CRD and Loop 1 (d1: {linkp1_crdp1_dist_2} "
                    f"d2: {linkp2_crdp2_dist_2} d3: {linkp2_crdp2_dist_2}) is below minimum "
                    f"allowed distance {min_dist} (2nd monomer)")
                return False
            else:
                LOGGER.info(f"Protein structure {name} doesn't have a knot in CRD-Loop 1 region "
                    f"Distances between parts of CRD and Loop 1 are d1:  {linkp1_crdp1_dist_2} "
                    f"d2: {linkp2_crdp2_dist_2}, and d3: {linkp3_crdp3_dist_2} which are above "
                    f"minimum allowed distance {min_dist} (2nd monomer)")
                return True

    @staticmethod
    def check_knot_hvr_rbd(u: mda.Universe, min_dist: float, name: str) -> bool:
        # Define a part (p1) in HVR and RBD that have potential to form a knot
        hvrp1 = u.select_atoms('index 417 420 422 425') # resid 180 181 182 183 and name BB
        rbdp1 = u.select_atoms('index 454 459 483 499 549 580 599 603') # resid 157 159 169 177 199 212 220 222 and name BB

        hvrp1_2 = u.select_atoms('index 3314 3317 3319 3322') # resid 180 181 182 183 and name BB
        rbdp1_2 = u.select_atoms('index 3351 3356 3380 3396 3446 3477 3496 3500') # resid 157 159 169 17

        # Get center of masses for each part of HVR and RBD
        hvrp1_com = hvrp1.center_of_mass(compound='group')
        rbdp1_com = rbdp1.center_of_mass(compound='group')

        hvrp1_com_2 = hvrp1_2.center_of_mass(compound='group')
        rbdp1_com_2 = rbdp1_2.center_of_mass(compound='group')

        # Calculate the distances between COMs of corresponding parts
        hvrp1_rbdp1_dist = distances.distance_array(hvrp1_com, rbdp1_com,
                                                    box=u.dimensions)
        hvrp1_rbdp1_dist_2 = distances.distance_array(hvrp1_com_2, rbdp1_com_2,
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
            if hvrp1_rbdp1_dist_2 < min_dist:
                LOGGER.info(f"Protein structure {name} has a potential knot formation as "
                    f"the distance between parts of HVR and RBD ({hvrp1_rbdp1_dist_2}) "
                    f"allowed distance {min_dist} (2nd monomer)")
                return False
            else:
                LOGGER.info(f"Protein structure {name} doesn't have a knot in HVR-RBD region "
                    f"Distance between parts of HVR and RBD is {hvrp1_rbdp1_dist_2} "
                    f"which is above minimum allowed distance {min_dist} (2nd monomer)")
                return True

    @staticmethod
    def check_rbd_crd_distance(u: mda.Universe, min_dist: float, name: str) -> bool:
        # Define RBD and CRD
        crd = u.select_atoms('name BB and index 637 to 740')
        rbd = u.select_atoms('name BB and index 452 to 618')

        crd_2 = u.select_atoms('name BB and index 3534 to 3637')
        rbd_2 = u.select_atoms('name BB and index 3349 to 3515')

        # Get center of masses for RBD and CRD
        crd_com = crd.center_of_mass(compound='group')
        rbd_com = rbd.center_of_mass(compound='group')

        crd_com_2 = crd_2.center_of_mass(compound='group')
        rbd_com_2 = rbd_2.center_of_mass(compound='group')

        # Calculate the distances between COMs of RBD and CRD
        rbd_crd_dist = distances.distance_array(rbd_com, crd_com,
                                            box=u.dimensions)
        rbd_crd_dist_2 = distances.distance_array(rbd_com_2, crd_com_2,
                                        box=u.dimensions)

        # Check if the COM distances is below minimum allowed distance to detect a bad configuration
        if rbd_crd_dist < min_dist:
            LOGGER.info(f"Protein structure {name} has a potentially bad RBD-CRD configuration "
                f"the COM distance {rbd_crd_dist} is below minimum allowed distance {min_dist}")
            return False
        else:
            LOGGER.info(f"Protein structure {name} has a proper RBD-CRD distance "
                f"the COM distance {rbd_crd_dist} is above minimum allowed distance {min_dist}")
            if rbd_crd_dist_2 < min_dist:
                LOGGER.info(f"Protein structure {name} has a potentially bad RBD-CRD configuration "
                    f"the COM distance {rbd_crd_dist_2} is below minimum allowed distance {min_dist} (2nd monomer)")
                return False
            else:
                LOGGER.info(f"Protein structure {name} has a proper RBD-CRD distance "
                    f"the COM distance {rbd_crd_dist_2} is above minimum allowed distance {min_dist} (2nd monomer)")
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
            LOGGER.debug(f"Protein structure {name} has a potential issue as at least "
                f"one distance between parts of RAF peptide and 14-3-3 helices "
                f"(d1: {peptiderafp1_helices1433p1_dist} d2: {peptiderafp2_helices1433p2_dist}) "
                f"is below minimum allowed distance {min_dist}")
            return False
        else:
            LOGGER.debug(f"Protein structure {name} doesn't have an issue in RAF peptide - 14-3-3 helix region "
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
        RBD_2 = u.select_atoms('name BB and index 3349 to 3515')
        RBDref_2 = ref.select_atoms('name BB and index 3349 to 3515')
        RBD_CRD_loop_2 = u.select_atoms('name BB and index 3516 to 3533')

        RBD_com = RBD.atoms.center_of_geometry()
        loop_com = RBD_CRD_loop.atoms.center_of_geometry()
        dist = mda.analysis.distances.distance_array(RBD_com,loop_com)
        RMSD = rmsd(RBD.positions,RBDref.positions,superposition=True)

        RBD_com_2 = RBD_2.atoms.center_of_geometry()
        loop_com_2 = RBD_CRD_loop_2.atoms.center_of_geometry()
        dist_2 = mda.analysis.distances.distance_array(RBD_com_2,loop_com_2)
        RMSD_2 = rmsd(RBD_2.positions,RBDref_2.positions,superposition=True)

        # Check if the RMSD is below the threshold to determine that the loop is correct
        # and check if the loop COM is far enough away from the RBD that it is not entangled
        if RMSD <= rmsd_cutoff and dist > dist_cutoff:
            LOGGER.info("Protein structure {} has a proper RBD structure".format(name))
            LOGGER.info("The RMSD {} is below the threshold {}".format(RMSD, rmsd_cutoff))
            LOGGER.info("The RBD loop is {} away from the RBD COM which is above the threshold {}".format(dist[0][0], dist_cutoff))
            if RMSD_2 <= rmsd_cutoff and dist_2 > dist_cutoff:
                LOGGER.info("Protein structure {} has a proper RBD structure (2nd monomer)".format(name))
                LOGGER.info("The RMSD {} is below the threshold {} (2nd monomer)".format(RMSD_2, rmsd_cutoff))
                LOGGER.info("The RBD loop is {} away from the RBD COM which is above the threshold {} (2nd monomer)".format(dist_2[0][0], dist_cutoff))
                return True
            elif RMSD_2 > rmsd_cutoff:
                LOGGER.info("Protein structure {} has a potentially bad RBD structure (2nd monomer)".format(name))
                LOGGER.info("The RMSD {} is above the threshold {} (2nd monomer)".format(RMSD_2, rmsd_cutoff))
                return False
            elif dist_2 < dist_cutoff:
                LOGGER.info("Protein structure {} has a potentially bad RBD structure (2nd monomer)".format(name))
                LOGGER.info("The RBD loop is {} away from the RBD COM which is below the threshold {} (2nd monomer)".format(dist_2[0][0], dist_cutoff))
                return False
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
        CRD_2 =   u.select_atoms('name BB and index 3534 to 3637')
        #refu = mda.Universe(str(ref))
        CRDref = ref.select_atoms('name BB and index 637 to 740')
        CRDref_2 = ref.select_atoms('name BB and index 3534 to 3637')

        RMSD = rmsd(CRD.positions,CRDref.positions,superposition=True)
        RMSD_2 = rmsd(CRD_2.positions,CRDref_2.positions,superposition=True)
        if RMSD <= rmsd_cutoff:
            LOGGER.info("Protein structure {} has a proper CRD structure".format(name))
            LOGGER.info("The RMSD {} is below the threshold {}".format(RMSD, rmsd_cutoff))
            if RMSD_2 <= rmsd_cutoff:
                LOGGER.info("Protein structure {} has a proper CRD structure (2nd monomer)".format(name))
                LOGGER.info("The RMSD {} is below the threshold {} (2nd monomer)".format(RMSD_2, rmsd_cutoff))
                return True
            else:
                LOGGER.info("Protein structure {} has a potentially bad CRD structure (2nd monomer)".format(name))
                LOGGER.info("The RMSD {} is above the threshold {} (2nd monomer)".format(RMSD_2, rmsd_cutoff))
                return False
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

        Gdom_2 =   u.select_atoms('name BB and index 2897 to 3284')
        #refu = mda.Universe(str(ref))
        Gdomref_2 = ref.select_atoms('name BB and index 2897 to 3284')
        HVRstart_2 = u.select_atoms('name BB and index 3285 to 3317')
        Gdom_com_2 = Gdom_2.atoms.center_of_geometry()
        HVRstart_com_2 = HVRstart_2.atoms.center_of_geometry()
        dist_2 = mda.analysis.distances.distance_array(Gdom_com_2,HVRstart_com_2)
        RMSD_2 = rmsd(Gdom_2.positions,Gdomref_2.positions,superposition=True)
        if RMSD <= rmsd_cutoff and dist > dist_cutoff:
            LOGGER.info("Protein structure {} has a proper G-domain structure".format(name))
            LOGGER.info("The RMSD {} is below the threshold {}".format(RMSD, rmsd_cutoff))
            LOGGER.info("The HVR is {} away from the G-domain COM which is above the threshold {}".format(dist[0][0], dist_cutoff))
            if RMSD_2 <= rmsd_cutoff and dist_2 > dist_cutoff:
                LOGGER.info("Protein structure {} has a proper G-domain structure (2nd monomer)".format(name))
                LOGGER.info("The RMSD {} is below the threshold {} (2nd monomer)".format(RMSD_2, rmsd_cutoff))
                LOGGER.info("The HVR is {} away from the G-domain COM which is above the threshold {} (2nd monomer)".format(dist_2[0][0], dist_cutoff))
                return True
            elif RMSD_2 > rmsd_cutoff:
                LOGGER.info("Protein structure {} has a potentially bad G-domain structure (2nd monomer)".format(name))
                LOGGER.info("The RMSD {} is above the threshold {} (2nd monomer)".format(RMSD_2, rmsd_cutoff))
                return False
            elif dist_2 < dist_cutoff:
                LOGGER.info("Protein structure {} has a potentially bad G-domain structure (2nd monomer)".format(name))
                LOGGER.info("The HVR is {} away from the G-domain COM which is below the threshold {} (2nd monomer)".format(dist_2[0][0], dist_cutoff))
                return False
        elif RMSD > rmsd_cutoff:
            LOGGER.info("Protein structure {} has a potentially bad G-domain structure".format(name))
            LOGGER.info("The RMSD {} is above the threshold {}".format(RMSD, rmsd_cutoff))
            return False
        elif dist < dist_cutoff:
            LOGGER.info("Protein structure {} has a potentially bad G-domain structure".format(name))
            LOGGER.info("The HVR is {} away from the G-domain COM which is below the threshold {}".format(dist[0][0], dist_cutoff))
            return False

    @staticmethod
    def LB1m1_rmsd_dimer(u: mda.Universe, ref: mda.Universe, rmsd_cutoff: float, name: str) -> bool:
        # Looks to see if the lobe1 of monomer 1 is close to the dimer structure
        
        LB_m1 =   u.select_atoms('name BB and index 1134:1467')
        LB_m1_ref = ref.select_atoms('name BB and index 716:1049')
        
        RMSD = rmsd(LB_m1.positions,LB_m1_ref.positions,superposition=True)

        if RMSD <= rmsd_cutoff:
            LOGGER.debug("The RMSD {} is below the dimer threshold {}".format(RMSD, rmsd_cutoff))
            LOGGER.debug("Protein structure {} has a lobe 1 m1 domain close to the dimer structure".format(name))
            return True
        else:
            LOGGER.debug("The RMSD {} is above the dimer threshold {}".format(RMSD, rmsd_cutoff))
            LOGGER.debug("Protein structure {} has a lobe 1 m1 domain close to the monomer structure".format(name))
            return False

    @staticmethod
    def LB1m1_rmsd_monomer(u: mda.Universe, ref: mda.Universe, rmsd_cutoff: float, name: str) -> bool:
        # Looks to see if the lobe1 of monomer 1 is close to the monomer structure

        LB_m1 =   u.select_atoms('name BB and index 1134:1467')
        LB_m1_ref = ref.select_atoms('name BB and index 1134:1467')
        RMSD = rmsd(LB_m1.positions,LB_m1_ref.positions,superposition=True)

        if RMSD <= rmsd_cutoff:
            LOGGER.debug("The RMSD {} is below the monomer threshold {}".format(RMSD, rmsd_cutoff))
            LOGGER.debug("Protein structure {} has a lobe 1 m1 close domain to the monomer structure".format(name))
            return True
        else:
            LOGGER.debug("The RMSD {} is above the monomer threshold {}".format(RMSD, rmsd_cutoff))
            LOGGER.debug("Protein structure {} has a lobe 1 m1 close domain to the dimer structure".format(name))
            return False

    @staticmethod
    def LB1m2_rmsd_dimer(u: mda.Universe, ref: mda.Universe, rmsd_cutoff: float, name: str) -> bool:
        # Looks to see if the lobe1 of monomer 2 is close to the dimer structure
        LB_m2 =   u.select_atoms('name BB and index 4031:4364')
        LB_m2_ref = ref.select_atoms('name BB and index 17:350')

        RMSD = rmsd(LB_m2.positions,LB_m2_ref.positions,superposition=True)

        if RMSD <= rmsd_cutoff:
            LOGGER.debug("The RMSD {} is below the dimer threshold {}".format(RMSD, rmsd_cutoff))
            LOGGER.debug("Protein structure {} has a lobe 1 m2 domain close to the dimer structure".format(name))
            return True
        else:
            LOGGER.debug("The RMSD {} is above the dimer threshold {}".format(RMSD, rmsd_cutoff))
            LOGGER.debug("Protein structure {} has a lobe 1 m2 domain close to the monomer structure".format(name))
        return False

    @staticmethod
    def LB1m2_rmsd_monomer(u: mda.Universe, ref: mda.Universe, rmsd_cutoff: float, name: str) -> bool:
        # Looks to see if the lobe1 of monomer 2 is close to the dimer structure

        LB_m2 =   u.select_atoms('name BB and index 4031:4364')
        LB_m2_ref = ref.select_atoms('name BB and index 1134:1467')
        RMSD = rmsd(LB_m2.positions,LB_m2_ref.positions,superposition=True)

        if RMSD <= rmsd_cutoff:
            LOGGER.debug("The RMSD {} is below the monomer threshold {}".format(RMSD, rmsd_cutoff))
            LOGGER.debug("Protein structure {} has a lobe 1 m2 domain close to the monomer structure".format(name))
            return True
        else:
            LOGGER.debug("The RMSD {} is above the monomer threshold {}".format(RMSD, rmsd_cutoff))
            LOGGER.debug("Protein structure {} has a lobe 1 m2 domain close to the dimer structure".format(name))
            return False

    @staticmethod
    def m1m2E(u: mda.Universe, ref: mda.Universe, name: str) -> int:

        m1_kd_565ASP = u.select_atoms('name BB and index 1389')
        m1_kd_516PHE = u.select_atoms('name BB and index 1270')
        m1_kd_530GLN = u.select_atoms('name BB and index 1303')

        m2_kd_477HSE = u.select_atoms('name BB and index 4080')
        m2_kd_509ARG = u.select_atoms('name BB and index 4150')
        m2_kd_510HSD = u.select_atoms('name BB and index 4153')

        # Get center of masses for each part
        m1_kd_565ASP_com = m1_kd_565ASP.center_of_mass(compound='group')
        m1_kd_516PHE_com = m1_kd_516PHE.center_of_mass(compound='group')
        m1_kd_530GLN_com = m1_kd_530GLN.center_of_mass(compound='group')

        m2_kd_477HSE_com = m2_kd_477HSE.center_of_mass(compound='group')
        m2_kd_509ARG_com = m2_kd_509ARG.center_of_mass(compound='group')
        m2_kd_510HSD_com = m2_kd_510HSD.center_of_mass(compound='group')


        m1_kd_565ASP_m2_kd_477HSE_dist = distances.distance_array(m1_kd_565ASP_com, m2_kd_477HSE_com,
                                            box=u.dimensions)

        m1_kd_516PHE_m2_kd_509ARG_dist = distances.distance_array(m1_kd_516PHE_com, m2_kd_509ARG_com,
                                                    box=u.dimensions)

        m1_kd_530GLN_m2_kd_510HSD_dist = distances.distance_array(m1_kd_530GLN_com, m2_kd_510HSD_com,
                                                                box=u.dimensions)

        m1_kd_565ASP_ref = ref.select_atoms('name BB and index 971')
        m1_kd_516PHE_ref = ref.select_atoms('name BB and index 852')
        m1_kd_530GLN_ref = ref.select_atoms('name BB and index 885')

        m2_kd_477HSE_ref = ref.select_atoms('name BB and index 66')
        m2_kd_509ARG_ref = ref.select_atoms('name BB and index 136')
        m2_kd_510HSD_ref = ref.select_atoms('name BB and index 139')

        # Get center of masses for each part
        m1_kd_565ASP_com_ref = m1_kd_565ASP_ref.center_of_mass(compound='group')
        m1_kd_516PHE_com_ref = m1_kd_516PHE_ref.center_of_mass(compound='group')
        m1_kd_530GLN_com_ref = m1_kd_530GLN_ref.center_of_mass(compound='group')

        m2_kd_477HSE_com_ref = m2_kd_477HSE_ref.center_of_mass(compound='group')
        m2_kd_509ARG_com_ref = m2_kd_509ARG_ref.center_of_mass(compound='group')
        m2_kd_510HSD_com_ref = m2_kd_510HSD_ref.center_of_mass(compound='group')

        m1_kd_565ASP_m2_kd_477HSE_dist_ref = distances.distance_array(m1_kd_565ASP_com_ref, m2_kd_477HSE_com_ref,
                                            box=u.dimensions)

        m1_kd_516PHE_m2_kd_509ARG_dist_ref = distances.distance_array(m1_kd_516PHE_com_ref, m2_kd_509ARG_com_ref,
                                                    box=u.dimensions)

        m1_kd_530GLN_m2_kd_510HSD_dist_ref = distances.distance_array(m1_kd_530GLN_com_ref, m2_kd_510HSD_com_ref,
                                                                box=u.dimensions)


        diff_dist1 = m1_kd_565ASP_m2_kd_477HSE_dist - m1_kd_565ASP_m2_kd_477HSE_dist_ref
        diff_dist2 = m1_kd_516PHE_m2_kd_509ARG_dist - m1_kd_516PHE_m2_kd_509ARG_dist_ref
        diff_dist3 = m1_kd_530GLN_m2_kd_510HSD_dist - m1_kd_530GLN_m2_kd_510HSD_dist_ref


        if diff_dist1 < 1 and diff_dist2 < 1 and diff_dist3 < 1:
            LOGGER.debug("Protein structure {} has a KD m1 domain and a KD m2 domain close to the dimer structure with distance differences dist1:{}, dist2:{}, dist3:{}".format(name,diff_dist1,diff_dist2,diff_dist3))
            # return True
            return 1 
        elif diff_dist1 > 10 and diff_dist2 > 10 and diff_dist3 > 10:
            LOGGER.debug("Protein structure {} does not have a KD m1 domain and a KD m2 domain close to the dimer structure with distance differences dist1:{}, dist2:{}, dist3:{}".format(name,diff_dist1,diff_dist2,diff_dist3))
            # return False
            return 0
        elif (diff_dist1 > 1 and diff_dist1 < 10) and (diff_dist2 > 1 and diff_dist2 < 10) and (diff_dist3 > 1 and diff_dist3 < 10):  
            LOGGER.debug("Protein structure {} has a KD m1 domain and a KD m2 domain slightly away from the dimer structure with distance differences dist1:{}, dist2:{}, dist3:{}".format(name,diff_dist1,diff_dist2,diff_dist3))
            # return 'TF'
            return -1
        elif diff_dist1 < 10 and diff_dist2 < 10 and diff_dist3 < 10:
            LOGGER.debug("Protein structure {} has a KD domain slightly away from and a KD domain close to the dimer structure with distance differences dist1:{}, dist2:{}, dist3:{}".format(name,diff_dist1,diff_dist2,diff_dist3))
            # return 'TF'
            return -1
        elif diff_dist1 > 10 or diff_dist2 > 10 or diff_dist3 > 10:
            LOGGER.debug("Protein structure {} does not have a KD m1 domain and a KD m2 domain close to the dimer structure with distance differences dist1:{}, dist2:{}, dist3:{}".format(name,diff_dist1,diff_dist2,diff_dist3))
            # return False
            return 0

    @staticmethod
    def m2_1433(u: mda.Universe, ref: mda.Universe, name: str) -> int:

        m1_1433_127ARG = u.select_atoms('name BB and index 2114')
        
        m2_kd_729SEP = u.select_atoms('name BB and index 4679')

        m1_1433_127ARG_com = m1_1433_127ARG.center_of_mass(compound='group')
        m2_kd_729SEP_com = m2_kd_729SEP.center_of_mass(compound='group')

        m1_1433_127ARG_m2_kd_729SEP_dist = distances.distance_array(m1_1433_127ARG_com, m2_kd_729SEP_com,
                                                    box=u.dimensions)
        

        m1_1433_127ARG_ref = ref.select_atoms('name BB and index 1786')

        m2_kd_729SEP_ref = ref.select_atoms('name BB and index 665')

        m1_1433_127ARG_com_ref = m1_1433_127ARG_ref.center_of_mass(compound='group')
        m2_kd_729SEP_com_ref = m2_kd_729SEP_ref.center_of_mass(compound='group')

        m1_1433_127ARG_m2_kd_729SEP_dist_ref = distances.distance_array(m1_1433_127ARG_com_ref, m2_kd_729SEP_com_ref,
                                                    box=u.dimensions)
        
        diff_dist = m1_1433_127ARG_m2_kd_729SEP_dist - m1_1433_127ARG_m2_kd_729SEP_dist_ref

        
        if diff_dist < 0:
            LOGGER.debug("Protein structure {} has a 14-3-3-365 and a KD m2 domain close to the dimer structure with a distance difference: {}".format(name,diff_dist))
            # return True
            return 1
        elif diff_dist > 8:
            LOGGER.debug("Protein structure {} doesn't have a 14-3-3-365 and a KD m2 domain close to the dimer structure with a distance difference: {}".format(name,diff_dist))
            # return False
            return 0
        elif diff_dist > 0 and diff_dist < 8:
            LOGGER.debug("Protein structure {} has a 14-3-3-365 and a KD m2 domain slightly away from the dimer structure with a distance difference: {}".format(name,diff_dist))
            # return 'TF'
            return -1
        
    @staticmethod 
    def generate_sequences(sequence: List[Any], index: int = 0, current_sequence: str = "", results: List[Any] =None) -> List:
        """
        Generate all possible sequences based on the rules:
        - 0 → 'N'
        - 1 → 'Y'
        - -1 → Both 'F' and 'Y'
        
        :param sequence: List of integers (e.g., [0, 1, -1, 1])
        :param index: Current index in the sequence being processed.
        :param current_sequence: The sequence being built so far.
        :param results: A list to store all final sequences.
        :return: A list of all possible sequences.
        """
        if results is None:
            results = []

        # Base case: If we've processed the entire sequence, store the result
        if index == len(sequence):
            results.append(current_sequence)
            return results

        # Recursive case: Process the current digit
        if sequence[index] == 0:
            # 0 → 'F'
            CGValidator_C3.generate_sequences(sequence, index + 1, current_sequence + "N", results)
        elif sequence[index] == 1:
            # 1 → 'Y'
            CGValidator_C3.generate_sequences(sequence, index + 1, current_sequence + "Y", results)
        elif sequence[index] == -1:
            # -1 → Both 'F' and 'Y'
            CGValidator_C3.generate_sequences(sequence, index + 1, current_sequence + "N", results)
            CGValidator_C3.generate_sequences(sequence, index + 1, current_sequence + "Y", results)

        return results

    # @staticmethod
    def _internal_validate(self, data: List[Any]) -> List:
        name, positions, cwd, resource_dir, complex_name, max_mf, max_pe, healing, cleanup, mini_mummi = data
        gromacs_file = None
        valid_status = False
        timeout_gmx = 120 # After 2 mins we kill GROMACS
        start = timer()

        # Make temp dir to run minimization in
        fullPathDir = os.path.join(cwd, name)
        gromacs_file = f"{name}.gro"
        old_cwd = os.getcwd()
        if not os.path.isdir(fullPathDir):
            os.makedirs(fullPathDir, exist_ok=True)

        try:
            os.chdir(fullPathDir)
            LOGGER.debug(f"Directory where files will be outputed {fullPathDir}")
        except Exception as e:
            LOGGER.error(f"chdir failed {fullPathDir} ({e.__class__.__name__} {e})")

        # LOGGER.debug(f"Switched dir to fullPathDir {fullPathDir} ({os.getcwd()})")

        # Save .gro file - @CHECK faster if MDAnalysis universe is in class and here only replace and save BUT maybe prolematic for later threaded version ???
        gromacs_ref_file = os.path.join(resource_dir, complex_name)
        gromacs_ref_file_crystal = os.path.join(resource_dir, 'dimer_ref_crystal.gro')
        gromacs_ref_file_monomer = os.path.join(resource_dir, 'monomer_ref.gro')
        pe = np.inf # Potential energy
        mf = np.inf # Maximum force

        try:
            u = mda.Universe(gromacs_ref_file)
            ref = u.copy() # This needs to be set before we update the positions
            u.atoms.positions = positions
            u_crystal = mda.Universe(gromacs_ref_file_crystal)
            u_monomer = mda.Universe(gromacs_ref_file_monomer)

            if not mini_mummi:
                ref_crystal = u_crystal.copy()
                ref_monomer = u_monomer.copy()
                if (CGValidator_C3.check_knot_crd_loop1(u, 8, name) and
                    CGValidator_C3.check_knot_hvr_rbd(u, 8, name) and
                    CGValidator_C3.check_peptide_1433_distance(u, 6, name) and
                    CGValidator_C3.check_rbd_crd_distance(u, 20, name) and
                    CGValidator_C3.Gdom_rmsd(u,ref,4,17,name) and
                    CGValidator_C3.CRD_rmsd(u,ref,5.5,name) and
                    CGValidator_C3.RBD_rmsd(u,ref,2,11,name)):
                    LOGGER.debug(f"Structure {name} passed the CRD-Loop1 knot, HVR-RBD knot, "
                        f"RAF peptide-14-3-3 helix distance, RBD-CRD COM distance checks, "
                        f"G-domain and HVR checks, CRD-mirror checks, and RBD-mirror and mis-threading checks"
                        f"Proceed with writing the gro file")
                else:
                    end = timer() - start
                    LOGGER.info(f"Structure {name} failed at one of the checks in {end:.3f} seconds")
                    self._cleanup_dir(old_cwd, fullPathDir)
                    return [[False], [name], [-1], [-1]]

            u.atoms.write(f"{name}.gro")
        except ValueError as e:
            end = timer() - start
            LOGGER.error(f"Writing atoms {name} pos=[min={positions.min()}, max={positions.max()}]: {e} (in {end:.3f} seconds)")
            self._cleanup_dir(old_cwd, fullPathDir)
            return [[False], [name], [-1], [-1]]

        # Check structure for needed FF changes
        added_mdp_options = ""
        name_suffixes = []
        if mini_mummi:
            added_mdp_options += " -DACYM280"
        else:
            if CGValidator_C3.zn_cutoff(u, 8, name):
                added_mdp_options += " -DACYM280"
            if CGValidator_C3.SEP365_cutoff(u, 8, name):
                added_mdp_options += " -DASEP365"
            if CGValidator_C3.LB1m1_rmsd_dimer(u,ref_crystal,2.9,name): #1.0
                name_suffixes.append(1) #Y
                if CGValidator_C3.LB1m1_rmsd_monomer(u,ref_monomer,3.5,name):
                    name_suffixes[0] = -1 #YN
            elif CGValidator_C3.LB1m1_rmsd_monomer(u,ref_monomer,3.5,name):
                name_suffixes.append(0) #N
            else:
                end = timer() - start
                LOGGER.info(f"Structure {name} failed at one of the FF checks for m1 in {end:.3f} seconds, exiting..")
                self._cleanup_dir(old_cwd, fullPathDir)
                return [[False], [name], [-1], [-1]]
            
            if CGValidator_C3.LB1m2_rmsd_dimer(u,ref_crystal,2.9,name):  #1.0
                name_suffixes.append(1) #Y
                if CGValidator_C3.LB1m2_rmsd_monomer(u,ref_monomer,3.5,name):
                    name_suffixes[1] = -1 #YN
            elif CGValidator_C3.LB1m2_rmsd_monomer(u,ref_monomer,3.5,name):
                name_suffixes.append(0) #N
            else:
                end = timer() - start
                LOGGER.info(f"Structure {name} failed at one of the FF checks for m2 in {end:.3f} seconds, exiting..")
                self._cleanup_dir(old_cwd, fullPathDir)
                return [[False], [name], [-1], [-1]]

            if CGValidator_C3.m1m2E(u,ref_crystal,name) == 1 and CGValidator_C3.LB1m1_rmsd_dimer(u,ref_crystal,2.9,name) and CGValidator_C3.LB1m2_rmsd_dimer(u,ref_crystal,2.9,name):
                name_suffixes.append(1) #Y
            elif CGValidator_C3.m1m2E(u,ref_crystal,name) == 0:
                name_suffixes.append(0) #N
            else: 
                #-1 value = YN
                if CGValidator_C3.LB1m1_rmsd_dimer(u,ref_crystal,2.9,name) and CGValidator_C3.LB1m2_rmsd_dimer(u,ref_crystal,2.9,name):
                    name_suffixes.append(-1) #YN
                else:
                    name_suffixes.append(0) #N

            if CGValidator_C3.m2_1433(u,ref_crystal,name) == 1 and CGValidator_C3.SEP365_cutoff(u, 8, name):
                name_suffixes.append(1) #Y
            elif CGValidator_C3.m2_1433(u,ref_crystal,name) == 0: 
                name_suffixes.append(0) #N
            else:
                #-1 value = YN
                if CGValidator_C3.SEP365_cutoff(u, 8, name):
                    name_suffixes.append(-1) #YN
                else:
                    name_suffixes.append(0) #N

            # Generate all the combinations
            all_combinations = CGValidator_C3.generate_sequences(name_suffixes)

            # Define distances for defining the force field
            m1_kd_565ASP = u.select_atoms('name BB and index 1389')
            m1_kd_508THR = u.select_atoms('name BB and index 1251')
            m1_kd_516PHE = u.select_atoms('name BB and index 1270')
            m1_1433_127ARG = u.select_atoms('name BB and index 2114')

            m2_kd_477HSE = u.select_atoms('name BB and index 4080')
            m2_kd_509ARG = u.select_atoms('name BB and index 4150')
            m2_kd_729SEP = u.select_atoms('name BB and index 4679')

            # Get center of masses for each part
            m1_kd_565ASP_com = m1_kd_565ASP.center_of_mass(compound='group')
            m1_kd_508THR_com = m1_kd_508THR.center_of_mass(compound='group')
            m1_kd_516PHE_com = m1_kd_516PHE.center_of_mass(compound='group')
            m1_1433_127ARG_com = m1_1433_127ARG.center_of_mass(compound='group')

            m2_kd_477HSE_com = m2_kd_477HSE.center_of_mass(compound='group')
            m2_kd_509ARG_com = m2_kd_509ARG.center_of_mass(compound='group')
            m2_kd_729SEP_com = m2_kd_729SEP.center_of_mass(compound='group')

            m1_kd_565ASP_m2_kd_477HSE_dist = distances.distance_array(m1_kd_565ASP_com, m2_kd_477HSE_com,
                                                box=u.dimensions)

            m1_kd_508THR_m2_kd_509ARG_dist = distances.distance_array(m1_kd_508THR_com, m2_kd_509ARG_com,
                                                box=u.dimensions)

            m1_kd_516PHE_m2_kd_509ARG_dist = distances.distance_array(m1_kd_516PHE_com, m2_kd_509ARG_com,
                                                        box=u.dimensions)

            m1_1433_127ARG_m2_kd_729SEP_dist = distances.distance_array(m1_1433_127ARG_com, m2_kd_729SEP_com,
                                                        box=u.dimensions)
            
            if m1_kd_565ASP_m2_kd_477HSE_dist < 10 and m1_kd_508THR_m2_kd_509ARG_dist < 14 and m1_kd_516PHE_m2_kd_509ARG_dist < 9 or m1_1433_127ARG_m2_kd_729SEP_dist < 15:
                LOGGER.debug(f"Structure {name} used ensemble C forcefield")
                force_field = "_eC"
            else:
                LOGGER.debug(f"Structure {name} used ensemble B* forcefield")
                force_field = "_eBstar"
        name_base = name

        results = {} # hashmap of results : key structure_name

        for combination in all_combinations:
            # print(combination)
            # Not accepted states
            if combination == "YNYY" or combination == "NYYY" or combination == "NNYY" or combination == "YNYN" or combination == "NYYN" or combination == "NNYN":
                continue

            name = name_base + "_" + combination

            # We downsample these two states, we only consider 25% of them
            if combination == "YYYY" or combination == "NNNN":
                if np.random.random() > 0.25:
                    LOGGER.debug(f"We do not consider {name} (downsampling)")
                    continue

            u.atoms.write(f"{name}.gro")

            results[name] = {"valid_status": False, "pe": -1, "mf": -1}

            if combination[0] == 'Y':
                added_mdp_options += " -DM1DIM"
            if combination[1] == 'Y':
                added_mdp_options += " -DM2DIM"
            if combination[2] == 'Y':
                added_mdp_options += " -DM1M2E"
            if combination[3] == 'Y':
                added_mdp_options += " -DM21433"

            subp.call(f"cp {resource_dir}/martini_v2.x_new-rf-em-sc-test.mdp .", shell=True)
            subp.call(f"cp {resource_dir}/martini_v2.x_new-rf-em-test.mdp .", shell=True)
            subp.call(f"sed -ie 's/XXX/{added_mdp_options}/g' martini_v2.x_new-rf-em-test.mdp", shell=True)
            subp.call(f"sed -ie 's/XXX/{added_mdp_options}/g' martini_v2.x_new-rf-em-sc-test.mdp", shell=True)

            # Run minimization
            try:
                LOGGER.debug(f"Run GROMACS minimization in dir {fullPathDir}/{name}")

                try:
                    # Run energy minimization (steepest descent, with protein constreins)
                    grompp = subp.run(
                        f'gmx grompp -c {name}.gro -f martini_v2.x_new-rf-em-sc-test.mdp -p {resource_dir}/system{force_field}-sc.top -o topol-sc_{combination}.tpr -maxwarn 5',
                        shell=True,
                        capture_output=True,
                        timeout=timeout_gmx
                    )
                    # @NOTE: some gmx versions are compiled without thread-MPI and do not support setting the number of threads "-nt" so fix OMP_NUM....
                    mdrun = subp.run(
                        f'export OMP_NUM_THREADS=1; gmx mdrun -nt 1 -ntmpi 1 -v -deffnm topol-sc_{combination} -c {name}-em-sc.gro  > md-sc.out 2>&1',
                        shell=True,
                        capture_output=True,
                        timeout=timeout_gmx
                    )
                    grompp2 = subp.run(
                        f'gmx grompp -c {name}-em-sc.gro -f martini_v2.x_new-rf-em-test.mdp -p {resource_dir}/system{force_field}.top -o topol_{combination}.tpr -maxwarn 5',
                        shell=True,
                        capture_output=True,
                        timeout=timeout_gmx
                    )
                    mdrun2 = subp.run(
                        f'export OMP_NUM_THREADS=1; gmx mdrun -nt 1 -ntmpi 1 -v -deffnm topol_{combination} -c {name}-em.gro  > md.out 2>&1',
                        shell=True,
                        capture_output=True,
                        timeout=timeout_gmx
                    )
                except subp.TimeoutExpired as e:
                    end = timer() - start
                    LOGGER.error(f"GROMACS timeout for {name}: {e} (in {end:.3f} seconds)")
                    continue

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
                            if healing:
                                ## This is a healing section, that's the new addition compared to classic CGValidator_C3 class
                                try:
                                    # gmx_output = subp.run(f'echo 0 | gmx trjconv -f {name}-em.gro -s topol.tpr -o {name}-em_whole.gro -pbc whole', shell=True, capture_output=True, timeout=timeout_gmx)
                                    # FA: Changed from whole to cluster to fix the pbc issues in structures for C3
                                    gmx_output = subp.run(
                                        f'echo 0 0 | gmx trjconv -f {name}-em.gro -s topol_{combination}.tpr -o {name}-em_whole.gro -pbc cluster',
                                        shell = True,
                                        capture_output=True,
                                        timeout=timeout_gmx
                                    )
                                except subp.TimeoutExpired as e:
                                    end = timer() - start
                                    LOGGER.error(f"GROMACS timeout for {name}: {e} (in {end:.3f} seconds)")
                                    results[name]["pe"] = pe
                                    results[name]["mf"] = mf
                                    continue

                                whole_gro_file = f"{name}-em_whole.gro"
                                try:
                                    u = mda.Universe(whole_gro_file)
                                    positions_new = u.atoms.positions
                                    positions_new[:,2] = positions_new[:,2]  - positions_new[430,2]
                                    positions_new[:,1] = positions_new[:,1]  - positions_new[430,1]
                                    if positions_new[430,0] < 10:
                                        positions_new[:,0] += 500 # Box size
                                    u.atoms.positions = positions_new

                                    cyf_z_distance = abs(positions_new[430,2] - positions_new[3327,2])  # Difference in z-coordinates
                                    LOGGER.info(f"Current z-axis distance: {cyf_z_distance} A")
                                    if cyf_z_distance >= 15:
                                        LOGGER.info(f"Need fix")


                                        # FA: Protocol for tilting protein to bring two CYF groups to a specific distance
                                        target_dist = 10.0 # Angstrom

                                        ts = u.trajectory.ts
                                        farnesyl_1 = u.select_atoms("resname CYF and index 430")
                                        farnesyl_2 = u.select_atoms("resname CYF and index 3327")

                                        positions = u.atoms.positions
                                        pos_far1 = positions[430,2]
                                        pos_far2 = positions[3327,2]

                                        # determine anchor point based on the location of CYF
                                        if pos_far2 > pos_far1:
                                            center = farnesyl_2.center_of_mass()
                                        else:
                                            center = farnesyl_1.center_of_mass()

                                        # rotate around the x axes
                                        axis_direction = np.array([1.0, 0.0, 0.0])

                                        cyf_z_distance = abs(positions[430,2] - positions[3327,2])
                                        cyf_com_yz_dist = math.sqrt((positions[430,1] - positions[3327,1])**2 + (positions[430,2] - positions[3327,2])**2)
                                        # FA: We can use the sign of this distance to determine the sign of rotation angle
                                        cyf_y_distance_init = positions[430,1] - positions[3327,1]
                                        cyf_y_distance_final = math.sqrt((cyf_com_yz_dist)**2 - (target_dist)**2)
                                        if cyf_y_distance_init < 0:
                                            cyf_y_distance_final = -1 * cyf_y_distance_final

                                        # get the rotation angle
                                        angle_deg_init = math.degrees(math.atan(cyf_z_distance/cyf_y_distance_init))
                                        angle_deg_final = math.degrees(math.atan(target_dist/cyf_y_distance_final))
                                        angle_deg_rot = angle_deg_init - angle_deg_final

                                        # FA: define the rotation and apply it to the structure in the universe
                                        rotation = transformations.rotate.rotateby(angle_deg_rot, direction=axis_direction, point=center)
                                        u.trajectory.ts = rotation(u.trajectory.ts)

                                        positions_new = u.atoms.positions
                                        cyf_z_distance = abs(positions_new[430,2] - positions_new[3327,2])
                                        LOGGER.info(f"New distance between two CYFs after tilting: {cyf_z_distance}")

                                        # bring the first CYF to its initial position
                                        old_coordinates = farnesyl_1[0].position
                                        new_coordinates = np.array([0.0, 0.0, 0.0])
                                        shift = new_coordinates - old_coordinates
                                        all_atoms = u.select_atoms("all")
                                        all_atoms.positions += shift

                                    # Checking the overlap of BB beads
                                    bb_beads = u.select_atoms("name BB")
                                        # Check if the z-coordinate of any BB bead is greater than 20 Å or lower than -120 Å
                                    z_positions = bb_beads.positions[:, 2]  # Extract the z-coordinates
                                    z_threshold_BB_in_membrane = 20
                                    z_threshold_BB_far_from_membrane = -140
                                    z_greater_than_threshold = z_positions > z_threshold_BB_in_membrane  # Boolean array: True if z > threshold
                                    z_lower_than_threshold = z_positions < z_threshold_BB_far_from_membrane

                                    # Print the result
                                    if z_greater_than_threshold.any() or z_lower_than_threshold.any():
                                        LOGGER.info(f"Structure {name} failed because it has {z_greater_than_threshold.sum()} BB beads in membrane (threshold: {z_threshold_BB_in_membrane}).")
                                        LOGGER.info(f"Structure {name} failed because it has {z_greater_than_threshold.sum()} BB beads in membrane and {z_lower_than_threshold.sum()} BB beads very far from the membrane (upper threshold: {z_threshold_BB_in_membrane}, and lower threshold: {z_threshold_BB_far_from_membrane})")

                                        results[name]["pe"] = pe
                                        results[name]["mf"] = mf
                                        continue

                                    u.atoms.write(whole_gro_file)
                                except ValueError as e:
                                    end = timer() - start
                                    LOGGER.error(f"Validating {name} (whole_gro_file) in {end:.3f} seconds pos=[min={positions_new.min()}, max={positions_new.max()}]: {e}")
                                    results[name]["pe"] = pe
                                    results[name]["mf"] = mf
                                    continue

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
                            
                            results[name]["pe"] = pe
                            results[name]["mf"] = mf
                            continue
                    else:
                        end = timer() - start
                        LOGGER.error(f"ERROR in GROMACS minimization of {name} - file md.out not found")
                        # print(grompp2.stderr.decode())
                        # print(grompp2.stdout.decode())
                        results[name]["pe"] = pe
                        results[name]["mf"] = mf
                        continue
                else:
                    LOGGER.error(f"ERROR in GROMACS minimization of {name} - file {name}-em.gro not found")
                    # print(grompp.stderr.decode())
                    # print(grompp.stdout.decode())
                    results[name]["pe"] = pe
                    results[name]["mf"] = mf
                    continue
            except Exception as error:
                LOGGER.error(f"ERROR in GROMACS minimization of {name} - error was: {error.__class__.__name__}: {error}")
                results[name]["pe"] = pe
                results[name]["mf"] = mf
                continue

            # Valid status
            results[name]["valid_status"] = True
            results[name]["pe"] = pe
            results[name]["mf"] = mf

        
        # Cleaning up
        self._cleanup_dir(old_cwd, fullPathDir)

        structure_names = list(results.keys())
        valid_status_comb = [] 
        pe_all = []
        mf_all = []

        for key, val in results.items():
            valid_status_comb.append(val["valid_status"])
            pe_all.append(val["pe"])
            mf_all.append(val["mf"])

        # Return all valid FF
        return [valid_status_comb, structure_names, mf_all, pe_all]
