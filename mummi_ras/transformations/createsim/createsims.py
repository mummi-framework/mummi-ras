#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2024, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
#
# @authors  Helgi I. Ingolfsson <ingolfsson1@llnl.gov>
#           Changhoan Kim <kimchang@us.ibm.com>
#           Francesco Di Natale <dinatale3@llnl.gov>
#           Loic Pottier <pottier1@llnl.gov>
# -----------------------------------------------------------------------------

"""Functions to setup, run and map between macro and CG resolutions."""

from argparse import ArgumentParser, RawTextHelpFormatter
import copy
import inspect
import logging
import os
import math
import sys
import pickle
import numpy as np
from random import randint
import fnmatch
import re
import shutil
import signal
import time
import traceback
from logging import getLogger
import MDAnalysis as mda
from MDAnalysis.analysis import distances

import mummi_core
import mummi_ras
from mummi_core.utils import Naming
from mummi_core.utils import utilities as utils
import mummi_ras.transformations.createsim.particlesim as particlesim
import mummi_ras.transformations.createsim.placeproteins as placeproteins
import mummi_ras.transformations.createsim.pdcglob as pdcglob
from mummi_ras.transformations.patch_creator import fetch_gc_patches
from mummi_ras.datastructures.patch import Patch as pfPatch
from mummi_ras.online.cg import get_tiltrot_z_state_multiple
from mummi_core.utils.logger import init_logger

# Program Globals
LOGGER = getLogger(__name__)
global gromacs, mdrun_mpi, mdrun_o


# Two functions from Tim C. to check structure and pick what FF to use 
# Note cutoff values should be in Angstroms 
# @TODO this analysis will also go into cgAnalysis and this test will also be used in createsims consolidate code
def zn_cutoff(u, cutoff, name):
    ZN740 = u.select_atoms('name ZN and index 1805')#and resid 740
    CYS280 = u.select_atoms('name SC1 and index 750')#and resid 280
    ZN_distances = mda.analysis.distances.distance_array(ZN740, CYS280)
    if ZN_distances < cutoff:
        LOGGER.info("Protein structure {} has ZN740 with {} of CYS280 less than cutoff {} use -DACYM280 parameters".format(name, ZN_distances, cutoff))
        return True
    else:
        LOGGER.info("Protein structure {} has ZN740 with {} of CYS280 more than cutoff {} dont use -DACYM280 parameters".format(name, ZN_distances, cutoff))
        return False

def SEP365_cutoff(u, cutoff, name):
    SEP365 = u.select_atoms('resname SEP and index 942 and name SC1')
    SEP365_ARG56 = u.select_atoms('resname ARG and index 1943 and name SC2')# and resid 56
    SEP365_ARG127 = u.select_atoms('resname ARG and index 2116 and name SC2')#and resid 127
    SEP_distances1 = mda.analysis.distances.distance_array(SEP365, SEP365_ARG56)
    SEP_distances2 = mda.analysis.distances.distance_array(SEP365, SEP365_ARG127)
    #print(SEP_distances1, SEP_distances2)
    if SEP_distances1 > cutoff or SEP_distances2 > cutoff:
        LOGGER.info("Protein structure {} has SEP365 {} and {} from ARG56 and 127, more than cutoff {} use -DASEP365 parameters".format(name, SEP_distances1, SEP_distances2, cutoff))
        return True
    else:
        LOGGER.info("Protein structure {} has SEP365 {} and {} from ARG56 and 127, within cutoff {} dont use -DASEP365 parameters".format(name, SEP_distances1, SEP_distances2, cutoff))
        return False


def setup_cg_sim(pSim, simpath, mini_mummi = False):
    """Set up and run eq for the CG sim."""
    LOGGER.info("---> Setup CG sim and run eq in dir {}".format(simpath))
    os.chdir(simpath)

    # Get all .mdp file locally and modify flags if needed 
    utils.sys_call("cp "+particlesim.dirMartini+"/martini_v2.x_new-rf*.mdp .") # copy all - should be 9 files 
    added_mdp_options = ""
    if (len(pSim.proteinUnits) > 0):
        added_mdp_options = pSim.proteinUnits[0]["ff_options"] # only one set of options (currently only tested with one protein)
    utils.sys_call(f"sed -ie 's/XXX/{added_mdp_options}/g' martini_v2.x_new-rf-*.mdp") # will change all x9

    # Get protein structures insertion string and add random rotation / placement as needed
    protein_string = ""
    # @TODO put this back in with alignment but also support no allignment if protein sturcture is specified (use below)

    if (len(pSim.proteinUnits) > 0):

        if mini_mummi:
            proteinRotated = [pSim.proteinUnits[0]["structure"]] # for mini-MuMMI don't change inital structure
        else:
            proteinRotated = ["struct_createsim.gro"]
        proteinPosInserted = [pSim.proteinUnits[0]["position"]]

        LOGGER.info("---> Proteins are pressent they are: {}".format(pSim.proteinUnits))
        LOGGER.info("---> After protein alignment: rotation {} posInserted {}".format(proteinRotated, proteinPosInserted))
        for i in range(len(pSim.proteinUnits)):
            cyfGroupX = -1
            cyfGroupY = -1
            cypGroupX = -1
            cypGroupY = -1
            if pSim.proteinUnits[i]["type"] == "RAS4A-RAF" or pSim.proteinUnits[i]["type"] == "RAS4A-ONLY":
                # Have kras4A need to fined location of CYF and CYP to feed to insane
                protein = mda.Universe(proteinRotated[i])
                farCOM = protein.select_atoms("resname CYF and name SC1").center_of_mass() / 10 # Rember to convert from A to nm
                palCOM = protein.select_atoms("resname CYP and name SC1").center_of_mass() / 10 # Rember to convert from A to nm
                LOGGER.info("---> RAS4A CYF = {},{} vs CYP = {},{} - COM = {},{}".format(farCOM[0], farCOM[1], palCOM[0], palCOM[1], proteinPosInserted[i][0], proteinPosInserted[i][1]))
                cyfGroupX = farCOM[0]
                cyfGroupY = farCOM[1]
                cypGroupX = palCOM[0]
                cypGroupY = palCOM[1]
            if pSim.proteinUnits[i]["type"] == "dimercomp":
                # Note a modified insane is used that knows the x2 locations are both FAR for the same complex - insane-patch-build-m3-c3.py
                # For Adm-C3 m1 FAR 
                cyfGroupX = proteinPosInserted[i][0]
                cyfGroupY = proteinPosInserted[i][1]
                # For Adm-C3 m2 FAR 
                cypGroupX = pSim.proteinUnits[0]["position_m2"][0]
                cypGroupY = pSim.proteinUnits[0]["position_m2"][1]
            else:
                # Have kras4B need to fined location of CYF only leave cyp as -1
                cyfGroupX = proteinPosInserted[i][0]
                cyfGroupY = proteinPosInserted[i][1]
                cypGroupX = -1
                cypGroupY = -1

                # Test only - remove after checking as MDAnalysis takes time to load
                if False:
                    protein = mda.Universe(proteinRotated[i])
                    farCOM = protein.select_atoms("resname CYF and name SC1").center_of_mass() / 10 # Rember to convert from A to nm
                    LOGGER.info("---> RAS4B CYF = {},{} - COM = {},{}".format(farCOM[0], farCOM[1], proteinPosInserted[i][0], proteinPosInserted[i][1]))

            protein_string += " -ras " + proteinRotated[i] + " l " + str(cyfGroupX) + " " + str(cyfGroupY) + " " + str(cypGroupX) + " " + str(cypGroupY) + " "

    LOGGER.info("---> Manual protein alignment insane protein string = {}".format(protein_string))

    # If inital setup / em fails tray again untill max tries
    maxTries = 10
    for i in range(maxTries):
        LOGGER.info("---> System setup/EM try number {}".format(i))

        # Make inital cordinates (file overwritten by insane so need to generate each time)
        utils.get_file(particlesim.dirMartini+"/"+pSim.defTopFile, 'system.top')
        utils.sys_call("sed -ie 's/xPATHx/"+particlesim.dirMartini.replace("/", "\/")+"/g' system.top")

        # WARNING if we have many CPUs compaired to the simulation size we will get into domain problems - expecially in the begining
        # so change -nt depending on system size e.g. for 15x15 nm system 12 CPUs is good initally
        ntString = ""  # default is empty ie let GROMACS use max CPUs
        if pSim.cgnt > 0:
            ntString = " -nt "+str(pSim.cgnt)+" "

        # Get total protein charge - this is now stored for each protein unit
        # Note: use -2 for RAS and -1 for RAS_RAF WARNING RAS-RAF is missing a carge of +4 this not included here so ions need to be addjusted in backmapping
        pCharge = int(np.sum([x['charge'] for x in pSim.proteinUnits]))

        # Setup simulation with modified insane
        # For the real system we want something like this size 31x31 nm or 3200 lipids (1600,1600 - 40x40 per leaflet
        # this will support x25 8x8 lipid subpatches).
        # For subpatches the lipid ratios have to fit a 64 lipid patch with no rounding error
        # safest to make consentration of all lipids add up to 64.
        # TEST - e.g. 16x16 nm give 800 (400,400 - 20x20 per leaflet) lipids - then we can try x4 (10x10 subpatches
        # with lipid cons. adding to 100)
        # WARNING with this version of KRAS the protein charge is +0 but it's detected by insane as +2 (so addjust acordingly)
        # WARNING as it's now ofset to lipid charge (auto calculated in insane) e.g. currently add x1 per KRAS molecule

        if pSim.lipidFullAsignment:
            # Assign lipids based on full assignment of each lipid as specified in -lgrid lipid_counts.npz
            utils.sys_call('python '+particlesim.dirMartini+'/insane-patch-build-m3-c3.py' +
                           ' -salt ' + str(pSim.ioncons) +
                           ' -x '+str(pSim.sizeX)+' -y '+str(pSim.sizeY)+' -z ' + str(pSim.sizeZ) +
                           ' -charge '+str(pCharge) +
                           ' -pbc square' +
                           ' -sol ' + pSim.water +
                           ' -asym ' + str(pSim.asym) +
                           ' -rand 0.095 -o lipids-water.gro' +
                           ' -lgrid lipid_counts.npz ' +
                           ' '+protein_string+' 2>&1 | tee -a system.top')
        else:
            # Assign lipids based on probabilities as specified in -uf and -lf files (always use Martini 3, Martini 2 is legacy)
            utils.sys_call('python '+particlesim.dirMartini+'/ insane-ck-m3.py' +
                           ' -salt ' + str(pSim.ioncons) +
                           ' -x '+str(pSim.sizeX)+' -y '+str(pSim.sizeY)+' -z ' + str(pSim.sizeZ) +
                           ' -charge '+str(pCharge) +
                           ' -pbc square' +
                           ' -sol ' + pSim.water +
                           ' -asym ' + str(pSim.asym) +
                           ' -rand 0.095 -o lipids-water.gro' +
                           ' -uf ' + particlesim.pSimSubPatchFileUpper +
                           ' -lf ' + particlesim.pSimSubPatchFileLower +
                           ' '+protein_string+' 2>&1 | tee -a system.top')

        # Turn on for testing lipid counts
        #save_insane_lipid_counts()
        #exit()

        # Run EM
        # If there are KRAS molecules: save .top for later, now construct one where each KRAS moleucle is split into x2 parts
        #   KRAS and a x3 bead FAR
        if (len(pSim.proteinUnits) > 0):
            # Prep pulling protein .top file
            #utils.sys_call("cp system.top system-pull.top") # not needed in M3 PULL is controled in .mdp options 

            # Update insane generated .top file to ref correct proteins
            #utils.sys_call("sed -ie 's/Protein          1/KRAS_pull        1/g' system-pull.top")
            #utils.sys_call("sed -ie 's/Protein          1/KRAS_GTP         1/g' system.top")
            #utils.sys_call("sed -ie 's/Protein          1/RASR_PULL        1/g' system-pull.top")
            #utils.sys_call("sed -ie 's/Protein          1/RAS_RAF         1/g' system.top")
            for j in range(len(pSim.proteinUnits)):
                if pSim.proteinUnits[j]["type"] == "RAS-ONLY":
                    utils.sys_call("sed -ie '0,/Protein          1/s//RAS_PULL         1/g' system-pull.top")
                    utils.sys_call("sed -ie '0,/Protein          1/s//RAS              1/g' system.top")
                elif pSim.proteinUnits[j]["type"] == "RAS-RAF":
                    # In Adm mini-MuMMI the pulling has been added as options into the default .itp
                    if mini_mummi:
                        utils.sys_call("sed -ie '0,/Protein          1/s//RAS_RAF          1/g' system.top")
                    else:
                        utils.sys_call("sed -ie '0,/Protein          1/s//RASR_PULL        1/g' system-pull.top")
                elif pSim.proteinUnits[j]["type"] == "RAS4A-ONLY":
                    utils.sys_call("sed -ie '0,/Protein          1/s//RAS4A_PULL       1/g' system-pull.top")
                    utils.sys_call("sed -ie '0,/Protein          1/s//RAS4A            1/g' system.top")
                elif pSim.proteinUnits[j]["type"] == "RAS4A-RAF":
                    utils.sys_call("sed -ie '0,/Protein          1/s//RAS4A_RAF_PULL   1/g' system-pull.top")
                    utils.sys_call("sed -ie '0,/Protein          1/s//RAS4A_RAF        1/g' system.top")
                elif pSim.proteinUnits[j]["type"] == "RAS_RAF_14_3_3":
                    utils.sys_call("sed -ie '0,/Protein          1/s//KRBR14           1/g' system.top")
                elif pSim.proteinUnits[j]["type"] == "dimercomp":
                    utils.sys_call("sed -ie '0,/Protein          1/s//c3complex        1/g' system.top")
                else:
                    error = "Protein type {} not supported".format(pSim.proteinUnits[j]["type"])
                    LOGGER.error(error)
                    raise AttributeError(error)

           # Run energy minimization (steepest descent, with protein constreins)
            utils.sys_call(gromacs+' grompp -c lipids-water.gro -r lipids-water.gro -f martini_v2.x_new-rf-em.mdp -p system.top -o topol.tpr -maxwarn 5')
            # @TODO fix this hack? ntString can't be used as we don't want to limit the CPU's of the proper md runs below
            #utils.sys_call(gromacs+' mdrun '+ntString+' -v -c lipids-water-em.gro >> mdrun-em.log 2>&1')
            # @TODO need to change this
            #utils.sys_call(gromacs+' mdrun -ntmpi 20 -v -c lipids-water-em.gro >> mdrun-em.log 2>&1')
            utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-em.gro >> mdrun-em.log 2>&1')
            utils.sys_call('echo 0 | ' + gromacs+' trjconv -f lipids-water-em.gro -o lipids-water-em.gro -pbc mol -s topol.tpr')
        else:  # No KRAS bilayer only simulations
            # Run energy minimization (steepest descent, with protein constreins)
            utils.sys_call(gromacs+' grompp -c lipids-water.gro -r lipids-water.gro -f martini_v2.x_new-rf-em.mdp -p system.top -o topol.tpr -maxwarn 5')
            utils.sys_call(gromacs+' mdrun -ntmpi 20 -v -c lipids-water-em.gro >> mdrun-em.log 2>&1')

        # If em sucesfull don't tray again
        if os.path.isfile("lipids-water-em.gro"):
            # GROMACS 2019.06 can write .gro files even with inf energy check if inf only if not break
            file = open("mdrun-em.log", "r")
            foundInf=False
            for line in file:
                if re.search("Norm of force", line):
                    engValue = line.split("=")[1].strip()
                    if engValue != "inf":
                        #print(line.split("=")[1].strip())
                        foundInf=False
                        break
                    else:
                        LOGGER.info("---> Energy minimization failed in try {}, inf force found, {}, retry".format(i, line))
                        foundInf=True
                        break
            if not foundInf:  # No error exit insne forloop
                break
        else:
            LOGGER.info("---> Energy minimization failed in try {}, no lipids-water-em.gro file, retry".format(i))

    # EM sucesfull continue with minimization
    # If there are KRAS molecules: save .top for later, now construct one where each KRAS moleucle is split into x2 parts
    #   KRAS and a x3 bead FAR
    if (len(pSim.proteinUnits) > 0):
        any4A = False # Flag to know if any of the proteins are 4A type
        if any("4A" in proteinUnit['type'] for proteinUnit in pSim.proteinUnits):
            any4A = True

        # Make index files for groups Lipids and Rest - Note the name of the water changes W, PW, NW
        tempOutput = open("index-selection.txt", 'w')
        tempOutput.write("del 1 - 200\n")
        tempOutput.write("r "+pSim.water+" | r NA | r CL\n")
        tempOutput.write("name 1 Solvent\n")
        tempOutput.write("!1\n")
        tempOutput.write("name 2 Rest\n")
        #tempOutput.write("r CYF & a SC1\n") # changed names in M3
        #tempOutput.write("r CYF & a BB\n")        
        tempOutput.write("a BB\n") # Adm-C2 change - pull protein COM using BB only 
        tempOutput.write("name 3 PULL1\n")
        #tempOutput.write("r CYF & a F1\n") # changed names in M3
        tempOutput.write("r CYF & a C1\n")       
        tempOutput.write("name 4 PULL2\n")
        # if any of the proteins is a RAS4A or RAS4A-RAF we allso need to pull CYP
        if any4A:
            tempOutput.write("r CYP & a SC1\n")
            tempOutput.write("name 5 PULL3\n")
            tempOutput.write("r CYP & a P1\n")
            tempOutput.write("name 6 PULL4\n")
        tempOutput.write("q\n")
        tempOutput.close()
        utils.sys_call(gromacs+' make_ndx -f lipids-water-em.gro -o index.ndx < index-selection.txt')

        # Count number of pulls for this sim one per RAS and RAS-RAF and two per RAS4A and RAS4A-RAF
        pullCount = 0
        for i in range(len(pSim.proteinUnits)):
            if "4A" in pSim.proteinUnits[i]["type"]:
                pullCount = pullCount + 2
            else:
                pullCount = pullCount + 1

        # @TODO here check all pull groups - get min dist - optimize pulling for min dist - then repeat in a next step

        # We need to split the PULL group's if more than one (also edit the .mdp file)
        if (len(pSim.proteinUnits) > 1):
            lineCount = 4
            if any4A:
                lineCount = 8

            utils.sys_call("mv index.ndx index-old.ndx")
            utils.sys_call("tail -n {} index-old.ndx > temp-index-end.txt".format(lineCount))

            # 2018.01.23 - Fixed to work on all BSD systems
            # WARNING head -n -X won't work on MAC - change this to something else
            # utils.sys_call("head -n -4 index-old.ndx > index.ndx")
            utils.sys_call("head -n $(( $(wc -l index-old.ndx | awk '{print $1}') - "+str(lineCount)+" )) index-old.ndx > index.ndx")

            tempIndexFile = open("temp-index-end.txt", "r")
            tempLines = tempIndexFile.readlines()
            p1 = tempLines[1].strip().split()
            p2 = tempLines[3].strip().split()
            if (len(p1) != len(p2) or len(p1) != len(pSim.proteinUnits)):
                error = "---> CG sim " + pSim.name + '/' + simpath + " should be " + len(pSim.proteinUnits) + \
                        " RAS but only found " + len(p1) + " in index file"
                LOGGER.error(error)
                raise ValueError(error)
            else:
                LOGGER.info("Found {} RAS pull groups for fransyl".format(len(p1)))
            if any4A:
                p3 = tempLines[5].strip().split()
                p4 = tempLines[7].strip().split()
                if (len(p3) != len(p4)):
                    error = "---> CG sim " + pSim.name + '/' + simpath + " somethign wrong with 4A pull groups"
                    LOGGER.error(error)
                    raise ValueError(error)
                else:
                    LOGGER.info("Found {} RAS 4A pull groups for palmitoyl".format(len(p3)))
            tempIndexFile.close()

            tempOutput = open("index.ndx", 'a')
            for i in range(len(p1)):
                tempOutput.write("[ PULL"+str(i*2+1)+" ] \n")
                tempOutput.write(str(p1[i])+" \n")
                tempOutput.write("[ PULL"+str(i*2+2)+" ] \n")
                tempOutput.write(str(p2[i])+" \n")
            if any4A:
                for i in range(len(p3)):
                    tempOutput.write("[ PULL"+str((i+len(p1))*2+1)+" ] \n")
                    tempOutput.write(str(p3[i])+" \n")
                    tempOutput.write("[ PULL"+str((i+len(p1))*2+2)+" ] \n")
                    tempOutput.write(str(p4[i])+" \n")
            tempOutput.close()

        # Copy in PULL .mdp file - we have to add pull groups
        pullMDPfile = "martini_v2.x_new-rf-eq3-PULL.mdp"
        #utils.sys_call("cp "+particlesim.dirMartini+"/"+pullMDPfile+" .") # now copy all .mdp files above


        tempOutput = open(pullMDPfile, 'a')
        tempOutput.write("; COM PULLING \n")
        tempOutput.write("pull                     = yes \n")
        tempOutput.write("pull_nstxout             = 100 \n")
        tempOutput.write("pull_nstfout             = 100 \n")
        tempOutput.write("pull_ngroups             = "+str(pullCount*2)+"    ; Number of pull groups \n")
        tempOutput.write("pull_ncoords             = "+str(pullCount)+" \n")
        tempOutput.write(" \n")
        # Add pulling on CYF for all RAS
        for i in range(len(pSim.proteinUnits)):
            tempOutput.write("; Groups - KRAS CYF #"+str(i+1)+" \n")
            tempOutput.write("pull_group"+str(i*2+1)+"_name         = PULL"+str(i*2+1)+" \n")
            #tempOutput.write("pull_group"+str(i*2+1)+"_pbcatom      = 0 \n")
            tempOutput.write("pull_group"+str(i*2+1)+"_pbcatom      = 431 \n") # Adm-C2 change @TODO fix should not be hardcoded value 
            tempOutput.write("pull_group"+str(i*2+2)+"_name         = PULL"+str(i*2+2)+" \n")
            tempOutput.write("pull_group"+str(i*2+2)+"_pbcatom      = 432 \n") # Adm-C2 change @TODO fix should not be hardcoded value
            #tempOutput.write("pull_group"+str(i*2+2)+"_pbcatom      = 0 \n") 
            tempOutput.write("pull_coord"+str(i+1)+"_type         = flat-bottom \n")
            tempOutput.write("pull_coord"+str(i+1)+"init          = 0.1 \n")  # was 0.4
            tempOutput.write("pull_coord"+str(i+1)+"_geometry     = distance  \n")
            tempOutput.write("pull_coord"+str(i+1)+"_groups       = "+str(i*2+1)+" "+str(i*2+2)+" \n")
            tempOutput.write("pull_coord"+str(i+1)+"_dim          = N N Y  \n")
            tempOutput.write("pull_coord"+str(i+1)+"_rate         = -0.0003  \n") # [nm/ps] - this is 1.2 nm on current length 400 000 * 0.01 = 4000 ps
            # Make less to pull slower?
            # 2021.07.02 400000 * 0.01 = 4000 ps, -0.00045 * 4000 = 1.8 nm. They should all start with a 2.1 nm diff
            # Could change and over pull using e.g. -0.0005 - might help when there is a CYF and CYP mismatch
            # in current martini_v2.x_new-rf-eq3-PW-PULL.mdp using 2,000,000 * 0.01 ps = 20,000 ps
            #   or -0.0003 nm/ps * 20,000 ps = -6 nm
            # This is a little bit more than is needed in the current test simulations
            # @TODO for production reduce box size have these somewhat closer and pull slower
            tempOutput.write("pull_coord"+str(i+1)+"_k            = 1000 \n")
            tempOutput.write("pull_coord"+str(i+1)+"_start        = yes  \n")
            tempOutput.write("pull-pbc-ref-prev-step-com          = yes \n") # Adm-C2 change (for the pdbatoms)
            tempOutput.write(" \n")

        # Add pulling on CYP for all 4A variants
        count4A=0
        for j in range(len(pSim.proteinUnits)):
            if "4A" in pSim.proteinUnits[j]["type"]:
                i = count4A + len(pSim.proteinUnits) # i is used as the offset after all CYF pulling
                count4A += 1
                tempOutput.write("; Groups - KRAS 4A CYP #"+str(i+1)+" \n")
                tempOutput.write("pull_group"+str(i*2+1)+"_name         = PULL"+str(i*2+1)+" \n")
                tempOutput.write("pull_group"+str(i*2+1)+"_pbcatom      = 0 \n")
                tempOutput.write("pull_group"+str(i*2+2)+"_name         = PULL"+str(i*2+2)+" \n")
                tempOutput.write("pull_group"+str(i*2+2)+"_pbcatom      = 0 \n")
                tempOutput.write("pull_coord"+str(i+1)+"_type         = flat-bottom \n")  # flat-bottom-high
                tempOutput.write("pull_coord"+str(i+1)+"init          = 0.0 \n") # 0.0 - this is added to the inital distance (with _start = yes) starting further appart
                tempOutput.write("pull_coord"+str(i+1)+"_geometry     = distance  \n")
                tempOutput.write("pull_coord"+str(i+1)+"_groups       = "+str(i*2+1)+" "+str(i*2+2)+" \n")
                tempOutput.write("pull_coord"+str(i+1)+"_dim          = N N Y  \n")
                tempOutput.write("pull_coord"+str(i+1)+"_rate         = -0.00045 \n") # [nm/ps], -0.0005 used for CYF
                tempOutput.write("pull_coord"+str(i+1)+"_k            = 1000 \n")
                tempOutput.write("pull_coord"+str(i+1)+"_start        = yes  \n")
                tempOutput.write(" \n")
        tempOutput.close()

        # Run short equilibration at 1 fs time setp + generate velocites
        utils.sys_call(gromacs+' grompp -c lipids-water-em.gro -f martini_v2.x_new-rf-eq1.mdp ' +
                       ' -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-em.gro')
        utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq1.gro >> mdrun-eq1.log 2>&1 ')

        # Run short equilibration at 5 fs time setp
        # WARNING for PW use martini_v2.x_new-rf-eq2-PW.mdp
        utils.sys_call(gromacs+' grompp -c lipids-water-eq1.gro -f martini_v2.x_new-rf-eq2.mdp ' +
                       ' -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-em.gro')
        utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq2.gro >> mdrun-eq2.log 2>&1')

        # Pull molecules to bilayer here, then change .top to while KRAS molecules - using local PULL .mdp
        utils.sys_call(gromacs+' grompp -c lipids-water-eq2.gro -f '+pullMDPfile +
                       ' -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-em.gro')
        utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq3.gro >> mdrun-eq3.log 2>&1')

    else:  # No protein, bilayer only simulations
        # Make index files for groups Lipids and Rest - Note the name of the water changes W, PW, NW
        tempOutput = open("index-selection.txt", 'w')
        tempOutput.write("del 1 - 200\n")
        tempOutput.write("r "+pSim.water+" | r NA | r CL\n")
        tempOutput.write("name 1 Solvent\n")
        tempOutput.write("!1\n")
        tempOutput.write("name 2 Rest\n")
        tempOutput.write("q\n")
        tempOutput.close()
        utils.sys_call(gromacs+' make_ndx -f lipids-water-em.gro -o index.ndx < index-selection.txt')

        # Run short equilibration at 1 fs time setp + generate velocites
        utils.sys_call(gromacs+' grompp -c lipids-water-em.gro -f martini_v2.x_new-rf-eq1.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-em.gro')
        utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq1.gro >> mdrun-eq1.log 2>&1')

        # Run short equilibration at 5 fs time setp
        # WARNING for PW use martini_v2.x_new-rf-eq2-PW.mdp
        utils.sys_call(gromacs+' grompp -c lipids-water-eq1.gro -f martini_v2.x_new-rf-eq2.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-em.gro')
        utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq2.gro >> mdrun-eq2.log 2>&1')

        # Run short equilibration at 10 fs time setp
        # WARNING for PW use martini_v2.x_new-rf-eq3-PW.mdp
        utils.sys_call(gromacs+' grompp -c lipids-water-eq2.gro -f martini_v2.x_new-rf-eq3.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-em.gro')
        utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq3.gro >> mdrun-eq3.log 2>&1')

    ## @TODO remove this in Campain 3 - this is for campain 1*, to fix the HVR allow it to fix the helix - also edit command below
    ## This run shold allow the HVR to move - fixing the few residues that need to form a helix
    #utils.sys_call(gromacs+' grompp -c lipids-water-eq3.gro -f '+particlesim.dirMartini +
    #               '/martini_v2.x_new-rf-eq4hvr.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-eq3.gro')
    ## To troubleshot topolgy add -pp temp-pp.top
    #utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq4hvr.gro >> mdrun-eq4hvr.log 2>&1')

    # Added a x2 eq steps (eq4 and eq5) to slowly relax protein on membrane - Note, HVR postres are also removed from here on
    utils.sys_call(gromacs+' grompp -c lipids-water-eq3.gro -f martini_v2.x_new-rf-eq4.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-eq3.gro')
    utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq4.gro >> mdrun-eq4.log 2>&1')
    utils.sys_call(gromacs+' grompp -c lipids-water-eq4.gro -f martini_v2.x_new-rf-eq5.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-eq4.gro')
    utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq5.gro >> mdrun-eq5.log 2>&1')

    # Run final equilibration at 20 fs time setp (same as pruduction - test stability + allow the box are to equilibrate
    # @TESTING for now do only 2 ns but should move this to 20-100 ns (for the 3000 lipid systesm)
    # WARNING for PW use martini_v2.x_new-rf-eq4-PW.mdp
    #utils.sys_call(gromacs+' grompp -c lipids-water-eq3.gro -f '+particlesim.dirMartini +
    #               '/martini_v2.x_new-rf-eq4.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-eq3.gro')
    #utils.sys_call(gromacs+' grompp -c lipids-water-eq4hvr.gro -f '+particlesim.dirMartini +
    #               '/martini_v2.x_new-rf-eq4.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-eq4hvr.gro')
    utils.sys_call(gromacs+' grompp -c lipids-water-eq5.gro -f martini_v2.x_new-rf-eq6.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-eq5.gro')
    utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq6.gro >> mdrun-eq6.log 2>&1')

    # Cleanup temp files
    utils.sys_call('rm \#* step*.pdb || true')
    utils.sys_call('rm system.tope state.cpt state_prev.cpt traj_comp.xtc traj.trr mdout.mdp md.log ener.edr index-old.ndx temp-index-end.txt || true')

    # Compress all cretsims files into createsims.tar.gz storage that are not needed for next steps and remove
    utils.sys_call('tar -czvf createsims.tar.gz psimSubPatch_lower.pat psimSubPatch_upper.pat protein-*-rotate.gro protein-ALL.gro ' +
                   'index-selection.txt lipids-water-em.gro lipids-water-eq1.gro lipids-water-eq2.gro lipids-water-eq3.gro ' +
                   'lipids-water.gro martini_v2.x_new-rf-eq3-PULL.mdp mdrun-em.log mdrun-eq1.log mdrun-eq2.log mdrun-eq3.log ' +
                   'mdrun-eq4.log pullf.xvg pullx.xvg lipid_counts.npz mdrun-eq5.log lipids-water-eq5.gro ' +
                   'mdrun-eq6.log lipids-water-eq4.gro macro_sim_dir/snapshot.* macro_sim_dir/macro_createsims.out struct_createsim.gro || true')
    utils.sys_call('rm psimSubPatch_lower.pat psimSubPatch_upper.pat protein-*-rotate.gro protein-ALL.gro ' +
                   'index-selection.txt lipids-water-em.gro lipids-water-eq1.gro lipids-water-eq2.gro lipids-water-eq3.gro ' +
                   'lipids-water.gro martini_v2.x_new-rf-e*.mdp martini_v2.x_new-rf-*.mdpe mdrun-em.log mdrun-eq1.log ' +
                   'mdrun-eq2.log mdrun-eq3.log traj_comp.step*.xtc ' +
                   'mdrun-eq4.log pullf.xvg pullx.xvg lipid_counts.npz mdrun-eq5.log lipids-water-eq5.gro ' +
                   'mdrun-eq6.log lipids-water-eq4.gro struct_createsim.gro || true')
    utils.sys_call('rm -r macro_sim_dir || true')


def run_multible_cg_sim(pSim, simpath, simrange, genvel):
    """Run multible production CG simulation."""
    for i in simrange:
        run_cg_sim(pSim, simpath, rName="-r"+str(i), genvel=genvel)


def prime_cg_sim(pSim, simpath, rName="", genvel=False):
    """Run production CG simulation."""
    LOGGER.info("---> Setup CG production sim in dir {}".format(simpath))
    os.chdir(simpath)

    # Copy over last eq frame and topology
    cgLastFrame = "lipids-water-eq6.gro"
    cgTopFile = "system.top"
    cgIndexFile = "index.ndx"
    #cg_eq_sim = pSim.name + '/cg-eq/'
    #utils.get_file(particlesim.get_sim_sub_dir(pSim, "cg-eq")+'/'+cgLastFrame, cgLastFrame)
    #utils.get_file(particlesim.get_sim_sub_dir(pSim, "cg-eq")+'/'+cgTopFile, cgTopFile)
    #utils.get_file(particlesim.get_sim_sub_dir(pSim, "cg-eq")+'/'+cgIndexFile, cgIndexFile)
    #utils.get_file(cg_eq_sim + cgLastFrame, cgLastFrame)
    #utils.get_file(cg_eq_sim + cgTopFile, cgTopFile)
    #utils.get_file(cg_eq_sim + cgIndexFile, cgIndexFile)

    # Make it less likely to have lipids on both side of the pbc
    utils.sys_call(gromacs+' editconf -f '+cgLastFrame+' -o '+cgLastFrame+' -translate 0 0 -0.5 ')

    # Set run for production (20 fs time step) - set new velocites e.g. if making multible repeates
    if genvel:
        # @WARNING -vel.mdp file is missing in resources @TODO add or remove this feature 
        utils.sys_call(gromacs+' grompp -c '+cgLastFrame+' -r '+cgLastFrame+' -f martini_v2.x_new-rf-prod-vel.mdp -p '+cgTopFile+' -o topol.tpr -maxwarn 5 -n '+cgIndexFile)
    else:
        utils.sys_call(gromacs+' grompp -c '+cgLastFrame+' -r '+cgLastFrame+' -f martini_v2.x_new-rf-prod.mdp -p '+cgTopFile+' -o topol.tpr -maxwarn 5 -n '+cgIndexFile)

    # Cleanup temp files
    utils.sys_call('rm \#* || true')
    utils.sys_call('rm fancy-slurm-syrah.she fancy.she lipids-water-eq4-fix.pdbe mdout.mdp system.tope molecule_1.data molecule_2.data molecule_3.data || true')
    # @TODO do further cleanup before main campaign add: tpr.log  others ??


def run_cg_sim_GROMACS(pSim, simpath, inpath, intype):
    """Run GROMACS in a alreaddy setup directory."""
    # @WARNING this on was used for testing and is currently not working @TODO fix
    LOGGER.info("---> Run CG sim in dir {}".format(simpath))
    os.chdir(simpath)

    # To troubleshot topolgy add -pp temp-pp.top
    utils.sys_call(mdrun_mpi+' '+mdrun_o+' -maxh 0.3 -c final.pdb >> mdrun-prod.log 2>&1')

    # Dump pdb file from Sim
    utils.sys_call('echo 0 | ' + gromacs+' trjconv -f traj_comp.xtc -o '+pSim.name+'_cg000_f.pdb -sep')

    # Cleanup temp
    utils.sys_call('rm \#* || true')


def backmapp_cg_sim(pSim, simpath, frame):
    """Backmapp a previusly run CG sim to AA resolution."""
    # @TODO HII fix so it uses specific frame from CG simulation
    LOGGER.info("---> Backmapp CG sim to AA in dir {}".format(simpath))
    os.chdir(simpath)

    # Find CG files
    # @TESTING - for now this is after the very short test CG simulations
    #  - here we should be able to specify what time/frame to use
    cgLastFrame = "confout-2ns.gro"
    cgTopFile = "system.top"
    aaTopFile = "system-init.top"
    cgLastFrame = particlesim.get_sim_sub_dir(pSim, "cg-prod")+"/" + cgLastFrame
    cgTopFile = particlesim.get_sim_sub_dir(pSim, "cg-prod")+"/" + cgTopFile
    utils.check_file(cgLastFrame)
    utils.check_file(cgTopFile)
    utils.get_file(particlesim.dirCHARMM+'/system-aa-default.top', aaTopFile)

    # @WARNING assumes x3 lines of solvent
    lines = np.sum(pSim.lipidLowerRatioArray != 0) + np.sum(pSim.lipidUpperRatioArray != 0)
    utils.sys_call('tail -n '+str(lines+3)+' '+cgTopFile+' | head -n ' + str(lines) + ' >> '+aaTopFile)

    # Run the actual backmapping and intal growing in of the AA force field
    utils.sys_call(particlesim.dirCHARMM+'/initram.sh -f '+cgLastFrame+' -p '+aaTopFile +
                   ' -o confout-AA-init.gro -po system.top -from martini -to charmm36')

    # Cleanup temp files
    utils.sys_call('rm \#*  || true')
    utils.sys_call('rm step*.pdb  || true')


def equilibrate_aa_sim(pSim, simpath, frame):
    """Equilibrate a AA simulation that has alreaddy been backmapped."""
    # @TODO HII fix so it uses specific frame from CG simulation
    LOGGER.info("---> Equilibrate AA sim in dir {}".format(simpath))
    os.chdir(simpath)

    # Get backmapped AA files
    lastFrame = "confout-AA-init.gro"
    topFile = "system.top"
    utils.get_file(particlesim.get_sim_sub_dir(pSim, "aa-backmapp") + "/" + lastFrame, lastFrame)
    utils.get_file(particlesim.get_sim_sub_dir(pSim, "aa-backmapp") + "/" + topFile, topFile)

    # Minimize and run AA simulation
    utils.sys_call(gromacs+" grompp -f "+particlesim.dirCHARMM+"/md-aa-em.mdp -o md-aa_minimization.tpr -c " +
                   lastFrame+" -p "+topFile+" -maxwarn -1")
    # WARNING here it might be good to use double pressision GROMACS
    utils.sys_call(gromacs+" mdrun -v -deffnm md-aa_minimization")
    utils.sys_call("cp md-aa_minimization.gro md-aa-eq6.0.gro")

    # Make index file
    output = open("ndx_input.txt", 'w')
    output.write("del 1 - 200 \n")
    output.write("r TIP3 | r POT | r CLA \n")
    output.write("name 1 SOL_ION \n")
    output.write("!1 \n")
    output.write("name 2 MEMB \n")
    output.write("q\n")
    output.close()
    utils.sys_call(gromacs+" make_ndx -f md-aa-eq6.0.gro -o index.ndx < ndx_input.txt")  # Make "MEMB" "SOL_ION" index groups

    # Run eq. sims
    command = ("for i in `seq 1 6`; do echo $i; clast=$(($i - 1)); "
               ""+gromacs+" grompp -f "+particlesim.dirCHARMM+"/md-aa-eq6.$i.mdp -o md-aa-eq6.$i.tpr " +
               "-c md-aa-eq6.$clast.gro -r md-aa_minimization.gro -n index.ndx -p "+topFile+" -maxwarn -1 ; "
               ""+gromacs+" mdrun -v -deffnm md-aa-eq6.$i; done")
    utils.sys_call(command)

    # Cleanup temp files
    utils.sys_call('rm \#* || true')


def run_aa_sim(pSim, simpath, frame):
    """Run production AA simulation."""
    # @TODO HII fix so it uses specific frame from CG simulation
    LOGGER.info("---> Run AA sim in dir {}".format(simpath))
    os.chdir(simpath)

    # Get backmapped AA files
    lastFrame = "md-aa-eq6.6.gro"
    topFile = "system.top"
    indexFile = "index.ndx"
    prodmdpFile = "md-aa-production-t310.mdp"
    utils.get_file(particlesim.get_sim_sub_dir(pSim, "aa-eq") + "/" + lastFrame, lastFrame)
    utils.get_file(particlesim.get_sim_sub_dir(pSim, "aa-eq") + "/" + topFile, topFile)
    utils.get_file(particlesim.get_sim_sub_dir(pSim, "aa-eq") + "/" + indexFile, indexFile)
    utils.get_file(particlesim.dirCHARMM + "/" + prodmdpFile, prodmdpFile)

    # Prep production run
    utils.sys_call(gromacs+" grompp -f "+prodmdpFile+" -o topol.tpr -c "+lastFrame+" -n index.ndx -p system.top -maxwarn -1")

    # Run production
    # @TESTING for now only do 100 ps this should be 100-5000-inf ns and put somewhere so it's easy for the simulation
    #  controler to stop / restart / extend
    utils.sys_call(gromacs+' mdrun -v -c confout-100ps.gro')

    # Cleanup temp files
    utils.sys_call('rm \#* || true')


def save_insane_lipid_counts():
    """Test function to read back up the simulation lipid counts - after insane call."""
    consArray = np.zeros([5,5,14])
    lipidTypes = ["POPX", "PAPC", "POPE", "DIPE", "DPSM", "CHOL", "POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"]
    lipidOrd = [8, 9, 10, 11, 12, 13, 0, 1, 2, 3, 4, 5, 6, 7]
    f = open("system.top")
    data = f.readlines()
    index = 0
    for i in range(len(data)):
        if lipidTypes[0] in data[i].strip().split(" ")[0]:
            index = i
            break
    #print(index) # first line with lipid

    for i in range(5):
        for j in range(5):
            for k in range(14):
                if lipidTypes[k] in data[index].strip().split(" ")[0]:
                    consArray[i,j,lipidOrd[k]] = int(data[index].strip().split(" ")[-1])
                    index += 1
                else: # lipid missing
                    consArray[i,j,lipidOrd[k]] = 0

    np.save('concentrations_cg.npy', consArray)
    f.close()


def correctPBCjump(deltaPost, box_size):
    """Correct for possible PBC wrapping of a bead - that would cause a jump in delta cords."""
    assert len(deltaPost) == len(box_size)
    for i in range(len(deltaPost)):
        chalfsize = box_size[i] / 2.0
        if deltaPost[i] > chalfsize:
            deltaPost[i] = deltaPost[i] - box_size[i]
        if deltaPost[i] <= -chalfsize:
            deltaPost[i] = deltaPost[i] + box_size[i]
    return deltaPost


def generate_probabilistic_patch(concentrations, particle_count):
    """
    Convert cell concentrations to particle counts.

    :params concentrations: A list of concentrations.
    :params particle_count: Number of particles contained in the cell.
    :returns: A list the same length as concentrations
    """
    # Compute the raw counts of particles in relation to their relative
    # concentration in the whole cell.
    raw_counts = concentrations*(particle_count/np.sum(concentrations))
    # Make an initial approximating distribution by rounding to the
    # largest integer number of particles for each specie that is less
    # than or equal to the prescribed (non-integral) number.
    int_counts = np.floor(raw_counts)
    int_sum = np.sum(int_counts)

    # Compute cumulative distribution function of remaining fractions
    cdf = np.cumsum(raw_counts - int_counts)
    if cdf[-1] == 0: # if they happen to match exactly no adjustment needed
        return int_counts.astype(int)
    cdf = cdf/cdf[-1]

    # Sample this distribution until we have as many
    # particles as needed. This ensures that:
    #   a) The number of particles for each specie is
    #      at least the largest whole number representing
    #      the correct particle count for that specie.
    #   b) on average (e.g. by sampling many times),
    #      the distribution is the prescribed one, i.e.
    #      the sampling has the correct expectation value.
    while int_sum < particle_count:
        j = 0
        r = np.random.random()
        while r > cdf[j]:
            j = j+1

        int_counts[j] += 1
        int_sum += 1

    # Convert count array to integer type before returning
    return int_counts.astype(int)


def creat_cg_patch(specifications, simnum, mini_mummi = False):
    """Creates a default mix8 one KRAS patch, then modifies based on specifications."""

    # Average Mix8b - no sub-composision
    # Always use Martini 3 lipid names Martini 2 are still here for backwards compatibility 
    if True:
        # M3
        lipidTypes = ["POPX", "POPC", "PAPC", "POPE", "DLPE", "SSM", "PAPS", "SAP6", "CHOL"]
        lipidUpperRatioArray = [[243,   0, 121, 20, 61, 242, 0, 0, 350]]
        lipidLowerRatioArray = [[0,   139, 75, 54, 161, 108, 161, 22, 243]]
    else:
        # M2 
        lipidTypes = ["POPX", "POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"]
        lipidUpperRatioArray = [[243,   0, 121, 20, 61, 242, 0, 0, 313]]
        lipidLowerRatioArray = [[0,   139, 75, 54, 161, 108, 161, 22, 280]]
    proteinUnits = []

    cProteinDict = {}
    cProteinDict["id"] = 0
    if "RAS4A_RAF" in simnum:
        # RAS4A-RAF
        cProteinDict["type"] = "RAS4A-RAF"
        cProteinDict["angle"] = 0
        cProteinDict["beads"] = -1
        cProteinDict["state"] = "N/A"
        cProteinDict["position"] =  [15,15]
        cProteinDict["structure_pre_assigned"] = True
        cProteinDict["structure"] = "KRAS4A_RAF_shift.gro"
        cProteinDict["charge"] = -1
    elif "RAS4A" in simnum:
        # RAS4A
        cProteinDict["type"] = "RAS4A-ONLY"
        cProteinDict["beads"] = -1
        cProteinDict["state"] = "N/A"
        cProteinDict["position"] =  [15,15]
        cProteinDict["structure_pre_assigned"] = True
        cProteinDict["structure"] = "KRAS4A_shift.gro"
        cProteinDict["charge"] = -2
    elif "RAS_RAF" in simnum:
        if mini_mummi:
            # Update for mini-MuMMI M3
            # RAS-RAF (with regular KRAS4b)
            cProteinDict["type"] = "RAS-RAF"
            cProteinDict["angle"] = "" # not set 
            cProteinDict["beads"] = -1
            cProteinDict["state"] = "N/A"
            cProteinDict["position"] =  [10,10]
            cProteinDict["structure_pre_assigned"] = True
            cProteinDict["structure"] = "n/a" # will set later 
            cProteinDict["charge"] = -1
        else:
            # RAS-RAF (with regular KRAS4b)
            cProteinDict["type"] = "RAS-RAF"
            cProteinDict["angle"] = 0
            cProteinDict["beads"] = -1
            cProteinDict["state"] = "N/A"
            cProteinDict["position"] =  [15,15]
            cProteinDict["structure_pre_assigned"] = True
            cProteinDict["structure"] = "RAS-RAF_state_ma_1.gro"
            cProteinDict["charge"] = -1
    elif "dimercomp" in simnum:
        # Full complex RAS-RAF (m1) 14-3-3^2 RAS-RAF (m2)
        cProteinDict["type"] = "dimercomp"
        cProteinDict["angle"] = "" # not set 
        cProteinDict["beads"] = -1
        cProteinDict["state"] = "N/A"
        cProteinDict["position"] =  [18,18]
        cProteinDict["structure_pre_assigned"] = True
        cProteinDict["structure"] = "n/a" # will set later 
        cProteinDict["charge"] = 0
    else:
        # RAS-only (regular KRAS4b)
        cProteinDict["type"] = "RAS-ONLY"
        cProteinDict["beads"] = -1
        cProteinDict["state"] = "N/A"
        cProteinDict["position"] =  [15,15]
        cProteinDict["structure_pre_assigned"] = True
        cProteinDict["structure"] = "RAS-ONLY_state_a_100.gro"
        cProteinDict["charge"] = -2

    proteinUnits.append(cProteinDict)

    # mini-MuMMI 20x20 patch in M3 
    if mini_mummi:
        patchSim = particlesim.ParticleSim(lipidTypes, lipidUpperRatioArray, lipidLowerRatioArray, lipidSuPatches=[1, 1], proteinUnits=proteinUnits, asym=-65, simName=simnum, lipidFullAsignment=False, mini_mummi = True)
    else:
        patchSim = particlesim.ParticleSim(lipidTypes, lipidUpperRatioArray, lipidLowerRatioArray, lipidSuPatches=[1, 1], proteinUnits=proteinUnits, asym=-118, simName=simnum, lipidFullAsignment=False)
    #patchSim.cgnt = 8 not supportet if you also specified in commond input

    # Alter default based on specifications
    for cAttr in specifications.split():
        [key,value] = cAttr.split('=')
        setattr(patchSim, key, value)

    return patchSim

def convert_gc_patch_to_cg_patch(p, simnum, wkpath):
    """Converts a GC patch to a CG patch."""

    if not p:
        LOGGER.error(f"The patch {simnum} {wkpath} appears to be empty")
        raise ValueError(f"The patch {simnum} {wkpath} appears to be empty")

    # Print patch info
    LOGGER.info(f"Macro patch is: {str(p)}")
    LOGGER.info("Patch id            = {}".format(p.id))
    LOGGER.info("Protein bead ids    = {}".format(p.protein_ids))
    LOGGER.info("Protein bead states = {}".format(p.protein_states))
    LOGGER.info("Protein bead pos.   = {}".format(p.protein_positions))
    LOGGER.info("Macro simname       = {}".format(p.config.simname))
    LOGGER.info("Macro filenames     = {}".format(p.config.filenames))
    LOGGER.info("Macro simname       = {} / {}".format(p.config.simtime, p.config.tunit))

    # Remove any suppriusly negative densities - this can happen in the macro model due to noise term
    p.concentrations = np.maximum(p.concentrations, 0)
    # For testing of lipid consentrations
    #np.save("concentrations_macro_all.npy", p.concentrations)

    # @TODO remove currently lconcentrations_full is used this should be removed here as well as in particle sim object.
    # Convert patch from native 37x37 to 5x5 patch size
    # @TOOD Harsh is looking at this now
    #lconcentrations = p.subsample()
    #lconcentrations = p._subsample_mean(p.concentrations, 30.0, 5)

    # @FIX - this was 73 x 73 now 72 x 72 so not working ???
    #lconcentrations = p.subsample_intg()
    lconcentrations = p.concentrations  # we don't really use this take it out soon
    # For testing of lipid consentrations
    #np.save("concentrations_macro.npy", lconcentrations)
    LOGGER.info("Lipid nat. grid     = {}".format(p.concentrations.shape))
    LOGGER.info("Lipid subs.grid     = {}".format(lconcentrations.shape))

    # Convert patch to CG input
    #   Info Patch.LIPID_NAMES =
    #   ['nInner_POPC', 'nInner_PAPC', 'nInner_POPE', 'nInner_DIPE', 'nInner_DPSM', 'nInner_PAPS', 'nInner_PAP6', 'nInner_CHOL',
    #    'nOuter_POPC', 'nOuter_PAPC', 'nOuter_POPE', 'nOuter_DIPE', 'nOuter_DPSM', 'nOuter_CHOL']
    lipidTypes = ["POPX", "POPC", "PAPC", "POPE", "DLPE", "SSM", "PAPS", "SAP6", "CHOL"]

    lInner = lconcentrations[:, :, 0:8]
    lOuter = lconcentrations[:, :, 8:14]
    sumOuter = np.sum(lOuter) / (lconcentrations.shape[0]*lconcentrations.shape[1])
    sumInner = np.sum(lInner) / (lconcentrations.shape[0]*lconcentrations.shape[1])
    asym = np.rint(2809 * (1 - (sumInner / sumOuter))).astype(int)
    # For testing of lipid consentrations
    #np.save("concentrations_asym.npy", asym)

    # Save full lipid type asignmet patch for insane to read
    lipid_types_inner = ["POPC", "PAPC", "POPE", "DLPE", "SSM", "PAPS", "SAP6", "CHOL"]
    lipid_types_outer = ["POPX", "PAPC", "POPE", "DLPE", "SSM", "CHOL"]
    #p.concentrations.shape # (37, 37, 14)  (for new patches)
    #lconcentrations_full = pdcglob.reinterp(37, 37, 14, p.concentrations, 40, 40)
    # To fit protein placement we need to use the transpose of x,y columns
    lconcentrations_full = pdcglob.reinterp(100, 100, 14, p.concentrations.transpose((1,0,2)), 53, 53)
    lInner_full = lconcentrations_full[:, :, 0:8]
    lOuter_full = lconcentrations_full[:, :, 8:14]
    sumOuter_full = np.sum(lOuter_full) / (lconcentrations_full.shape[0]*lconcentrations_full.shape[1])
    sumInner_full = np.sum(lInner_full) / (lconcentrations_full.shape[0]*lconcentrations_full.shape[1])

    # @WARNING should not be asym_full always 0
    asym_full = np.rint(2809 * (1 - (sumInner_full / sumOuter_full))).astype(int)
    LOGGER.info("Both asym(s) shoudl be the same asym {} and asym_full {} outer {} inner {}".format(asym, asym_full, sumOuter_full, sumInner_full))
    lipid_counts_inner = pdcglob.pdcglob(53, 53, 8, lInner_full)
    lipid_counts_outer = pdcglob.pdcglob(53, 53, 6, lOuter_full)
    lipid_counts_file = os.path.join(wkpath, "lipid_counts.npz")
    np.savez_compressed(lipid_counts_file, lipid_counts_inner = lipid_counts_inner, lipid_counts_outer = lipid_counts_outer, asym = asym_full, lipid_types_inner = lipid_types_inner, lipid_types_outer = lipid_types_outer)
    LOGGER.info("Full patch lipid type assignment saved to {}, asym {}, shape {}".format(lipid_counts_file, asym_full, lconcentrations_full.shape))


    if (lInner.min() < 0 or lOuter.min() < 0):
        error = "---> Negative lipid consentrations found "
        LOGGER.error(error)
        raise ValueError(error)

    if (asym > 600 or asym < -600):
        error = "---> Bilayer asymmetry is all to high asym = " + str(asym)
        LOGGER.error(error)
        raise ValueError(error)

    # if p.protein_ids not set make it a empyt list so len = 0
    if p.protein_ids is None:
        p.protein_ids = []

    # @WARNING this should not be here - remove after C3 - a temp fix to not build patches with more than 4 protein - used in last part of C3
    if (len(p.protein_ids) > 4):
        error = "---> To many protein beads in path - stop build. Current p.protein_ids count = " + str(len(p.protein_ids))
        LOGGER.error(error)
        raise ValueError(error)

    # Adm-C2 if 25 subpatches (25*64=1600, not we have 25*112 = 2800 (should be 2809) 
    # this is not exact but we dont' use this anymore so not an issue 
    # @TODO change this to based on 64 lipids per subpatch (to explore rounding)
    # @TODO depdning on variations in cons from PF maybe to "proper" pre rounding here ???
    #old_lOuter = np.rint(64 * lOuter / sumOuter).astype(int)
    #old_lInner = np.rint(64 * lInner / sumInner).astype(int)

    # Compute a probablistic construction of the particles in a cell.
    # We use a cumlative distribution function to
    # TODO: We should really assert or something here to verify that the
    # inner and outer leaflets are the same size.
    # Iterate through each grid point and convert to the number of particles
    # based on the probability.
    for i in range(0, lInner.shape[0]):
        for j in range(0, lInner.shape[1]):
            lInner[i, j] = generate_probabilistic_patch(lInner[i, j], 112)
    lInner = lInner.astype(int)
    for i in range(0, lOuter.shape[0]):
        for j in range(0, lOuter.shape[1]):
            lOuter[i, j] = generate_probabilistic_patch(lOuter[i, j], 112)
    lOuter = lOuter.astype(int)

    # @TODO Use this to explor lipid consetnrations - excel magic ( a 0 values will be left out of the list but 1/1000 is ok)
    # not needed for production and makes other consentrations a little strange
    # lOuter[lOuter == 0] = 1
    # lInner[lInner == 0] = 1

    # For testing of lipid consentrations
    # @TODO add this to a extra debug / analysis flag
    #saveCons = np.zeros(lconcentrations.shape)
    #saveCons[:, :, 0:8] = lInner
    #saveCons[:, :, 8:14] = lOuter
    #np.save("concentrations_macro_int_"+simnum+".npy", saveCons)

    # Convet into the right shape (adding empty lipids as 0 and chansing 5x5 to 25)
    # lipidUpperRatioArray = [[243,   0, 121, 20, 61, 242, 0, 0, 313]]
    # lipidLowerRatioArray = [[0,   139, 75, 54, 161, 108, 161, 22, 280]]
    lipidUpperRatioArray = np.insert(np.insert(np.insert(lOuter, 1, 0, axis=2), 6, 0, axis=2), 7, 0, axis=2)
    lipidLowerRatioArray = np.insert(lInner, 0, 0, axis=2)
    lipidUpperRatioArray = np.reshape(lipidUpperRatioArray, (lconcentrations.shape[0]*lconcentrations.shape[1], len(lipidTypes)))
    lipidLowerRatioArray = np.reshape(lipidLowerRatioArray, (lconcentrations.shape[0]*lconcentrations.shape[1], len(lipidTypes)))

    # Convert protein bead list to protein "unit" list - where each unit is inserted a one piece later in the pipeline
    # For each protein unit - get inital structure, insertion properties (location, angle etc), and charge

    LOGGER.info("Protein bead ids    = {}".format(p.protein_ids))
    LOGGER.info("Protein bead states = {}".format(p.protein_states))
    LOGGER.info("Protein bead pos.   = {}".format(p.protein_positions))
    LOGGER.info("Protein complex ids = {}".format(p.complex_ids))
    #LOGGER.info("Protein comple dis. = {}".format(p.complex_dists)) # removed from storage!

    # Initial shift not needed anymore as first beads is at center pixel (not 0,0 but within 1/2 pixel of that)
    # posArray = p.protein_positions - p.protein_positions[0]  # Center at 0,0

    # Now I assume they are all in the right box (0,0 to max,max)
    posArray = p.protein_positions

    # if p.complex_ids not set make it a empyt list so len = 0
    if p.complex_ids is None:
        p.complex_ids = []

    # Change from RAS4B to RAS4A for C4 hearling run (code below in if for both RAS-only and RAS-RAF)
    #  If False this will be a normal 4B setup if True change to 4A
    c4_healing_convert_to_4A = False
    if (len(simnum) > 0 and "4A" in simnum):
        c4_healing_convert_to_4A = True
        LOGGER.warning("Convert 4B to 4A")

    # @TODO fix this? this is an evil hack to get the protein unites form current complexes - Harsh make RAS also complex?
    # works as we only have x2 types of proteins
    proteinUnits = []
    proteinUnitCount = len(p.protein_ids) - len(p.complex_ids) # Currently each complex has x2 beads
    LOGGER.info("Protein unit count  = {}".format(proteinUnitCount))
    remainig_ids = p.protein_ids
    for i in range(proteinUnitCount):
        cProteinDict = {}
        cProteinDict["id"] = i

        # check if current bead belongs to a complex
        if p.protein_states[i] == 'rasraf2D':
            cComplex = []
            for x in p.complex_ids:
                if remainig_ids[0] in x:
                    cComplex = x
                    break
            cProteinDict["beads"] = cComplex
            cProteinDict["state"] = 'rasraf2D'
            cProteinDict["type"] = "RAS-RAF"
            cProteinDict["structure"] = "N/A"
            cProteinDict["charge"] = -1

            # RAS-RAF is places usign RAS farnesyl x,y pos
            cPosRAS = posArray[0][0:2] # @WARNING there are x3 cords last used for GCMM analysis data not relevant here and not z value
            cProteinDict["position"] = cPosRAS
            assert cProteinDict["position"].shape[0] == 2

            LOGGER.info("Protein unit {} adding rasraf2D complex for bead ids {} full info: {}".format(i, cComplex, str(cProteinDict)))

        elif (len(p.complex_ids) > 0 and any(remainig_ids[0] in string for string in p.complex_ids)):
            # Adding a RAS-RAF protein (either 4B or 4A)
            cComplex = []
            for x in p.complex_ids:
                if remainig_ids[0] in x:
                    cComplex = x
                    break
            ras_bead_id = cComplex[0]
            raf_bead_id = cComplex[1]
            # Flip beads if RAS/RAF bead order was reverse
            if "RAF" in p.protein_states[np.where(p.protein_ids == ras_bead_id)][0]:
                ras_bead_id = cComplex[1]
                raf_bead_id = cComplex[0]
            cProteinDict["beads"] = cComplex
            ras_state = p.protein_states[np.where(p.protein_ids == ras_bead_id)][0]
            cProteinDict["state"] = ras_state[0] + ras_state[-1]

            # Get the angle orientation of RAS-RAF
            cPosRAS = posArray[np.where(p.protein_ids == ras_bead_id)][0][0:2] # @WARNING there are x3 cords last used for GCMM analysis data not relevant here and not z value
            cPosRAF = posArray[np.where(p.protein_ids == raf_bead_id)][0][0:2] # @WARNING there are x3 cords last used for GCMM analysis data not relevant here and not z value

            deltaPost = cPosRAF - cPosRAS
            deltaPost = correctPBCjump(deltaPost, [40, 40]) # in current GCMM schema beads are allways within PBC here we don't want that but the RAF that is closest to RAS
            rafAngle = math.degrees(math.atan2(deltaPost[1], deltaPost[0])) % 360
            # atan2 returns in radians -pi to pi and takes in y,x so flip, convert and set to 0-360
            cProteinDict["angle"] = rafAngle
            LOGGER.info("RAF angle RAS.pos {} RAF.pos {}, delta.pos {}, angle {}".format(cPosRAS, cPosRAF, deltaPost, rafAngle))

            # RAS-RAF is places usign RAS farnesyl x,y pos
            cProteinDict["position"] = cPosRAS
            assert cProteinDict["position"].shape[0] == 2

            if "RAS4A" in p.protein_states[np.where(p.protein_ids == ras_bead_id)][0] or c4_healing_convert_to_4A: # c4_healing_convert_to_4A is to convert 4B to 4A for specific runs e.g. C4 hearling runs
                # Now we are adding a RAS4A-RAF (with RAS4A)
                cProteinDict["type"] = "RAS4A-RAF"
                cProteinDict["structure"] = "N/A"
                cProteinDict["charge"] = -1
                LOGGER.info("Protein unit {} adding RAS4A-RAF (4A RAS) complex for bead ids {} full info: {}".format(i, cComplex, str(cProteinDict)))
            else:
                # Now we are adding a RAS-RAF (with RAS4B)
                cProteinDict["type"] = "RAS-RAF"
                cProteinDict["structure"] = "N/A"
                cProteinDict["charge"] = -1
                LOGGER.info("Protein unit {} adding RAS-RAF (4B RAS) complex for bead ids {} full info: {}".format(i, cComplex, str(cProteinDict)))
            # Remove both complex beads from list
            remainig_ids = remainig_ids[~np.isin(remainig_ids, cComplex)]

        else:
            # Adding a RAS protein (either 4B or 4A)
            cProteinDict["beads"] = [remainig_ids[0]]
            cProteinDict["state"] = p.protein_states[np.where(p.protein_ids == remainig_ids[0])][0][-1]
            cProteinDict["position"] = posArray[np.where(p.protein_ids == remainig_ids[0])][0]
            cProteinDict["position"] = cProteinDict["position"][0:2] # @WARNING there are x3 cords used for GCMM analysis data not relevant here and not z value
            assert cProteinDict["position"].shape[0] == 2

            if "RAS4A" in p.protein_states[np.where(p.protein_ids == remainig_ids[0])][0] or c4_healing_convert_to_4A: # c4_healing_convert_to_4A is to convert 4B to 4A for specific runs e.g. C4 hearling runs
                # Now we are adding a RAS4A-only
                cProteinDict["type"] = "RAS4A-ONLY"
                cProteinDict["structure"] = "N/A"
                # Test for inital structure library
                cProteinDict["charge"] = -2
                LOGGER.info("Protein unit {} adding RAS4A-ONLY (4A RAS) for bead id {} full info: {}".format(i, remainig_ids[0], str(cProteinDict)))
            else:
                # Now we are adding a RAS-only (with RAS4B)
                cProteinDict["type"] = "RAS-ONLY"
                cProteinDict["structure"] = "N/A"
                cProteinDict["charge"] = -2
                LOGGER.info("Protein unit {} adding RAS-ONLY (4B RAS) for bead id {} full info: {}".format(i, remainig_ids[0], str(cProteinDict)))
            # remove current bead from list
            remainig_ids = remainig_ids[~np.isin(remainig_ids, remainig_ids[0])]

        proteinUnits.append(cProteinDict)


    # Print CG patch input
    LOGGER.info("CG_patch input:")
    LOGGER.info("macro patch summ_outer = {}".format(sumOuter))
    LOGGER.info("macro patch summ_inner = {}".format(sumInner))
    LOGGER.info("patch asym             = {}".format(asym))
    LOGGER.info("Lipids outer cons.     = {}".format(str(lipidUpperRatioArray)))
    LOGGER.info("Lipids inner cons.     = {}".format(str(lipidLowerRatioArray)))
    LOGGER.info("Protein units all      = {}".format(proteinUnits))

    # Add simnum after name, if provided (shoud be empty in MuMMI runs)
    simName=p.id
    if len(simnum) > 0:
        simName += "_" + simnum

    patchSim = particlesim.ParticleSim(lipidTypes, lipidUpperRatioArray, lipidLowerRatioArray, [lconcentrations.shape[0], lconcentrations.shape[1]],
                                       proteinUnits, asym=-asym, simName=simName, lipidFullAsignment=True)
    return patchSim


def convert_pf_patch_to_cg_patch(p, simnum):
    """Converts a pfpatch p to a CG patch."""

    # Print patch info
    LOGGER.info("Macro patch is:")
    LOGGER.info(str(p))
    LOGGER.info("Patch id            = {}".format(p.id))
    LOGGER.info("Protein bead ids    = {}".format(p.protein_ids))
    LOGGER.info("Protein bead states = {}".format(p.protein_states))
    LOGGER.info("Protein bead pos.   = {}".format(p.protein_positions))
    LOGGER.info("Macro simname       = {}".format(p.config.simname))
    LOGGER.info("Macro filenames     = {}".format(p.config.filenames))
    LOGGER.info("Macro simname       = {} / {}".format(p.config.simtime, p.config.tunit))

    # Remove any suppriusly negative densities - this can happen in the macro model due to noise term
    p.concentrations = np.maximum(p.concentrations, 0)
    # For testing of lipid consentrations
    #np.save("concentrations_macro_all.npy", p.concentrations)

    # Convert patch from native 37x37 to 5x5 patch size
    # @TOOD Harsh is looking at this now
    #lconcentrations = p.subsample()
    #lconcentrations = p._subsample_mean(p.concentrations, 30.0, 5)
    #lconcentrations = p.subsample_intg()
    lconcentrations = p.concentrations
    # For testing of lipid consentrations
    #np.save("concentrations_macro.npy", lconcentrations)

    LOGGER.info("Lipid nat. grid     = {}".format(p.concentrations.shape))
    LOGGER.info("Lipid subs.grid     = {}".format(lconcentrations.shape))

    # Get protein bead posisions - these should all be within the patch
    # Get RAS posision within the patch - RAS is now back to be patch centric
    #localRasPos = p.rasPositions - p.extents[0]
    #LOGGER.info("RAS pos         = {}".format(p.rasPositions))
    #LOGGER.info("RAS local pos   = {}".format(localRasPos))
    #LOGGER.info("Protein bead pos.   = {}".format(p.protein_positions))

    # Convert patch to CG input
    #   Info Patch.LIPID_NAMES =
    #   ['nInner_POPC', 'nInner_PAPC', 'nInner_POPE', 'nInner_DIPE', 'nInner_DPSM', 'nInner_PAPS', 'nInner_PAP6', 'nInner_CHOL',
    #    'nOuter_POPC', 'nOuter_PAPC', 'nOuter_POPE', 'nOuter_DIPE', 'nOuter_DPSM', 'nOuter_CHOL']
    lipidTypes = ["POPX", "POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"]

    lInner = lconcentrations[:, :, 0:8]
    lOuter = lconcentrations[:, :, 8:14]
    sumOuter = np.sum(lOuter) / (lconcentrations.shape[0]*lconcentrations.shape[1])
    sumInner = np.sum(lInner) / (lconcentrations.shape[0]*lconcentrations.shape[1])
    asym = np.rint(1600 * (1 - (sumInner / sumOuter))).astype(int)
    # For testing of lipid consentrations
    #np.save("concentrations_asym.npy", asym)

    if (lInner.min() < 0 or lOuter.min() < 0):
        error = "---> Negative lipid consentrations found "
        LOGGER.error(error)
        raise ValueError(error)

    if (asym > 500 or asym < -500):
        error = "---> Bilayer asymmetry is all to high asym = " + str(asym)
        LOGGER.error(error)
        raise ValueError(error)

    # @WARNING this should not be here - remove after C3 - a temp fix to not build patches with more than 4 protein - used in last part of C3
    if (len(p.protein_ids) > 4):
        error = "---> To many protein beads in patch - stop build. Current p.protein_ids count = " + str(len(p.protein_ids))
        LOGGER.error(error)
        raise ValueError(error)

    # @TODO change this to based on 64 lipids per subpatch (to explore rounding)
    # @TODO depdning on variations in cons from PF maybe to "proper" pre rounding here ???
    #old_lOuter = np.rint(64 * lOuter / sumOuter).astype(int)
    #old_lInner = np.rint(64 * lInner / sumInner).astype(int)

    # Compute a probablistic construction of the particles in a cell.
    # We use a cumlative distribution function to
    # TODO: We should really assert or something here to verify that the
    # inner and outer leaflets are the same size.
    # Iterate through each grid point and convert to the number of particles
    # based on the probability.
    for i in range(0, lInner.shape[0]):
        for j in range(0, lInner.shape[1]):
            lInner[i, j] = generate_probabilistic_patch(lInner[i, j], 64)
    lInner = lInner.astype(int)
    for i in range(0, lOuter.shape[0]):
        for j in range(0, lOuter.shape[1]):
            lOuter[i, j] = generate_probabilistic_patch(lOuter[i, j], 64)
    lOuter = lOuter.astype(int)

    # @TODO Use this to explor lipid consetnrations - excel magic ( a 0 values will be left out of the list but 1/1000 is ok)
    # not needed for production and makes other consentrations a little strange
    # lOuter[lOuter == 0] = 1
    # lInner[lInner == 0] = 1

    # For testing of lipid consentrations
    # @TODO add this to a extra debug / analysis flag
    #saveCons = np.zeros(lconcentrations.shape)
    #saveCons[:, :, 0:8] = lInner
    #saveCons[:, :, 8:14] = lOuter
    #np.save("concentrations_macro_int_"+simnum+".npy", saveCons)

    # Convet into the right shape (adding empty lipids as 0 and chansing 5x5 to 25)
    # lipidUpperRatioArray = [[243,   0, 121, 20, 61, 242, 0, 0, 313]]
    # lipidLowerRatioArray = [[0,   139, 75, 54, 161, 108, 161, 22, 280]]
    lipidUpperRatioArray = np.insert(np.insert(np.insert(lOuter, 1, 0, axis=2), 6, 0, axis=2), 7, 0, axis=2)
    lipidLowerRatioArray = np.insert(lInner, 0, 0, axis=2)
    lipidUpperRatioArray = np.reshape(lipidUpperRatioArray, (lconcentrations.shape[0]*lconcentrations.shape[1], len(lipidTypes)))
    lipidLowerRatioArray = np.reshape(lipidLowerRatioArray, (lconcentrations.shape[0]*lconcentrations.shape[1], len(lipidTypes)))

    # Save full lipid type asignmet patch for insane to read
    lipid_types_inner = ["POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"]
    lipid_types_outer = ["POPX", "PAPC", "POPE", "DIPE", "DPSM", "CHOL"]
    #p.concentrations.shape # (37, 37, 14)  (for new patches)
    #lconcentrations_full = pdcglob.reinterp(37, 37, 14, p.concentrations, 40, 40)
    # To fit protein placement we need to use the transpose of x,y columns
    lconcentrations_full = pdcglob.reinterp(37, 37, 14, p.concentrations.transpose((1,0,2)), 40, 40)
    lInner = lconcentrations_full[:, :, 0:8]
    lOuter = lconcentrations_full[:, :, 8:14]
    sumOuter = np.sum(lOuter) / (lconcentrations_full.shape[0]*lconcentrations_full.shape[1])
    sumInner = np.sum(lInner) / (lconcentrations_full.shape[0]*lconcentrations_full.shape[1])
    asym_full = np.rint(1600 * (1 - (sumInner / sumOuter))).astype(int)
    LOGGER.info("Both asym(s) shoudl be the same asym {} and asym_full {}".format(asym, asym_full))
    lipid_counts_inner = pdcglob.pdcglob(40, 40, 8, lInner)
    lipid_counts_outer = pdcglob.pdcglob(40, 40, 6, lOuter)
    lipid_counts_file = os.path.join(wkpath, "lipid_counts.npz")
    np.savez_compressed(lipid_counts_file, lipid_counts_inner = lipid_counts_inner, lipid_counts_outer = lipid_counts_outer, asym = asym_full, lipid_types_inner = lipid_types_inner, lipid_types_outer = lipid_types_outer)
    LOGGER.info("Full patch lipid type assignment saved to {}, asym {}, shape".format(lipid_counts_file, asym_full, lconcentrations_full.shape))

    # Convert protein bead list to protein "unit" list - where each unit is inserted a one piece later in the pipeline
    # For each protein unit - get inital structure, insertion properties (location, angle etc), and charge

    LOGGER.info("Protein bead ids    = {}".format(p.protein_ids))
    LOGGER.info("Protein bead states = {}".format(p.protein_states))
    LOGGER.info("Protein bead pos.   = {}".format(p.protein_positions))
    LOGGER.info("Protein types       = {}".format(p.type_c3()))
    LOGGER.info("Protein complex ids = {}".format(p.complex_ids))
    #LOGGER.info("Protein comple dis. = {}".format(p.complex_dists)) # removed from storage!

    # Initial shift not needed anymore as first beads is at center pixel (not 0,0 but within 1/2 pixel of that)
    # posArray = p.protein_positions - p.protein_positions[0]  # Center at 0,0

    # Shift from box center at 0,0 to box lower/left edge at 0,0 as used in the CG simulations
    posArray = []
    if len(p.protein_positions) > 0:
        box_dim = 31.0 # @WARNING hared coded, also in particesim.py
        half_box_dim = box_dim / 2.0 # Get box center to move protein 0,0 to box edge
        posArray = p.protein_positions # First bead in the patche should be RAS and now centerd (close to 0,0)
        posArray = posArray + [half_box_dim, half_box_dim]
        if (posArray.max() > 31 or posArray.min() < 0):
            error = "---> Protein posisions are outside the CG box size; pos are: " + np.array2string(krasPosArray)
            LOGGER.error(error)
            raise ValueError(error)
        LOGGER.info("Protein bead pos ce = {}".format(posArray))

    # Change from RAS4B to RAS4A for C4 hearling run (code below in if for both RAS-only and RAS-RAF)
    #  If False this will be a normal 4B setup if True change to 4A
    c4_healing_convert_to_4A = True

    # @TODO fix this? this is an evil hack to get the protein unites form current complexes - Harsh make RAS also complex?
    # works as we only have x2 types of proteins
    proteinUnits = []
    proteinUnitCount = len(p.protein_ids) - len(p.complex_ids) # Currently each complex has x2 beads
    LOGGER.info("Protein unit count  = {}".format(proteinUnitCount))
    remainig_ids = p.protein_ids
    for i in range(proteinUnitCount):
        cProteinDict = {}
        cProteinDict["id"] = i

        # check if current bead belongs to a complex
        if (len(p.complex_ids) > 0 and any(remainig_ids[0] in string for string in p.complex_ids)):
            cComplex = []
            if c4_healing_convert_to_4A: # Convert 4B to 4A for C4 hearling runs
                # Now we are adding a RAS4A-RAF
                cProteinDict["type"] = "RAS4A-RAF"
                for x in p.complex_ids:
                    if remainig_ids[0] in x:
                        cComplex = x
                        break
                ras_bead_id = cComplex[0]
                raf_bead_id = cComplex[1]
                # Flip beads if RAS/RAF bead order was reverse
                if "RAF" in p.protein_states[np.where(p.protein_ids == ras_bead_id)][0]:
                    ras_bead_id = cComplex[1]
                    raf_bead_id = cComplex[0]
                cProteinDict["beads"] = cComplex
                ras_state = p.protein_states[np.where(p.protein_ids == ras_bead_id)][0]
                cProteinDict["state"] = ras_state[0] + ras_state[-1]

                cPosRAS = posArray[np.where(p.protein_ids == ras_bead_id)][0]
                cPosRAF = posArray[np.where(p.protein_ids == raf_bead_id)][0]
                deltaPost = cPosRAF - cPosRAS
                rafAngle = math.degrees(math.atan2(deltaPost[1], deltaPost[0])) % 360
                # atan2 returns in radians -pi to pi and takes in y,x so flip, convert and set to 0-360
                cProteinDict["angle"] = rafAngle
                LOGGER.info("RAF angle RAS.pos {} RAF.pos {}, delta.pos {}, angle {}".format(cPosRAS, cPosRAF, deltaPost, rafAngle))

                # RAS-RAF is places usign RAS farnesyl x,y pos
                cProteinDict["position"] = cPosRAS
                assert cProteinDict["position"].shape[0] == 2

                cProteinDict["structure"] = "N/A"
                # @TODO change to N/A and - after we fix library for 4A - cProteinDict["structure"] = "N/A"
                cProteinDict["charge"] = -1
                LOGGER.info("Protein unit {} converting 4B to 4A - adding RAS4A-RAF complex for bead ids {} full info: {}".format(i, cComplex, str(cProteinDict)))
            else:
                # Now we are adding a RAS-RAF
                cProteinDict["type"] = "RAS-RAF"
                for x in p.complex_ids:
                    if remainig_ids[0] in x:
                        cComplex = x
                        break
                ras_bead_id = cComplex[0]
                raf_bead_id = cComplex[1]
                # Flip beads if RAS/RAF bead order was reverse
                if "RAF" in p.protein_states[np.where(p.protein_ids == ras_bead_id)][0]:
                    ras_bead_id = cComplex[1]
                    raf_bead_id = cComplex[0]
                cProteinDict["beads"] = cComplex
                ras_state = p.protein_states[np.where(p.protein_ids == ras_bead_id)][0]
                cProteinDict["state"] = ras_state[0] + ras_state[-1]

                # @TODO get the angle orientation of RAS-RAF
                cPosRAS = posArray[np.where(p.protein_ids == ras_bead_id)][0]
                cPosRAF = posArray[np.where(p.protein_ids == raf_bead_id)][0]
                deltaPost = cPosRAF - cPosRAS
                rafAngle = math.degrees(math.atan2(deltaPost[1], deltaPost[0])) % 360
                # atan2 returns in radians -pi to pi and takes in y,x so flip, convert and set to 0-360
                cProteinDict["angle"] = rafAngle
                LOGGER.info("RAF angle RAS.pos {} RAF.pos {}, delta.pos {}, angle {}".format(cPosRAS, cPosRAF, deltaPost, rafAngle))

                # RAS-RAF is places usign RAS farnesyl x,y pos
                cProteinDict["position"] = cPosRAS
                assert cProteinDict["position"].shape[0] == 2

                # Now structure is assigned in placeprotein
                '''
                # Numer of diffrent configurations saved for each RAS-RAF state
                # @WARNING hardcoded  _0, _1, ... _librarySize-1
                librarySize = 5000
                cInt = randint(0, librarySize - 1)
                cStructureName = "{}/{}".format(Naming.dir_res('structures'), Naming.protein_structure(cProteinDict["type"], cProteinDict["state"], cInt))
                LOGGER.info("Get {} structure: state {}, randint {}, filename {}".format(cProteinDict["type"], cProteinDict["state"], cInt, cStructureName))
                if not os.path.isfile(cStructureName):
                    error = "---> Protein structure file {} not found.".format(cStructureName)
                    LOGGER.error(error)
                    raise ValueError(error)
                cProteinDict["structure"] = cStructureName
                '''
                cProteinDict["structure"] = "N/A"
                cProteinDict["charge"] = -1
                LOGGER.info("Protein unit {} adding RAS-RAF complex for bead ids {} full info: {}".format(i, cComplex, str(cProteinDict)))

            # Remove both complex beads from list
            remainig_ids = remainig_ids[~np.isin(remainig_ids, cComplex)]
        else:
            if c4_healing_convert_to_4A: # Convert 4B to 4A for C4 hearling runs
                # Now we are adding a RAS4A-only
                cProteinDict["type"] = "RAS4A-ONLY"
                cProteinDict["beads"] = [remainig_ids[0]]
                cProteinDict["state"] = p.protein_states[np.where(p.protein_ids == remainig_ids[0])][0][3:]
                cProteinDict["position"] =  posArray[np.where(p.protein_ids == remainig_ids[0])][0]
                assert cProteinDict["position"].shape[0] == 2
                cProteinDict["structure"] = "N/A"
                # @TODO change to N/A and - after we fix library for 4A - cProteinDict["structure"] = "N/A"
                cProteinDict["charge"] = -2
                LOGGER.info("Protein unit {} converting 4B to 4A - adding RAS4A-ONLY for bead id {} full info: {}".format(i, remainig_ids[0], str(cProteinDict)))
            else:
                # Now we are adding a RAS-only
                cProteinDict["type"] = "RAS-ONLY"
                cProteinDict["beads"] = [remainig_ids[0]]
                cProteinDict["state"] = p.protein_states[np.where(p.protein_ids == remainig_ids[0])][0][3:]
                cProteinDict["position"] =  posArray[np.where(p.protein_ids == remainig_ids[0])][0]
                assert cProteinDict["position"].shape[0] == 2

                # Now structure is assigned in placeprotein
                '''
                # Numer of diffrent configurations saved for each RAS-ONLY state
                # @WARNING hardcoded  _0, _1, ... _librarySize-1
                librarySize = 5000
                cInt = randint(0, librarySize - 1)
                cStructureName = "{}/{}".format(Naming.dir_res('structures'), Naming.protein_structure(cProteinDict["type"], cProteinDict["state"], cInt))
                LOGGER.info("Get {} structure: state {}, randint {}, filename {}".format(cProteinDict["type"], cProteinDict["state"], cInt, cStructureName))
                if not os.path.isfile(cStructureName):
                    error = "---> Protein structure file {} not found.".format(cStructureName)
                    LOGGER.error(error)
                    raise ValueError(error)
                cProteinDict["structure"] = cStructureName
                '''
                cProteinDict["structure"] = "N/A"

                cProteinDict["charge"] = -2
                LOGGER.info("Protein unit {} adding RAS-ONLY for bead id {} full info: {}".format(i, remainig_ids[0], str(cProteinDict)))
            # remove current bead from list
            remainig_ids = remainig_ids[~np.isin(remainig_ids, remainig_ids[0])]

        #LOGGER.debug("Ids {}".format(remainig_ids))
        proteinUnits.append(cProteinDict)


    # Print CG patch input
    LOGGER.info("CG_patch input:")
    LOGGER.info("macro patch summ_outer = {}".format(sumOuter))
    LOGGER.info("macro patch summ_inner = {}".format(sumInner))
    LOGGER.info("patch asym             = {}".format(asym))
    LOGGER.info("Lipids outer cons.     = {}".format(str(lipidUpperRatioArray)))
    LOGGER.info("Lipids inner cons.     = {}".format(str(lipidLowerRatioArray)))
    LOGGER.info("Protein units all      = {}".format(proteinUnits))

    # Add simnum after name, if provided (shoud be empty in MuMMI runs)
    simName=p.id
    if len(simnum) > 0:
        simName += "_" + simnum

    patchSim = particlesim.ParticleSim(lipidTypes, lipidUpperRatioArray, lipidLowerRatioArray, [lconcentrations.shape[0], lconcentrations.shape[1]],
                                       proteinUnits, asym=-asym, simName=simName, lipidFullAsignment=True)
    return patchSim


def convert_pf_patch_to_cg_patch_old_c1(p):
    """Converts a pfpatch p to a CG patch - This is the old fucntion from Campain 1."""

    # Print patch info
    LOGGER.info("Macro patch is:")
    LOGGER.info(str(p))
    LOGGER.info("Patch id        = {}".format(p.id))
    LOGGER.info("RAS ids         = {}".format(p.rasIds))
    LOGGER.info("RAS states      = {}".format(p.rasStates))
    LOGGER.info("RAS pos         = {}".format(p.rasPositions))
    LOGGER.info("Macro simname   = {}".format(p.config.simname))
    LOGGER.info("Macro filenames = {}".format(p.config.filenames))
    LOGGER.info("Macro simname   = {} / {}".format(p.config.simtime, p.config.tunit))

    # Convert patch to CG input
    #   Info Patch.LIPID_NAMES =
    #   ['nInner_POPC', 'nInner_PAPC', 'nInner_POPE', 'nInner_DIPE', 'nInner_DPSM', 'nInner_PAPS', 'nInner_PAP6', 'nInner_CHOL',
    #    'nOuter_POPC', 'nOuter_PAPC', 'nOuter_POPE', 'nOuter_DIPE', 'nOuter_DPSM', 'nOuter_CHOL']
    lipidTypes = ["POPX", "POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"]

    lInner = p.concentrations[:, :, 0:8]
    lOuter = p.concentrations[:, :, 8:14]
    sumOuter = np.sum(lOuter) / (p.concentrations.shape[0]*p.concentrations.shape[1])
    sumInner = np.sum(lInner) / (p.concentrations.shape[0]*p.concentrations.shape[1])
    asym = np.rint(1600 * (1 - (sumInner / sumOuter))).astype(int)

    # Compute a probablistic construction of the particles in a cell.
    # We use a cumlative distribution function to
    # TODO: We should really assert or something here to verify that the
    # inner and outer leaflets are the same size.
    # Iterate through each grid point and convert to the number of particles
    # based on the probability.
    for i in range(0, lInner.shape[0]):
        for j in range(0, lInner.shape[1]):
            lInner[i, j] = generate_probabilistic_patch(lInner[i, j], 64)
    lInner = lInner.astype(int)
    for i in range(0, lOuter.shape[0]):
        for j in range(0, lOuter.shape[1]):
            lOuter[i, j] = generate_probabilistic_patch(lOuter[i, j], 64)
    lOuter = lOuter.astype(int)

    # @TODO Use this to explor lipid consetnrations - excel magic ( a 0 values will be left out of the list but 1/1000 is ok)
    # not needed for production and makes other consentrations a little strange
    # lOuter[lOuter == 0] = 1
    # lInner[lInner == 0] = 1

    # Convet into the right shape (adding empty lipids as 0 and chansing 5x5 to 25)
    # lipidUpperRatioArray = [[243,   0, 121, 20, 61, 242, 0, 0, 313]]
    # lipidLowerRatioArray = [[0,   139, 75, 54, 161, 108, 161, 22, 280]]
    lipidUpperRatioArray = np.insert(np.insert(np.insert(lOuter, 1, 0, axis=2), 6, 0, axis=2), 7, 0, axis=2)
    lipidLowerRatioArray = np.insert(lInner, 0, 0, axis=2)
    lipidUpperRatioArray = np.reshape(lipidUpperRatioArray, (p.concentrations.shape[0]*p.concentrations.shape[1], len(lipidTypes)))
    lipidLowerRatioArray = np.reshape(lipidLowerRatioArray, (p.concentrations.shape[0]*p.concentrations.shape[1], len(lipidTypes)))

    half_box_dim = 31.0 / 2 # This is fixed here and in particlesim.py
    #krasPosArray = p.rasPositions - p.rasPositions[0]  # Center at 0,0
    krasPosArray = p.rasPositions # KRAS in patches are now centerd at the first KRAS
    krasPosArray = krasPosArray + [half_box_dim, half_box_dim]

    if (krasPosArray.max() > (half_box_dim * 2) or krasPosArray.min() < 0):
        error = "---> KRAS posisions are outside the CG box size; pos are: " + np.array2string(krasPosArray)
        LOGGER.error(error)
        raise ValueError(error)

    # Print CG patch input
    LOGGER.info("CG_patch input:")
    LOGGER.info("macro patch summ_outer = {}".format(sumOuter))
    LOGGER.info("macro patch summ_inner = {}".format(sumInner))
    LOGGER.info("patch asym             = {}".format(asym))
    LOGGER.info("Lipids outer cons.     = {}".format(str(lipidUpperRatioArray)))
    LOGGER.info("Lipids inner cons.     = {}".format(str(lipidLowerRatioArray)))
    LOGGER.info("RAS init pos.          = {}".format(krasPosArray))
    LOGGER.info("RAS states             = {}".format(p.rasStates))

    # Select KRAS structures acording to states - pfPatch state to saved files + random value (e.g. 1-1000)
    # "KRAS-s"+str(cS)+"-r"+str(randint(1, librarySize))+".gro"
    librarySize = 1000  # numer of diffrent configurations saved for each state
    krasStructs = []
    for i in range(len(p.rasStates)):
        cInt = randint(1, librarySize)
        cState = int(p.rasStates[i])
        if (cState != 1 and cState != 2):
            error = "---> KRAS state {} is not supported: ".format(cState)
            LOGGER.error(error)
            raise ValueError(error)
        krasStruct = "{}/random_state{}.{}.KRAS.pbc.mol.gro".format(particlesim.dirKRASStructures, cState, cInt)
        LOGGER.info("Get KRAS structure file: KRAS # {}, in state {} and rand {}, name {}".format(i, cState, cInt, krasStruct))
        krasStructs.append(krasStruct)

    patchSim = particlesim.ParticleSim(lipidTypes, lipidUpperRatioArray, lipidLowerRatioArray, [p.concentrations.shape[0], p.concentrations.shape[1]],
                                       krasStructs, krasPosArray, asym=-asym, simName=p.id, lipidFullAsignment=False)
    return patchSim


def convert_macro_to_cg_patch(ptachID, patch_concentrations, patch_rstates, patch_rpositions):
    """Converts a numpy macro to a CG patch - This is a new function similar to convert_pf_patch_to_cg_patch except creats a sim from the numpy patch data."""

    # Print patch info
    LOGGER.info("Raw macro patch data is:")
    LOGGER.info("Patch id        = {}".format(ptachID))
    LOGGER.info("RAS states      = {}".format(patch_rstates))
    LOGGER.info("RAS pos         = {}".format(patch_rpositions))

    # @TODO fix this - now manually centered
    patch_rpositions = np.array([[0, 0]])
    LOGGER.info("RAS centerd pos         = {}".format(patch_rpositions))

    # Convert patch to CG input
    #   Info Patch.LIPID_NAMES =
    #   ['nInner_POPC', 'nInner_PAPC', 'nInner_POPE', 'nInner_DIPE', 'nInner_DPSM', 'nInner_PAPS', 'nInner_PAP6', 'nInner_CHOL',
    #    'nOuter_POPC', 'nOuter_PAPC', 'nOuter_POPE', 'nOuter_DIPE', 'nOuter_DPSM', 'nOuter_CHOL']
    lipidTypes = ["POPX", "POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"]

    lInner = patch_concentrations[:, :, 0:8]
    lOuter = patch_concentrations[:, :, 8:14]
    sumOuter = np.sum(lOuter) / (patch_concentrations.shape[0]*patch_concentrations.shape[1])
    sumInner = np.sum(lInner) / (patch_concentrations.shape[0]*patch_concentrations.shape[1])
    asym = np.rint(1600 * (1 - (sumInner / sumOuter))).astype(int)

    # @TODO change this to based on 64 lipids per subpatch (to explore rounding)
    # @TODO depdning on variations in cons from PF maybe to "proper" pre rounding here ???
    #old_lOuter = np.rint(64 * lOuter / sumOuter).astype(int)
    #old_lInner = np.rint(64 * lInner / sumInner).astype(int)

    # Compute a probablistic construction of the particles in a cell.
    # We use a cumlative distribution function to
    # TODO: We should really assert or something here to verify that the
    # inner and outer leaflets are the same size.
    # Iterate through each grid point and convert to the number of particles
    # based on the probability.
    for i in range(0, lInner.shape[0]):
        for j in range(0, lInner.shape[1]):
            #print(lInner[i, j])
            lInner[i, j] = generate_probabilistic_patch(lInner[i, j], 64)
    lInner = lInner.astype(int)
    for i in range(0, lOuter.shape[0]):
        for j in range(0, lOuter.shape[1]):
            lOuter[i, j] = generate_probabilistic_patch(lOuter[i, j], 64)
    lOuter = lOuter.astype(int)

    # @TODO Use this to explor lipid consetnrations - excel magic ( a 0 values will be left out of the list but 1/1000 is ok)
    # not needed for production and makes other consentrations a little strange
    # lOuter[lOuter == 0] = 1
    # lInner[lInner == 0] = 1

    # Convet into the right shape (adding empty lipids as 0 and chansing 5x5 to 25)
    # lipidUpperRatioArray = [[243,   0, 121, 20, 61, 242, 0, 0, 313]]
    # lipidLowerRatioArray = [[0,   139, 75, 54, 161, 108, 161, 22, 280]]
    lipidUpperRatioArray = np.insert(np.insert(np.insert(lOuter, 1, 0, axis=2), 6, 0, axis=2), 7, 0, axis=2)
    lipidLowerRatioArray = np.insert(lInner, 0, 0, axis=2)
    lipidUpperRatioArray = np.reshape(lipidUpperRatioArray, (patch_concentrations.shape[0]*patch_concentrations.shape[1], len(lipidTypes)))
    lipidLowerRatioArray = np.reshape(lipidLowerRatioArray, (patch_concentrations.shape[0]*patch_concentrations.shape[1], len(lipidTypes)))

    half_box_dim = 31.0 / 2 # This is fixed here and in particlesim.py
    #krasPosArray = p.rasPositions - p.rasPositions[0]  # Center at 0,0
    krasPosArray = patch_rpositions # KRAS in patches are now centerd at the first KRAS
    krasPosArray = krasPosArray + [half_box_dim, half_box_dim]

    if (krasPosArray.max() > (half_box_dim * 2) or krasPosArray.min() < 0):
        error = "---> KRAS posisions are outside the CG box size; pos are: " + np.array2string(krasPosArray)
        LOGGER.error(error)
        raise ValueError(error)

    # Print CG patch input
    LOGGER.info("CG_patch input:")
    LOGGER.info("macro patch summ_outer = {}".format(sumOuter))
    LOGGER.info("macro patch summ_inner = {}".format(sumInner))
    LOGGER.info("patch asym             = {}".format(asym))
    LOGGER.info("Lipids outer cons.     = {}".format(str(lipidUpperRatioArray)))
    LOGGER.info("Lipids inner cons.     = {}".format(str(lipidLowerRatioArray)))
    LOGGER.info("RAS init pos.          = {}".format(krasPosArray))
    LOGGER.info("RAS states             = {}".format(patch_rstates))

    # Select KRAS structures acording to states - pfPatch state to saved files + random value (e.g. 1-1000)
    # "KRAS-s"+str(cS)+"-r"+str(randint(1, librarySize))+".gro"
    librarySize = 1000  # numer of diffrent configurations saved for each state
    krasStructs = []
    for i in range(len(patch_rstates)):
        cInt = randint(1, librarySize)
        cState = int(patch_rstates[i])
        if (cState != 1 and cState != 2):
            error = "---> KRAS state {} is not supported: ".format(cState)
            LOGGER.error(error)
            raise ValueError(error)
        krasStruct = "{}/random_state{}.{}.KRAS.pbc.mol.gro".format(particlesim.dirKRASStructures, cState, cInt)
        LOGGER.info("Get KRAS structure file: KRAS # {}, in state {} and rand {}, name {}".format(i, cState, cInt, krasStruct))
        krasStructs.append(krasStruct)


    patchSim = particlesim.ParticleSim(lipidTypes, lipidUpperRatioArray, lipidLowerRatioArray, [patch_concentrations.shape[0], patch_concentrations.shape[1]],
                                       krasStructs, krasPosArray, asym=-asym, simName=ptachID, lipidFullAsignment=False)
    return patchSim


def read_patches_from_combined(ppath, patch_ids):
    """Function from Harsh - https://code.ornl.gov/fdinatal/campaign1s_inputs/-/blob/bhatia/scripts_for_macro/read_patches.py."""
    assert isinstance(ppath, str)
    assert isinstance(patch_ids, (np.ndarray,list))

    fpattern = 'patches_{}.npz'

    # TODO: group all the patches of the same frame so
    # no need to read the file again
    patch_concs = []
    patch_rstates = []
    patch_rpositions = []
    for p in patch_ids:
        fid = p//300
        pid = p%300

        d = np.load(os.path.join(ppath, fpattern.format(fid)))
        patch_concs.append(d['patches'][pid])
        patch_rstates.append(d['ras_states'][pid])
        patch_rpositions.append(d['ras_positions'][pid])

    return np.array(patch_concs), np.array(patch_rstates), np.array(patch_rpositions)


# modified from https://github.com/LLNL/maestrowf/blob/develop/maestrowf/maestro.py by Francesco Di Natale, dinatale3@llnl.gov.
def setup_argparser():
    """Set up the program's argument parser."""
    parser = ArgumentParser(prog="CreateSims", description="CreateSims: creates and equilibrates a Martini particle based simulation \n"
                            "from on a macro model converted pf_patch (a few other creation methods are also supported see below). ",
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument('--fstype', type=str, default="mummi", choices=mummi_ras.get_interfaces(),
                        help='Chose how to write files/data (default: mummi)')

    parser.add_argument("-o", "--outpath", type=str, default="",
                        help="Output path, shoudl be a full path and if not set use env $MUMMI_ROOT/sims-cg/[sim name] (default in MuMMI)")

    parser.add_argument("-ol", "--outlocal", type=str, default="",
                        help="Local output path, shoudl be a full path and if set run everything in this path and only copy to --outpath in the end \n"
                        "WARNING, if --outpath not set (default in MuMMI) /sims-cg/[sim name]_[timestamp] will be added to this path.")

    parser.add_argument("-path", "--inpath", type=str, required=True,
                       help="Input path for patch data (db or fs path)")

    parser.add_argument("-d", "--loglevel", type=int, default=2,
                        help="Level of logging messages to be output:\n"
                             "5 - Critical\n"
                             "4 - Error\n"
                             "3 - Warning\n"
                             "2 - Info (Default)\n"
                             "1 - Debug")

    parser.add_argument("-l", "--logpath", type=str, default="./",
                        help="Path for log information")

    parser.add_argument("-c", "--logstdout", action="store_true",
                        help="Log to stdout in addition to a file")

    parser.add_argument("--logfile", type=str, default="createsims",
                        help="Name of createsims logger")

    parser.add_argument("-p", "--patch", type=str, required=True,
                        help="Depending of prefix pfpatch, load or create will:\n"
                             "pfpatch: Name of patch to crate sim for (db key or fs pickle file name without the ending).\n"
                             "load: Load in CG patch psim.npz in path 'laod <path>' \n"
                             "create: Create a CG patch based on provided input")

    parser.add_argument("-g", "--gromacs", type=str, default="/usr/local/gromacs/bin/",
                        help="GROMACS bin directory, used for e.g. grompp, mdrun minimization, make_ndx")

    parser.add_argument("-m", "--mpi", type=str, default="",
                        help="mpi run command used for all GROMACS MD mdrun's")

    parser.add_argument("-r", "--mdrunopt", type=str, default="",
                        help="Options to add to GROMACS mpi MD mdrun's")

    parser.add_argument("--simnum", type=str, default="",
                        help="Only used in offline runs and testing adds to the end of sim name in some use cases")

    parser.add_argument("--mini", action='store_true',
                        help="If set we will use mini-MuMMI setting")

    return parser


def store_files(outpath: str, sim_path: str, failure: bool):
    if sim_path != outpath and os.path.exists(sim_path):
        dest_path = f'{outpath}/{utils.add_timestamp("failed")}' if failure else outpath
        LOGGER.debug("---> Copy local %s to global folder %s", sim_path, dest_path)
        os.makedirs(dest_path, exist_ok=True)
        utils.sys_call("cp -r {}/* {}/".format(sim_path, dest_path))


def cleanup_files(path: str):
    if path and os.path.exists(path):
        #return # @TODO put this back in commented out for local debug
        LOGGER.debug(f"Cleaning up {path}")
        utils.sys_call(f"rm -r {path}")


def sigterm_wrapper(outpath: str, local_path: str):

    def sigterm_handler(signum, frame):
        LOGGER.info(f"Caught {signum} -- Aborting.")
        store_files(outpath, local_path, True)
        cleanup_files(local_path)
        # terminate all the way so not exicute finally block
        os._exit(1)

    return sigterm_handler


def get_npz_as_gro(coords, fname, gro_ref_file):
    """
    Return the coordinates from numpy arrays to gro file on disk
    """
    if "mda" not in dir():
        import MDAnalysis as mda

    # load reference
    protein = mda.Universe(gro_ref_file)  # read in reference gro file
    protein.atoms.positions = coords
    protein.atoms.write(fname)

    return protein

def get_protein_contacts_from_z(ppos, atomids_dict, cutoff, farnesyl_id):
    """Return an array with protein counts"""

    domains = atomids_dict.keys()

    ncontacts = {}
    ref_z = ppos[farnesyl_id,2]  # get farnesly height 

    for dom in domains:
        bead_zs = ppos[atomids_dict[dom].indices,2]
        diff = np.abs(bead_zs - ref_z)
        bary = diff <= cutoff
        ncontacts[dom] = bary.sum()

    return ncontacts

def get_state_from_protein_coordinates(ppos,
                                       atomids_dict,
                                       codes_to_state,
                                       cutoff,
                                       contacts_threshold,
                                       farnesyl_id,
                                       COM_1,
                                       COM_2):
    """
    Return state code for this protein structure.

    checks if we are in state A or B

    cutoff is the z distance to the CYF BB required and contacts_threshold is the
    number of beads within cutoff required for the domain to be "touching"
    """
    M1_domains = ["KD1", "CRD1", "loop11", "loop21", "b1433", "RAS1"]
    M2_domains = ["KD2", "CRD2", "loop12", "loop22", "RAS2"]

    _M1_dict = {key : atomids_dict[key] for key in M1_domains}
    _M2_dict = {key : atomids_dict[key] for key in M2_domains}

    # do M1 contacts
    contacts = get_protein_contacts_from_z(ppos, _M1_dict, cutoff, farnesyl_id[0])

    # make an array of all contacts
    code1 = ""
    # now accumulate all the codes into the dictionaries  
    for dom in M1_domains:

        # we are always assigning RAS to 1
        if dom == "RAS1":
            code1 += "1"
            continue

        if contacts[dom] >= contacts_threshold:
            code1 += "1"
        else:
            code1 += "0"

    # tag this with KD orientation 
    kdpos = [atomids_dict['KD1'].positions[:,0].mean(),
             atomids_dict['KD1'].positions[:,1].mean()]
    crdpos = [atomids_dict['CRD1'].positions[:,0].mean(),
              atomids_dict['CRD1'].positions[:,1].mean()]

    kdvec = np.array([kdpos[0] - COM_1[0], kdpos[1] - COM_1[1]])
    crdang = np.arctan2(crdpos[0] - COM_1[0], crdpos[1] - COM_1[1])

    kdrot = np.array([kdvec[0] * np.cos(-crdang) - kdvec[1] * np.sin(-crdang),
                      kdvec[0] * np.sin(-crdang) + kdvec[1] * np.cos(-crdang)]) 

    ymean = np.arctan2(kdrot[0], kdrot[1])

    if ymean >= 0:
        ytag = "ypos"
    else:
        ytag = "yneg"

    if code1 in codes_to_state.keys():
        state = code1
    else:
        state = "-".join([code1,ytag])

    M1_state = codes_to_state[state]

    # do M2 contacts
    contacts = get_protein_contacts_from_z(ppos, _M2_dict, cutoff, farnesyl_id[1])

    code2 = ""
    # now accumulate all the codes into the dictionaries  
    for dom in M2_domains:

        # we are always assigning RAS to 1
        if dom == "RAS2":
            code2 += "1"
            continue

        if contacts[dom] >= contacts_threshold:
            code2 += "1"
        else:
            code2 += "0"

    # tag this with KD orientation 
    kdpos = [atomids_dict['KD2'].positions[:,0].mean(),
             atomids_dict['KD2'].positions[:,1].mean()]
    crdpos = [atomids_dict['CRD2'].positions[:,0].mean(),
              atomids_dict['CRD2'].positions[:,1].mean()]

    kdvec = np.array([kdpos[0] - COM_2[0], kdpos[1] - COM_2[1]])
    crdang = np.arctan2(crdpos[0] - COM_2[0], crdpos[1] - COM_2[1])

    kdrot = np.array([kdvec[0] * np.cos(-crdang) - kdvec[1] * np.sin(-crdang),
                      kdvec[0] * np.sin(-crdang) + kdvec[1] * np.cos(-crdang)]) 

    ymean = np.arctan2(kdrot[0], kdrot[1])

    if ymean >= 0:
        ytag = "ypos"
    else:
        ytag = "yneg"

    if code2 in codes_to_state.keys():
        state = code2
    else:
        state = "-".join([code2,ytag])

    M2_state = codes_to_state[state]

    return ("_".join(["state",M1_state]), "_".join(["state", M2_state]))

def run_macro_code(prot_state, tmpl_dir, sim_dir, influx=False):
    """
    Setup and launch a macro simulaiton using a fixed RAS state.
    tmlp_dir -> directory in mummi_resources that contains the template files for a macro run
    sim_dir -> directory in localpath in which to run macro sim 
    """

    # call subprocess to launch macro
    # TODO: replace here with flux mini-like call... not using srun / jsrun
    #utils.sys_call("pushd {}; srun -n4 {}/gridsim2dras run-2Dforce.cfg >> macro_createsims.out  2>&1; popd".format(sim_dir,
    #                                                                              tmpl_dir))
    #utils.sys_call("pushd {}; mpirun -np 2 ./gridsim2dras run-2Dforce.cfg >> macro_createsims.out  2>&1; popd".format(sim_dir))
    #utils.sys_call("flux mini run -n 4 -o mpi=spectrum gridsim2dras run-2Dforce.cfg >> macro_createsims.out  2>&1 ", run_dir=sim_dir)
    #utils.sys_call('export OMPI_MCA_btl="self"; gridsim2dras run-2Dforce.cfg >> macro_createsims.out  2>&1 ', run_dir=sim_dir)
    # utils.sys_call('env -i PATH="$PATH"  OMPI_MCA_btl="self" OMPI_COMM_WORLD_RANK=0 OMPI_COMM_WORLD_SIZE=1 gridsim2dras run-2Dforce.cfg >> macro_createsims.out  2>&1', run_dir=sim_dir)
    if influx:
        # was -np 1, tried -np 8
        # utils.sys_call('OMPI_MCA_schizo=^flux OMPI_MCA_pmix=^flux mpirun -np 1 gridsim2dras run-gc.cfg >> macro_createsims.out  2>&1', run_dir=sim_dir)
        #utils.sys_call('flux run -n1 -c1 gridsim2dras run-gc.cfg >> macro_createsims.out  2>&1', test=True, run_dir=sim_dir) # in current dir
        utils.sys_call('flux run -n1 -c1 gridsim2dras run-gc.cfg >> macro_createsims.out  2>&1', run_dir=sim_dir)
        # utils.sys_call2('OMPI_MCA_schizo=^flux OMPI_MCA_pmix=^flux gridsim2dras run-gc.cfg >> macro_createsims.out  2>&1', run_dir=sim_dir)
        # mpi=none
    else:
        #utils.sys_call('/lustre/orion/lrn005/proj-shared/spack/bin/gridsim2dras run-gc.cfg >> macro_createsims.out  2>&1', run_dir=sim_dir) # in current dir
        utils.sys_call('gridsim2dras run-gc.cfg >> macro_createsims.out  2>&1', test=True, run_dir=sim_dir) # in current dir
        utils.sys_call('gridsim2dras run-gc.cfg >> macro_createsims.out  2>&1', run_dir=sim_dir) # in current dir

    #utils.sys_call('./gridsim2dras run-gc.cfg >> macro_createsims.out  2>&1', run_dir=sim_dir) # in current dir 
    #utils.sys_call('/usr/gapps/mummi/lassen/macro-gc/c1-iso/gridsim2dras run-gc.cfg >> macro_createsims.out  2>&1', run_dir=sim_dir) # full path 

    # confirm macro ran and we have output, if not return raise error
    # use LOGGER.error here.  
    assert os.path.exists("{}/snapshot.00000001".format(sim_dir)), "Macro code failed to generate outout."

    return

def setup_macro_sim_dir(tldir, sdir, prot_state):
    """
    Setup macro simulation.

    prot_state is now a dictionary of parameters:
    """

    # copy cfg, atom and bonds files to sim dir
    #utils.sys_call(f"cp {tldir}/bonds.init {sdir}")
    utils.sys_call(f"cp {tldir}/RASRAF-* {sdir}")
    utils.sys_call(f"cp {tldir}/run-gc.cfg {sdir}")
    utils.sys_call(f"cp {tldir}/atoms.init {sdir}")
    utils.sys_call(f"cp {tldir}/bonds.init {sdir}")
    utils.sys_call(f"cp {tldir}/cofr-14lips-martini3-0.dat {sdir}")
    utils.sys_call(f"cp {tldir}/rep12_pot.dat {sdir}")

    #utils.sys_call(f"cp {tldir}/gridsim2dras {sdir}")

    # replace RSEED
    for i in range(4):
        utils.sys_call('sed -zi "s/_RSEED/{}/" {}/run-gc.cfg'.format(hex(np.random.randint(2**60)), sdir))

    # replace M1 parameters
    for i in range(2):
        utils.sys_call('sed -zi "s/_M1/{}/" {}/atoms.init'.format(prot_state['M1'], sdir))
    utils.sys_call('sed -zi "s/_X1/{}/" {}/atoms.init'.format(prot_state['X1'], sdir))
    utils.sys_call('sed -zi "s/_Y1/{}/" {}/atoms.init'.format(prot_state['Y1'], sdir))
    utils.sys_call('sed -zi "s/_A1/{}/" {}/atoms.init'.format(prot_state['A1'], sdir))

    # replace M2
    for i in range(2):
        utils.sys_call('sed -zi "s/_M2/{}/" {}/atoms.init'.format(prot_state['M2'], sdir))
    utils.sys_call('sed -zi "s/_X2/{}/" {}/atoms.init'.format(prot_state['X2'], sdir))
    utils.sys_call('sed -zi "s/_Y2/{}/" {}/atoms.init'.format(prot_state['Y2'], sdir))
    utils.sys_call('sed -zi "s/_A2/{}/" {}/atoms.init'.format(prot_state['A2'], sdir))

    return


def main():
    """
    Launcher for testing and setting up simulations.

    Can be launced from the workflow or standalone.
    Example for a standalone version run:
        export $MUMMI_CODE=/g/g90/helgi/mummi
        source $MUMMI_CODE/setup/setup.env.sh
        python $MUMMI_CODE/mummi/transformations/createsim/createsims.py \
          --gromacs gmx --loglevel 1 --fstype mummi \
          --logpath ./ -c --mpi "gmx mdrun" \
          --patch pfpatch_000000384620 --mdrunopt " -nt 96 -rdd 2.0 -ntomp 4 -dd 4 3 2" \
          --inpath '/p/gpfs1/bhatia4/pilot2/campaign1star/patches/' --outpath ./ --simnum test_run_1
    """

    # Set up the necessary base data structures to begin study set up.
    parser = setup_argparser()
    args = parser.parse_args()
    print (f'Launching Createsim from {os.getcwd()}: {args}')

    global gromacs
    global mdrun_o
    global mdrun_mpi

    gpu_id = os.environ.get("CUDA_VISIBLE_DEVICES", "UNSET")
    binding_to_gpu = ""
    if gpu_id != "UNSET":
        gpu_id = gpu_id.split(',')
        if len(gpu_id) > 1:
            LOGGER.warning(f"Detected more than 1 GPUs per GROMACS {gpu_id}")
            gpu_id = '-'.join(gpu_id)
        else:
            gpu_id = int(gpu_id[0])
        # On Frontier there is one L3 cache per GPU
        binding_to_gpu = f"hwloc-bind $(hwloc-calc L3:{gpu_id}) -- "

    gromacs = f"{binding_to_gpu} {args.gromacs}"
    mdrun_o = args.mdrunopt
    mdrun_mpi = f"{binding_to_gpu} {args.mpi}"

    # If no mpi mdrun given use default gromacs mdrun
    if (len(mdrun_mpi) < 1):
        mdrun_mpi = gromacs + " mdrun"

    # Sets up MuMMI and the logger
    mummi_core.init()
    mummi_core.create_root()
    mummi_core.init_logger(argparser=args)

    # Get args
    fstype = args.fstype
    inpath = args.inpath
    patchName = args.patch
    outpath = args.outpath  # this is gpfs path
    localpath = args.outlocal  # this is localdisk path
    mini_mummi = args.mini # True if we use Mini MuMMI

    # Show options
    LOGGER.info("Running createsims, options are:")
    LOGGER.info(args)
    LOGGER.info("GROMACS       = {}".format(gromacs))
    LOGGER.info("GROMACS single precision path " + str(shutil.which(args.gromacs)))
    LOGGER.info("mdrun_o       = {}".format(mdrun_o))
    LOGGER.info("mdrun_mpi     = {}".format(mdrun_mpi))
    LOGGER.info("init outpath  = {}".format(outpath))
    LOGGER.info("init outlocal = {}".format(localpath))

    # This is the working dir for the setup
    simpath = outpath if not localpath else localpath
    # Set up the SIGTERM handler
    signal.signal(signal.SIGTERM, sigterm_wrapper(outpath, simpath))
    signal.signal(signal.SIGINT, sigterm_wrapper(outpath, simpath))
    signal.signal(signal.SIGUSR1, sigterm_wrapper(outpath, simpath))
    signal.signal(signal.SIGUSR2, sigterm_wrapper(outpath, simpath))

    # @TODO This statement is to kill existing GROMACS runs on the node. This is
    # a temporary fix for flux overscheduling to some ranks.
    # @TODO removed as we definitly can't kill gmx as many creatsims or cganalysis running on the same note
    #utils.sys_call('pkill gmx_mpi*', return_codes=[0,1])

    # Initilize ouput DBR/files etc
    iointerface = mummi_ras.get_io(fstype)
    LOGGER.info('Using IO_Interface [{}] for all output'.format(iointerface))

    # Flags used to set if creatims is successful or failed with error and should not be resterted
    flag_success, flag_failure = Naming.status_flags('createsim')
    marked_failure = True

    # Time creatsims
    cs_begin = time.time()

    # Create or get particlesim and macro model patch
    patchSim = []
    p = None
    try:
        if (len(patchName) == 0):
            error = "createsims was called without -p <patchName>"
            LOGGER.error(error)
            raise ValueError(error)
        elif (patchName[:6] == "gc-sim" or patchName[:6] == "oneRas" or patchName[:2] == "mu"):
            # This is the C4 method used
            LOGGER.info("Get GC macro patch {} and convert it to a CG patch".format(patchName))
            p = fetch_gc_patches(inpath, patchName)
            patchSim = convert_gc_patch_to_cg_patch(p, args.simnum)
        elif (patchName[:7] == "pfpatch"):
            # This is the C3 method used
            LOGGER.info("Get macro patch {} and convert it to a CG patch".format(patchName))
            p = iointerface.load_patches(inpath, [patchName])[0]
            patchSim = convert_pf_patch_to_cg_patch(p, args.simnum)
        elif (patchName[:4] == "load"):
            # Loads a particlesim object that has been created and saved
            LOGGER.info("Load CG patch {}".format(patchName))
            patchSim = particlesim.load_sim(patchName[5:].strip())
        elif (patchName[:6] == "create"):
            # This method is for creating specific types of simulations - both hardcoded and from comand line
            LOGGER.info("Create CG patch based on {} with name {}".format(patchName, args.simnum))
            patchSim = creat_cg_patch(patchName[6:].strip(), args.simnum, mini_mummi)
        elif (patchName[:8] == "test_raw"):
            LOGGER.info("Testing raw patch '%s'", patchName)
            ppath = "/p/gpfs1/bhatia4/pilot2/campaign1star/patches_37"
            pids = np.array([1062730])
            patch_concentrations, patch_rstates, patch_rpositions = read_patches_from_combined(ppath, pids)
            #patchSim = convert_macro_to_cg_patch("rawpatch_1062730", patch_concentrations[0], patch_rstates[0], patch_rpositions[0])
            patchSim = convert_macro_to_cg_patch("rawpatch_1062730", patch_concentrations[0], patch_rstates, patch_rpositions)

        elif (patchName[0:9] == "structure_pc3"):    # This is what was used for "structure_pc3" in Pilot 2 pre-campaing3

            # -------------------------------------
            # files and directories needed
            # -------------------------------------

            # temporary directory to run macro code. 
            macro_sim_dir = os.path.join(localpath, "macro_sim_dir")
            # location of macro template
            #tmpl_dir = os.path.join(Naming.dir_res('macro'), "createsim_template")
            os.path.join(Naming.dir_res('macro_c2'))
            # location to write gro file from npz input
            gro_file = os.path.join(localpath,"struct_createsim.gro")
            # ml absolute input path
            inPathInteration=os.path.join(inpath, patchName.split('_')[1])

            # look up a reference gro file from the library 
            gro_ref = os.path.join(Naming.dir_res('structures'),
                                   "RAS-RAF",
                                   "RAS-RAF_state_ma_0.gro")

            # should make an isolated sub dir in localpath in which to write gro and run macro sim
            os.makedirs(macro_sim_dir, exist_ok=True)

            struct = iointerface.load_npz(inPathInteration, [patchName])[0]["data"]

            # convert npz to gro and write to disk
            if "mda" not in dir():
                import MDAnalysis as mda

            # load reference
            _prot = mda.Universe(gro_ref)  # read in reference gro file
            _prot.atoms.positions = struct
            _prot.atoms.write("struct_createsim_init.gro")

            # instatiate an MDanalysis universe for analysing gro etc
            #kras_indices = mummi_ras.online.cg.get_tiltrot_z_state_multiple.get_protein_info(_prot,'resname ACE1 and name BB')
            _kr_indx = _prot.select_atoms('resname ACE1 and name BB')
            kras_indices = [(k.resid, 'RAS-RAF') for k in _kr_indx]
            statesList = mummi_ras.online.cg.get_tiltrot_z_state_multiple.getKRASstates(_prot, kras_indices)

            # hardcode state here, will eventually need routine to determine this from gro file:
            # import this library from RAS states? --> could test
            # coded: 1 -> ma
            #        2 -> mb
            #        3 -> za

            # hardcode this to RASRAF 
            _statenames = {1 : "mRASa_mRAFa",
                           2 : "mRASb_mRAFb",
                           3 : "zRASa_zRAFa"}

            _state_index = statesList[0][0][-1]
            state = _statenames[_state_index]

            # set up macro simulation dir
            # run this on outlocal where Mummi holds temp local files
            # should make an isolated sub dir in localpath. 
            setup_macro_sim_dir(tmpl_dir, macro_sim_dir, state)

            # launch macro job
            run_macro_code(state, tmpl_dir, macro_sim_dir, influx=True)
            # breakpoint()

            # pulls out the last state from from the above run
            p = fetch_gc_patches(localpath, "macro_sim_dir_000000000001", "simple")

            # sets coordinates in patchSim 
            patchSim = convert_gc_patch_to_cg_patch(p, args.simnum, localpath)

            # Pre-Campaign 3 specific - Fix rotation of protein 
            _protBB = _prot.select_atoms('name BB') # only one protein in pre-campaign 3
            RESIDXS_GDOMAIN = np.arange(0, 165, dtype=int)
            RESIDXS_CRD_4B = np.array([276,277,278,279, 289,290,291,292])
            RESIDXS_FARNESYL_4B = np.array([183])

            resids = ' '.join([str(i) for i in RESIDXS_GDOMAIN])
            com_GDOMAIN = _protBB.select_atoms(f'resid {resids}').center_of_mass()
            resids = ' '.join([str(i) for i in RESIDXS_CRD_4B])
            com_CRD_4B = _protBB.select_atoms(f'resid {resids}').center_of_mass()
            resids = ' '.join([str(i) for i in RESIDXS_FARNESYL_4B])
            com_FARNESYL_4B = _protBB.select_atoms(f'resid {resids}').center_of_mass()
            LOGGER.info("Pre-campaign3 - set rotation: gdomain = {} A, crd_4B = {} A, farnesyl = {} A".format(com_GDOMAIN, com_CRD_4B, com_FARNESYL_4B))

            dxy = np.array([com_GDOMAIN[0], com_GDOMAIN[1], 0])
            dvec = com_CRD_4B - com_GDOMAIN
            angle = np.rad2deg(np.arctan2(dvec[1], dvec[0]))
            LOGGER.info("Pre-campaign3 - rotate by angle  = {}".format(angle))

            # rotate
            _prot.atoms.translate(-dxy)
            _prot.atoms.rotateby(-angle, [0.,0.,1.], point=[0,0,0])
            _prot.atoms.translate(dxy)

            # translate z-axes to correct location (somewhat below the membrane)
            #_prot.atoms.translate([0, 0, com_FARNESYL_4B[2] - (10*5.455)]) - not needed was old 

            # Pre-Campaign 3 specific - Change protein location from RAS G-domain to farnesyl (Note, MDAnalysis is in A, macro model in nm) 
            patchSim.proteinUnits[0]["position"] = [patchSim.proteinUnits[0]["position"][0] - (com_GDOMAIN[0] - com_FARNESYL_4B[0]) / 10, patchSim.proteinUnits[0]["position"][1] - (com_GDOMAIN[1] - com_FARNESYL_4B[1]) / 10]

            # back to normal stuff (before pre-campaign 3) 
            _prot.atoms.write(gro_file)
            patchSim.proteinUnits[0]["structure"] = gro_file
            # enforce that it's always of type "RAS-RAF"
            patchSim.proteinUnits[0]["type"] = "RAS-RAF"

        elif (patchName[0:9] == "structure_c2_mini"):
            # this is for C2 and mini_mummi - change back to e.g. "structure" and modify C3 below

            # convert npz to gro and write to disk
            if "mda" not in dir():
                import MDAnalysis as mda

            # location to write gro file from npz input
            gro_file = os.path.join(localpath,"struct_createsim.gro")
            # ml absolute input path
            inPathInteration = os.path.join(inpath, patchName.split('_')[1])
            gro_input = inPathInteration + "/" + patchName + ".gro"

            if not os.path.exists(gro_input):
                LOGGER.warning(f"File {gro_input} does not exist.")
                gro_input = patchName + ".gro"
                if not os.path.exists(gro_input):
                    raise ValueError(f"File {gro_input} does not exist. Aborting")


            LOGGER.info(f"Read in generated input structure after healing = {gro_input}")
            _prot = mda.Universe(gro_input)

            if mini_mummi:
                LOGGER.info(f"Setting up patch {patchName} in mini-MuMMI style")
                patchSim = creat_cg_patch("", "RAS_RAF", mini_mummi)
                # Update FF parameters
                added_mdp_options = " -DACYM280"
                patchSim.proteinUnits[0]["ff_options"] = added_mdp_options # note these are used in .mdp files to change the protein FF parameters
                # Change protein name 
                patchSim.proteinUnits[0]["structure"] = f"{patchName}.gro" # assuming file is copied into sim dir
                _prot.atoms.write(gro_file)
                patchSim.proteinUnits[0]["structure"] = gro_file
            else:
                # -------------------------------------
                # files and directories needed
                # -------------------------------------

                # temporary directory to run macro code. 
                # location of macro template
                tmpl_dir = Naming.dir_res('macro_c2')
                # look up a reference gro file from the library 
                gro_ref = os.path.join(Naming.dir_res('ref_structures'), "ras-raf-14-3-3-ref-CG.gro")
                macro_sim_dir = os.path.join(localpath, "macro_sim_dir")

                # should make an isolated sub dir in localpath in which to write gro and run macro sim
                if not os.path.exists(macro_sim_dir):
                    LOGGER.info(f"macro_sim_dir not found. Creating a new directory on {localpath}")
                    utils.sys_call(f"mkdir {macro_sim_dir}")

                # load reference
                #_prot = mda.Universe(gro_ref)  # read in reference gro file
                #_prot.atoms.positions = struct 
                #_prot.atoms.write("struct_createsim_init.gro")

                # load generated healead .gro file directly 


                # The stucte posistion are from the .npz - ML generated before healing run so don't use
                #struct = iointerface.load_npz(inPathInteration, [patchName])[0]["data"]
                #_prot.atoms.positions = struct 

                # Check structure for needed FF changes
                added_mdp_options = ""
                # zn_cutoff returns True if it is close enough to coordiante (B state)
                if zn_cutoff(_prot, 8, patchName):
                    added_mdp_options += " -DACYM280"
                # 365 cutoff returns True if it is far enough away to not be bound (B state)
                if SEP365_cutoff(_prot, 8, patchName):
                    added_mdp_options += " -DASEP365"

                inBstate = zn_cutoff(_prot, 8, patchName) and SEP365_cutoff(_prot, 8, patchName)

                # canonical ordering of domains for B state
                domains_dict = {"KD": _prot.select_atoms('index 1136 to 1768'),
                                "CRD": _prot.select_atoms('index 637 to 740'),
                                "loop1": _prot.select_atoms('index 753 to 928'),
                                "loop2": _prot.select_atoms('index 955 to 1135'),
                                "b1433": _prot.select_atoms('index 2358 to 2896'),
                                "RAS": _prot.select_atoms('index 0 to 437')}

                """
                Ordering for domainsi n A state.
                domains_dict = {"KD" : _prot.select_atoms('index 1136 to 1768'),
                                "loop2" : _prot.select_atoms('index 955 to 1135'),
                                "b1433" : _prot.select_atoms('index 2358 to 2896'),
                                "RAS" : _prot.select_atoms('index 0 to 437')}
                """


                # this is for B state classification only
                codes_to_state = {"000001": "a",
                                "000011": "a",
                                "000101": "a",
                                "000111": "a",
                                "010001": "b",
                                "010011": "b",
                                "010101": "b",
                                "010111": "b",
                                "100001-ypos": "c",
                                "100011-ypos": "c",
                                "100101-ypos": "c",
                                "100111-ypos": "c",
                                "110001-ypos": "d",
                                "110011-ypos": "d",
                                "110101-ypos": "d",
                                "110111-ypos": "d",
                                "100001-yneg": "e",
                                "100011-yneg": "e",
                                "100101-yneg": "e",
                                "100111-yneg": "e",
                                "110001-yneg": "f",
                                "110011-yneg": "f",
                                "110101-yneg": "f",
                                "110111-yneg": "f",
                                "001001": "g",
                                "001011": "g",
                                "001101": "g",
                                "001111": "g",
                                "011001": "b",
                                "011011": "b",
                                "011101": "b",
                                "011111": "b",
                                "101001-ypos": "i",
                                "101011-ypos": "i",
                                "101101-ypos": "i",
                                "101111-ypos": "i",
                                "111001-ypos": "d",
                                "111011-ypos": "d",
                                "111101-ypos": "d",
                                "111111-ypos": "d",
                                "101001-yneg": "k",
                                "101011-yneg": "k",
                                "101101-yneg": "k",
                                "101111-yneg": "k",
                                "111001-yneg": "f",
                                "111011-yneg": "f",
                                "111101-yneg": "f",
                                "111111-yneg": "f"}

                state = get_state_from_protein_coordinates(_prot.atoms.positions,
                                                           domains_dict,
                                                           codes_to_state,
                                                           15,
                                                           5,
                                                           430)  # farnesyl id is hardcoded as BB

                setup_macro_sim_dir(tmpl_dir, macro_sim_dir, state)

                # launch macro job
                run_macro_code(state, tmpl_dir, macro_sim_dir, influx=True) # True for Frontier right now (April 2025)

                # pulls out the last state from from the above run
                p = fetch_gc_patches(localpath, "macro_sim_dir_000000000001", "simple")

                # sets coordinates in patchSim 
                patchSim = convert_gc_patch_to_cg_patch(p, args.simnum, localpath)

                # Pre-Campaign 3 specific - Fix rotation of protein 
                _protBB = _prot.select_atoms('name BB') # only one protein in pre-campaign 3
                # Martini 2 - pre-campaign 3 numbers 
                '''
                RESIDXS_GDOMAIN = np.arange(0, 165, dtype=int)
                RESIDXS_CRD_4B = np.array([276,277,278,279, 289,290,291,292])
                RESIDXS_FARNESYL_4B = np.array([183])
                resids = ' '.join([str(i) for i in RESIDXS_GDOMAIN])
                com_GDOMAIN = _protBB.select_atoms(f'resid {resids}').center_of_mass()
                resids = ' '.join([str(i) for i in RESIDXS_CRD_4B])
                com_CRD_4B = _protBB.select_atoms(f'resid {resids}').center_of_mass()
                resids = ' '.join([str(i) for i in RESIDXS_FARNESYL_4B])
                com_FARNESYL_4B = _protBB.select_atoms(f'resid {resids}').center_of_mass()
                '''

                # Martini 3 - campaign 3 indexes
                com_GDOMAIN = _protBB.select_atoms('index 0 to 380').center_of_mass() # should be KRAS G-domain 
                com_CRD_4B = _protBB.select_atoms('index 649 to 660 or index 685 to 694').center_of_mass() # this is resid 241 to 244 (THR, PHE, PHE, THR) and resid 254 to 257 (LUE LEU PHE GLN) - CRD membrane interacting part
                com_FARNESYL_4B = _protBB.select_atoms('index 430 to 431').center_of_mass() # should be BB of CYF resnum 185 in KRAS 
                LOGGER.info("M3 - set rotation: gdomain = {} A, crd_4B = {} A, farnesyl = {} A".format(com_GDOMAIN, com_CRD_4B, com_FARNESYL_4B))

                dxy = np.array([com_GDOMAIN[0], com_GDOMAIN[1], 0])
                dvec = com_CRD_4B - com_GDOMAIN
                angle = np.rad2deg(np.arctan2(dvec[1], dvec[0]))
                LOGGER.info("M3 - rotate by angle  = {}".format(angle))

                # rotate
                _prot.atoms.translate(-dxy)
                _prot.atoms.rotateby(-angle, [0.,0.,1.], point=[0,0,0])
                _prot.atoms.translate(dxy)

                # M3 specific - Change protein location from RAS G-domain to farnesyl (Note, MDAnalysis is in A, macro model in nm) 
                newx = patchSim.proteinUnits[0]["position"][0] - (com_GDOMAIN[0] - com_FARNESYL_4B[0])/10
                newy = patchSim.proteinUnits[0]["position"][1] - (com_GDOMAIN[1] - com_FARNESYL_4B[1])/10
                patchSim.proteinUnits[0]["position"] = [newx, newy]
                LOGGER.info("M3 - x,y location shifted to FAR BB, now is = {}".format(patchSim.proteinUnits[0]["position"]))

                # translate z-axes to correct location (somewhat below the membrane) - we want CYF BB to be at -5.455 nm a little below the membrane
                # transkate protein x,y of CYF BB to the correct location (same as lipid location given to insane) 
                #_prot.atoms.translate([com_FARNESYL_4B[0] - 10*newx, com_FARNESYL_4B[1] - 10*newy,  (-1*com_FARNESYL_4B[2]) - (10*5.455)]) 
                _prot.atoms.translate([10*newx - com_FARNESYL_4B[0], 10*newy - com_FARNESYL_4B[1],  10*5.455 - com_FARNESYL_4B[2]])

                LOGGER.info("M3 - farnesyl location after centering should be ({},{},5.455) and is {}".format(newx, newy, _protBB.select_atoms('index 430 to 431').center_of_mass()/10))

                # back to normal stuff (before pre-campaign 3) 
                _prot.atoms.write(gro_file)
                patchSim.proteinUnits[0]["structure"] = gro_file
                # enforce that it's always of type "RAS-RAF"
                patchSim.proteinUnits[0]["type"] = "RAS_RAF_14_3_3"
                # and set charge - note can also be done above but here overwirtten 
                patchSim.proteinUnits[0]["charge"] = -11 # Note insane cacluates the charge but not corectly so this is the offset to get the charge = 0 - actual protein charge is -7
                patchSim.proteinUnits[0]["ff_options"] = added_mdp_options # note these are used in .mdp files to change the protein FF parameters


        elif (patchName[0:9] == "structure"):
            # New for ADMIRRAL C3 - note C2 and mini MuMMI are above 

            # Get ML generated structure 
            inPathInteration = os.path.join(inpath, patchName.split('_')[1]) # ML input + iter dir 
            gro_input = inPathInteration + "/" + patchName + ".gro" # absolute input path
            if not os.path.exists(gro_input):
                LOGGER.warning(f"File {gro_input} does not exist.")
                gro_input = patchName + ".gro"
                if not os.path.exists(gro_input):
                    raise ValueError(f"File {gro_input} does not exist. Aborting")

            LOGGER.info(f"Read in generated input structure after healing = {gro_input}")
            if "mda" not in dir():
                import MDAnalysis as mda
            _prot = mda.Universe(gro_input)
            _protBB = _prot.select_atoms('name BB') 
            #_prot.atoms.write("test1.gro")

            # temporary directory to run macro code. 
            # location of macro template
            tmpl_dir = Naming.dir_res('macro_c3')
            macro_sim_dir = os.path.join(localpath, "macro_sim_dir")

            # should make an isolated sub dir in localpath in which to write gro and run macro sim
            if not os.path.exists(macro_sim_dir):
                LOGGER.info(f"macro_sim_dir not found. Creating a new directory on {localpath}")
                utils.sys_call(f"mkdir {macro_sim_dir}")

            # Check structure for needed FF changes
            # In C3 zn_cutoff and SEP365_cutoff always true (lags set but function not checked)
            added_mdp_options = ""
            # zn_cutoff returns True if it is close enough to coordiante (B state)
            #if zn_cutoff(_prot, 8, patchName):
            if True:
                added_mdp_options += " -DACYM280"
            # 365 cutoff returns True if it is far enough away to not be bound (B state)
            #if SEP365_cutoff(_prot, 8, patchName):
            if True:
                added_mdp_options += " -DASEP365"
            # Get C3 states 
            combination = patchName.split('_')[3]
            if combination[0] == 'Y':
                added_mdp_options += " -DM1DIM"
            if combination[1] == 'Y':
                added_mdp_options += " -DM2DIM"
            if combination[2] == 'Y':
                added_mdp_options += " -DM1M2E"
            if combination[3] == 'Y':
                added_mdp_options += " -DM21433"
            LOGGER.info(f"Options to be added to .mdp file {added_mdp_options}")

            # zero center m1 FAR and align m2 FAR to x-axes 
            com_m1_FARNESYL_4B = _protBB.select_atoms('index 430 to 431').center_of_mass() # should be BB of CYF resnum 185 in m1 KRAS 
            com_m2_FARNESYL_4B = _protBB.select_atoms('index 3327 to 3328').center_of_mass() # should be BB of CYF resnum 185 in m2 KRAS 
            _prot.atoms.translate(-1 * com_m1_FARNESYL_4B)
            dvec = com_m2_FARNESYL_4B - com_m1_FARNESYL_4B
            angle = np.rad2deg(np.arctan2(dvec[1], dvec[0]))
            _prot.atoms.rotateby(-angle, [0.,0.,1.], point=[0,0,0])
            #_prot.atoms.write("test2.gro")

            # center to box middle 
            com_m1_FARNESYL_4B = _protBB.select_atoms('index 430 to 431').center_of_mass() # should be BB of CYF resnum 185 in m1 KRAS 
            com_m2_FARNESYL_4B = _protBB.select_atoms('index 3327 to 3328').center_of_mass() # should be BB of CYF resnum 185 in m2 KRAS 
            xadd = 410/2 - (com_m2_FARNESYL_4B[0] - com_m1_FARNESYL_4B[0]) / 2 
            yadd = 410/2 
            zadd = 10*5.455 # Note this is default inital offset shoudl not be needed  
            # translate z-axes to correct location (somewhat below the membrane) - we want CYF BB to be at -5.455 nm a little below the membrane            
            _prot.atoms.translate((xadd, yadd, zadd))
            #_prot.atoms.write("test3.gro")

            # Get protein posisions for macro model 
            com_m1_GDOMAIN = _protBB.select_atoms('index 0 to 380').center_of_mass() # should be KRAS G-domain 
            com_m1_CRD_4B = _protBB.select_atoms('index 649 to 660 or index 685 to 694').center_of_mass() # this is resid 241 to 244 (THR, PHE, PHE, THR) and resid 254 to 257 (LUE LEU PHE GLN) - CRD membrane interacting part
            com_m1_FARNESYL_4B = _protBB.select_atoms('index 430 to 431').center_of_mass() # should be BB of CYF resnum 185 in KRAS 
            com_m2_GDOMAIN = _protBB.select_atoms('index 2897 to 3277').center_of_mass() # should be KRAS G-domain 
            com_m2_CRD_4B = _protBB.select_atoms('index 3546 to 3557 or index 3582 to 3591').center_of_mass() # this is resid 241 to 244 (THR, PHE, PHE, THR) and resid 254 to 257 (LUE LEU PHE GLN) - CRD membrane interacting part
            com_m2_FARNESYL_4B = _protBB.select_atoms('index 3327 to 3328').center_of_mass() # should be BB of CYF resnum 185 in KRAS 
            dvec = com_m1_CRD_4B - com_m1_GDOMAIN
            angle_m1 = np.rad2deg(np.arctan2(dvec[1], dvec[0]))
            dvec = com_m2_CRD_4B - com_m2_GDOMAIN
            angle_m2 = np.rad2deg(np.arctan2(dvec[1], dvec[0]))
            LOGGER.info(f'Protein far post, com_m1_FARNESYL_4B = {com_m1_FARNESYL_4B}, com_m2_FARNESYL_4B = {com_m2_FARNESYL_4B}')
            LOGGER.info(f'Protein pos for macro model, com_m1_GDOMAIN = {com_m1_GDOMAIN}, angle_m1 = {angle_m1}, com_m2_GDOMAIN = {com_m2_GDOMAIN}, angle_m2 = {angle_m2}')

            # Note, need to change below here for C3 macro model 
            inBstate = True

            # canonical ordering of domains for B state
            domains_dict = {"KD1": _prot.select_atoms('index 1136 to 1768'),
                            "CRD1": _prot.select_atoms('index 637 to 740'),
                            "loop11": _prot.select_atoms('index 753 to 928'),
                            "loop21": _prot.select_atoms('index 955 to 1135'),
                            "b1433": _prot.select_atoms('index 2358 to 2896'),
                            "RAS1": _prot.select_atoms('index 0 to 437'),
                            "KD2": _prot.select_atoms('index 4034 to 4664'),
                            "CRD2": _prot.select_atoms('index 3535 to 3638'),
                            "loop12": _prot.select_atoms('index 3651 to 3826'),
                            "loop22": _prot.select_atoms('index 3853 to 4033'),
                            "RAS2": _prot.select_atoms('index 2898 to 3335')}

            # this is for B state classification only
            codes_to_state = {"000001": "a",
                              "000011": "a",
                              "000101": "a",
                              "000111": "a",
                              "010001": "b",
                              "010011": "b",
                              "010101": "b",
                              "010111": "b",
                              "100001-ypos": "c",
                              "100011-ypos": "c",
                              "100101-ypos": "c",
                              "100111-ypos": "c",
                              "110001-ypos": "d",
                              "110011-ypos": "d",
                              "110101-ypos": "d",
                              "110111-ypos": "d",
                              "100001-yneg": "e",
                              "100011-yneg": "e",
                              "100101-yneg": "e",
                              "100111-yneg": "e",
                              "110001-yneg": "f",
                              "110011-yneg": "f",
                              "110101-yneg": "f",
                              "110111-yneg": "f",
                              "001001": "g",
                              "001011": "g",
                              "001101": "g",
                              "001111": "g",
                              "011001": "b",
                              "011011": "b",
                              "011101": "b",
                              "011111": "b",
                              "101001-ypos": "i",
                              "101011-ypos": "i",
                              "101101-ypos": "i",
                              "101111-ypos": "i",
                              "111001-ypos": "d",
                              "111011-ypos": "d",
                              "111101-ypos": "d",
                              "111111-ypos": "d",
                              "101001-yneg": "k",
                              "101011-yneg": "k",
                              "101101-yneg": "k",
                              "101111-yneg": "k",
                              "111001-yneg": "f",
                              "111011-yneg": "f",
                              "111101-yneg": "f",
                              "111111-yneg": "f",
                              "00001": "a_Bstar",
                              "00011": "a_Bstar",
                              "01001": "b_Bstar",
                              "01011": "b_Bstar",
                              "10001-ypos": "c_Bstar",
                              "10011-ypos": "c_Bstar",
                              "11001-ypos": "d_Bstar",
                              "11011-ypos": "d_Bstar",
                              "10001-yneg": "e_Bstar",
                              "10011-yneg": "e_Bstar",
                              "11001-yneg": "f_Bstar",
                              "11011-yneg": "f_Bstar",
                              "00101": "g_Bstar",
                              "00111": "g_Bstar",
                              "01101": "b_Bstar",
                              "01111": "b_Bstar",
                              "10101-ypos": "i_Bstar",
                              "10111-ypos": "i_Bstar",
                              "11101-ypos": "d_Bstar",
                              "11111-ypos": "d_Bstar",
                              "10101-yneg": "k_Bstar",
                              "10111-yneg": "k_Bstar",
                              "11101-yneg": "f_Bstar",
                              "11111-yneg": "f_Bstar"}

            # protein atom indexes go as M1 is 0 to 2897 inclusive, then M2 is 2898 to 4716
            states_M12 = get_state_from_protein_coordinates(_prot.atoms.positions,
                                                            domains_dict,
                                                            codes_to_state,
                                                            15,
                                                            5,
                                                            [430,3327], 
                                                            com_m1_GDOMAIN,
                                                            com_m2_GDOMAIN)  # farnesyl id is hardcoded as BB


            state = {"M1" : states_M12[0],
                     "X1" : com_m1_GDOMAIN[0], 
                     "Y1" : com_m1_GDOMAIN[1], 
                     "A1" : angle_m1,
                     "M2" : states_M12[1],
                     "X2" : com_m2_GDOMAIN[0], 
                     "Y2" : com_m2_GDOMAIN[1], 
                     "A2" : angle_m2}

            setup_macro_sim_dir(tmpl_dir, macro_sim_dir, state)

            # launch macro job
            run_macro_code(state, tmpl_dir, macro_sim_dir, influx=True) # False for Frontier right now

            # pulls out the last state from from the above run
            p = fetch_gc_patches(localpath, "macro_sim_dir_000000000001", "simple")

            # sets coordinates in patchSim 
            patchSim = convert_gc_patch_to_cg_patch(p, args.simnum, localpath)

            # remove superfluos 2nd ras
            if len(patchSim.proteinUnits) > 1:
                patchSim.proteinUnits = [patchSim.proteinUnits[0]]

            patchSim.proteinUnits[0]["position"] = [com_m1_FARNESYL_4B[0]/10, com_m1_FARNESYL_4B[1]/10]
            patchSim.proteinUnits[0]["position_m2"] = [com_m2_FARNESYL_4B[0]/10, com_m2_FARNESYL_4B[1]/10] # new in C3 for m2 FAR location 
            # location to write gro file from npz input
            gro_file = os.path.join(localpath,"struct_createsim.gro")
            _prot.atoms.write(gro_file)
            patchSim.proteinUnits[0]["structure"] = gro_file
            patchSim.proteinUnits[0]["type"] = "dimercomp" # Full complex RAS-RAF (m1) 14-3-3^2 RAS-RAF (m2)
            patchSim.proteinUnits[0]["charge"] = -4 # Note insane cacluates the charge but not corectly so this is the offset to get the charge = 0
            patchSim.proteinUnits[0]["ff_options"] = added_mdp_options # note these are used in .mdp files to change the protein FF parameters

        else:
            error = "createsims was called with -p {} - value not supported".format(patchName)
            LOGGER.error(error)
            raise ValueError(error)

        # If outpath not set (default in MuMMI) get full paths
        if outpath == "":
            outpath = Naming.dir_sim('cg', simname=patchSim.name)
            LOGGER.info("Output path not set for patch '%s', setting to '%s'", patchSim.name, outpath)
            # Also add full local patch if using local
            #if localpath != "":
            #    localpath = localpath + "/" + utils.add_timestamp(patchSim.name)

        LOGGER.info("final outpath = {}".format(outpath))
        LOGGER.info("final outlocal = {}".format(localpath))

        # Check if previus simulatiosn exist, if so and creatsims_success there, resave sucess and exit with warning
        # if sim folder exist with "many" files but no sucess flag re-name and move - this SHOULD not happen
        if os.path.isdir(outpath):
            if iointerface.test_signal(outpath, flag_success):
                LOGGER.warning("---> FOUND previous sim setup with sim name {} and marked success, will exit".format(patchSim.name))
                marked_failure = False
                exit() # We jump directly to finally block
            else:
                # This should be the normal behavior in MuMMI - but check if there is unexpected files in the folder
                # We do expect at least the current log file
                if len(os.listdir(outpath)) > 8:
                    mvFolder = utils.add_timestamp(outpath + "/{}_OLD_not_success".format(patchSim.name))
                    LOGGER.warning("---> FOUND previous sim setup with sim name {} and NOT marked success, move to {}".format(patchSim.name, mvFolder))
                    #utils.sys_call("rm -r {}".format(outpath))
                    utils.sys_call('mkdir -p ' + mvFolder, test=False)
                    utils.sys_call("mv {}/* {}/ || true".format(outpath, mvFolder))
        else:
            # This should not happen in a MuMMI run as folder is created by worfklow (and log of this file stored there)
            # This is allowed due to testing
            LOGGER.warning("---> Sim {} folder didn't exist and was created - this is not expected within a MuMMI run".format(outpath))
            utils.sys_call('mkdir -p ' + outpath, test=False)

        LOGGER.info("Creating PG Particle simulation...")
        patchSim.create(simpath) # Creates folder and sub-grid concentration files
        patchSim.save() # Saves patchSim object
        if p != None:
            # This is saved locally (local node) in current woking dir and then coped to gloabl file system at the end
            pfPatch.save_npz("ppatch.npz", p)

        setup_cg_sim(patchSim, simpath, mini_mummi)

        LOGGER.info("Running conversion from GROMACS to ddcMD...")
        prime_cg_sim(patchSim, simpath)

        # Setup success
        if iointerface.test_signal(outpath, flag_failure):
            LOGGER.warning(f"Trying to set flag to success ({flag_success}) and found an already existing failure flag ({flag_failure}).")
        else:
            iointerface.send_signal(outpath, flag_success)
        marked_failure = False
        LOGGER.info("---> Createsims finish successfully for sim {}".format(patchSim.name))

        # Log time
        cs_end = time.time()
        LOGGER.info("Creatsims takes = {} seconds".format(cs_end - cs_begin))

    except Exception as error:
        LOGGER.error(error)
        traceback.print_exc()

        '''
        iointerface.send_signal(outpath, flag_failure)
        # Copy back from local and cleanup
        try:
            LOGGER.info("---> Copied form local to global - {}/* {}/".format(simpath, outpath))
            store_files(outpath, simpath, marked_failure)
            cleanup_files(localpath)
        except Exception as error:
            # @TODO - add better error handling and saving of creatsims_error file if always faital
            LOGGER.error("---> ERROR in copying form local to global. {}".format(error))
            iointerface.send_signal(outpath, flag_failure)
            cleanup_files(localpath)
            raise error
        #utils.sys_call(f'cp {os.path.join(simpath, "*.log")} {outpath}/{patchSim.name}_failed.log', test=False)
        '''
        raise error
    finally:
        # Was an error that didn't get cought in an Exception
        if marked_failure:
            LOGGER.error("---> In createsim finally but din't reach success - most likely a gromacs error/exit.")
            if iointerface.test_signal(outpath, flag_success):
                LOGGER.warning(f"Trying to set flag ({outpath}) to failure ({flag_failure}) and found an already existing success flag ({flag_success}). Aborting.")
            else:
                iointerface.send_signal(outpath, flag_failure)
        # Copy back from local and cleanup
        try:
            LOGGER.info("---> Copied form local to global - {}/* {}/".format(simpath, outpath))
            store_files(outpath, simpath, marked_failure)
            LOGGER.info("---> Copy to {} completed.".format(outpath))
            cleanup_files(localpath)
        except Exception as error:
            # @TODO - add better error handling and saving of createsims_error file if always fatal
            LOGGER.error("---> ERROR in copying form local to global. {}".format(error))
            # For now don't delete local
            #if localpath != "":
            #    utils.sys_call("rm -r {}".format(simpath))
            if iointerface.test_signal(outpath, flag_success):
                LOGGER.warning(f"Got error: trying to set flag to failure ({flag_failure}) and found an already existing success flag ({flag_success}).")
            else:
                iointerface.send_signal(outpath, flag_failure)
            cleanup_files(localpath)
            raise error

if __name__ == "__main__":
    main()
