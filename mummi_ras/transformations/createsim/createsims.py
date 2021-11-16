###############################################################################
# @todo add Pilot2-splash-app disclaimer
# @authors  Helgi I. Ingolfsson <ingolfsson1@llnl.gov>
#           Changhoan Kim <kimchang@us.ibm.com>
#           Francesco Di Natale <dinatale3@llnl.gov>
#           @todo add other contributers
###############################################################################

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

import mummi_core
import mummi_ras
from mummi_core.utils import Naming
from mummi_core.utils import utilities as utils
import mummi_ras.transformations.createsim.particlesim as particlesim
import mummi_ras.transformations.createsim.placeproteins as placeproteins
import mummi_ras.transformations.createsim.pdcglob as pdcglob
from mummi_ras.transformations.patch_creator import fetch_gc_patches
from mummi_ras.datastructures.patch import Patch as pfPatch
from mummi_core.utils.logger import init_logger

# Program Globals
LOGGER = getLogger(__name__)
global gromacs, mdrun_mpi, mdrun_o


def setup_cg_sim(pSim, simpath):
    """Set up and run eq for the CG sim."""
    LOGGER.info("---> Setup CG sim and run eq in dir {}".format(simpath))
    os.chdir(simpath)

    # Get protein structures insertion string and add random rotation / placement as needed
    protein_string = ""
    # @TODO put this back in with alignment but also support no allignment if protein sturcture is specified (use below)

    if (len(pSim.proteinUnits) > 0):
        proteinRotated, proteinPosInserted = placeproteins.protein_align(pSim.proteinUnits)
        LOGGER.info("---> Proteins are pressent they are: {}".format(pSim.proteinUnits))
        LOGGER.info("---> After protein alignment: rotation {} posInserted {}".format(proteinRotated, proteinPosInserted))
        for i in range(len(pSim.proteinUnits)):
            cyfGroupX = -1
            cyfGroupY = -1
            cypGroupX = -1
            cypGroupY = -1
            if pSim.proteinUnits[i]["type"] == "RAS4A-RAF" or pSim.proteinUnits[i]["type"] == "RAS4A-ONLY":
                # Have kras4A need to fined location of CYF and CYP to feed to insane
                import MDAnalysis as mda
                protein = mda.Universe(proteinRotated[i])
                farCOM = protein.select_atoms("resname CYF and name SC1").center_of_mass() / 10 # Rember to convert from A to nm
                palCOM = protein.select_atoms("resname CYP and name SC1").center_of_mass() / 10 # Rember to convert from A to nm
                LOGGER.info("---> RAS4A CYF = {},{} vs CYP = {},{} - COM = {},{}".format(farCOM[0], farCOM[1], palCOM[0], palCOM[1], proteinPosInserted[i][0], proteinPosInserted[i][1]))
                cyfGroupX = farCOM[0]
                cyfGroupY = farCOM[1]
                cypGroupX = palCOM[0]
                cypGroupY = palCOM[1]
            else:
                # Have kras4B need to fined location of CYF only leave cyp as -1
                cyfGroupX = proteinPosInserted[i][0]
                cyfGroupY = proteinPosInserted[i][1]
                cypGroupX = -1
                cypGroupY = -1

                # Test only - remove after checking as MDAnalysis takes time to load
                if False:
                    import MDAnalysis as mda
                    protein = mda.Universe(proteinRotated[i])
                    farCOM = protein.select_atoms("resname CYF and name SC1").center_of_mass() / 10 # Rember to convert from A to nm
                    LOGGER.info("---> RAS4B CYF = {},{} - COM = {},{}".format(farCOM[0], farCOM[1], proteinPosInserted[i][0], proteinPosInserted[i][1]))


            # utils.sys_call("cp "+pSim.rasList[i]+" KRAS-"+str(i)+".gro")
            # utils.sys_call(gromacs+" editconf -f KRAS-"+str(i)+".gro -o KRAS-"+str(i)+"-rotate.gro -rotate 0 0 "+str(randint(1, 360)))
            # kras_string += " -ras KRAS-"+str(i)+"-rotate.gro l " + str(pSim.rasPos[i][0]) + " " + str(pSim.rasPos[i][1]) + " "
            # kras_string += " -ras " + pSim.rasList[i] + " l " + str(pSim.rasPos[i][0]) + " " + str(pSim.rasPos[i][1]) + " "
            # kras_string += " -ras " + krasRotated[i] + " l " + str(krasPosInserted[i][0]) + " " + str(krasPosInserted[i][1]) + " "

            # @TODO fix this so it will work generalily for C4 - either both insane's or merge support into x1 insane
            #protein_string += " -ras " + proteinRotated[i] + " l " + str(proteinPosInserted[i][0]) + " " + str(proteinPosInserted[i][1]) + " "
            protein_string += " -ras " + proteinRotated[i] + " l " + str(cyfGroupX) + " " + str(cyfGroupY) + " " + str(cypGroupX) + " " + str(cypGroupY) + " "
    '''
    for i in range(len(pSim.proteinUnits)):
        protein_string += " -ras " + pSim.proteinUnits[i]["structure"] + " l " + str(pSim.proteinUnits[i]["position"][0]) + " " + str(pSim.proteinUnits[i]["position"][1]) + " "
    '''
    LOGGER.info("---> Manual protein alignment insane protein string = {}".format(protein_string))

    # If inital setup / em fails tray again untill max tries
    maxTries = 10
    for i in range(maxTries):
        LOGGER.info("---> System setup/EM try number {}".format(i))

        # Make inital cordinates (file overwritten by insane so need to generate each time)
        utils.get_file(particlesim.dirMartini+"/"+pSim.defTopFile, 'system.top')
        utils.sys_call("sed -ie 's/xPATHx/"+particlesim.dirMartini.replace("/", "\/")+"/g' system.top")
        utils.sys_call("sed -ie 's/YYY/"+pSim.rasItp.replace("/", "\/")+"/g' system.top")
        utils.sys_call("sed -ie 's/ZZZ/"+pSim.rasrafItp.replace("/", "\/")+"/g' system.top")
        utils.sys_call("sed -ie 's/WWW/"+pSim.ras4aItp.replace("/", "\/")+"/g' system.top")
        utils.sys_call("sed -ie 's/XXX/"+pSim.ras4arafItp.replace("/", "\/")+"/g' system.top")

        # WARNING if we have many CPUs compaired to the simulation size we will get into domain problems - expecially in the begining
        # so change -nt depending on system size e.g. for 15x15 nm system 12 CPUs is good initally
        ntString = ""  # default is empty ie let GROMACS use max CPUs
        if pSim.cgnt > 0:
            ntString = " -nt "+str(pSim.cgnt)+" "

        # Get total protein charge - this is now stored for each protein unit
        # Note: use -2 for RAS and -1 for RAS_RAF WARNING RAS-RAF is missing a carge of +4 this not included here so ions need to be addjusted in backmapping
        pCharge = int(np.sum([x['charge'] for x in pSim.proteinUnits]))
        print(pCharge)

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
            utils.sys_call('python '+particlesim.dirMartini+'/insane-patch-build.py' +
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
            # Assign lipids based on probabilities as specified in -uf and -lf files
            utils.sys_call('python '+particlesim.dirMartini+'/insane-ck.py' +
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
            utils.sys_call("cp system.top system-pull.top")

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
                    utils.sys_call("sed -ie '0,/Protein          1/s//RASR_PULL        1/g' system-pull.top")
                    utils.sys_call("sed -ie '0,/Protein          1/s//RAS_RAF          1/g' system.top")
                elif pSim.proteinUnits[j]["type"] == "RAS4A-ONLY":
                    utils.sys_call("sed -ie '0,/Protein          1/s//RAS4A_PULL       1/g' system-pull.top")
                    utils.sys_call("sed -ie '0,/Protein          1/s//RAS4A            1/g' system.top")
                elif pSim.proteinUnits[j]["type"] == "RAS4A-RAF":
                    utils.sys_call("sed -ie '0,/Protein          1/s//RAS4A_RAF_PULL   1/g' system-pull.top")
                    utils.sys_call("sed -ie '0,/Protein          1/s//RAS4A_RAF        1/g' system.top")
                else:
                    error = "Protein type {} not supported".format(pSim.proteinUnits[j]["type"])
                    LOGGER.error(error)
                    raise AttributeError(error)

           # Run energy minimization (steepest descent, with protein constreins)
            utils.sys_call(gromacs+' grompp -c lipids-water.gro -r lipids-water.gro -f '+particlesim.dirMartini +
                           '/martini_v2.x_new-rf-em.mdp -p system-pull.top -o topol.tpr -maxwarn 5')
            # @TODO fix this hack? ntString can't be used as we don't want to limit the CPU's of the proper md runs below
            #utils.sys_call(gromacs+' mdrun '+ntString+' -v -c lipids-water-em.gro >> mdrun-em.log 2>&1')
            # @TODO need to change this
            utils.sys_call(gromacs+' mdrun -ntmpi 20 -v -c lipids-water-em.gro >> mdrun-em.log 2>&1')
            #utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-em.gro >> mdrun-em.log 2>&1 ')
            utils.sys_call('echo 0 | ' + gromacs+' trjconv -f lipids-water-em.gro -o lipids-water-em.gro -pbc mol -s topol.tpr')
        else:  # No KRAS bilayer only simulations
            # Run energy minimization (steepest descent, with protein constreins)
            utils.sys_call(gromacs+' grompp -c lipids-water.gro -r lipids-water.gro -f '+particlesim.dirMartini +
                           '/martini_v2.x_new-rf-em.mdp -p system.top -o topol.tpr -maxwarn 5')
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
        tempOutput.write("r CYF & a SC1\n")
        tempOutput.write("name 3 PULL1\n")
        tempOutput.write("r CYF & a F1\n")
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
        # @TODO later we migth also have to addjust pull length based on KRAS structure ?
        # pullMDPfile = "martini_v2.x_new-rf-eq3-PW-PULL.mdp"
        pullMDPfile = "martini_v2.x_new-rf-eq3-PULL.mdp"
        utils.sys_call("cp "+particlesim.dirMartini+"/"+pullMDPfile+" .")


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
            tempOutput.write("pull_group"+str(i*2+1)+"_pbcatom      = 0 \n")
            tempOutput.write("pull_group"+str(i*2+2)+"_name         = PULL"+str(i*2+2)+" \n")
            tempOutput.write("pull_group"+str(i*2+2)+"_pbcatom      = 0 \n")
            tempOutput.write("pull_coord"+str(i+1)+"_type         = flat-bottom \n")
            tempOutput.write("pull_coord"+str(i+1)+"init          = 0.0 \n")  # was 0.4
            tempOutput.write("pull_coord"+str(i+1)+"_geometry     = distance  \n")
            tempOutput.write("pull_coord"+str(i+1)+"_groups       = "+str(i*2+1)+" "+str(i*2+2)+" \n")
            tempOutput.write("pull_coord"+str(i+1)+"_dim          = N N Y  \n")
            tempOutput.write("pull_coord"+str(i+1)+"_rate         = -0.00045  \n") # [nm/ps]
            # Make less to pull slower?
            # 2021.07.02 400000 * 0.01 = 4000 ps, -0.00045 * 4000 = 1.8 nm. They should all start with a 2.1 nm diff
            # Could change and over pull using e.g. -0.0005 - might help when there is a CYF and CYP mismatch
            # in current martini_v2.x_new-rf-eq3-PW-PULL.mdp using 2,000,000 * 0.01 ps = 20,000 ps
            #   or -0.0003 nm/ps * 20,000 ps = -6 nm
            # This is a little bit more than is needed in the current test simulations
            # @TODO for production reduce box size have these somewhat closer and pull slower
            tempOutput.write("pull_coord"+str(i+1)+"_k            = 1000 \n")
            tempOutput.write("pull_coord"+str(i+1)+"_start        = yes  \n")
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
        utils.sys_call(gromacs+' grompp -c lipids-water-em.gro -f '+particlesim.dirMartini+'/martini_v2.x_new-rf-eq1.mdp ' +
                       ' -p system-pull.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-em.gro')
        utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq1.gro >> mdrun-eq1.log 2>&1 ')

        # Run short equilibration at 5 fs time setp
        # WARNING for PW use martini_v2.x_new-rf-eq2-PW.mdp
        utils.sys_call(gromacs+' grompp -c lipids-water-eq1.gro -f '+particlesim.dirMartini+'/martini_v2.x_new-rf-eq2.mdp ' +
                       ' -p system-pull.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-em.gro')
        utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq2.gro >> mdrun-eq2.log 2>&1')

        # Pull molecules to bilayer here, then change .top to while KRAS molecules - using local PULL .mdp
        utils.sys_call(gromacs+' grompp -c lipids-water-eq2.gro -f '+pullMDPfile +
                       ' -p system-pull.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-em.gro')
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
        utils.sys_call(gromacs+' grompp -c lipids-water-em.gro -f '+particlesim.dirMartini +
                       '/martini_v2.x_new-rf-eq1.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-em.gro')
        utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq1.gro >> mdrun-eq1.log 2>&1')

        # Run short equilibration at 5 fs time setp
        # WARNING for PW use martini_v2.x_new-rf-eq2-PW.mdp
        utils.sys_call(gromacs+' grompp -c lipids-water-eq1.gro -f '+particlesim.dirMartini +
                       '/martini_v2.x_new-rf-eq2.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-em.gro')
        utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq2.gro >> mdrun-eq2.log 2>&1')

        # Run short equilibration at 10 fs time setp
        # WARNING for PW use martini_v2.x_new-rf-eq3-PW.mdp
        utils.sys_call(gromacs+' grompp -c lipids-water-eq2.gro -f '+particlesim.dirMartini +
                       '/martini_v2.x_new-rf-eq3.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-em.gro')
        utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq3.gro >> mdrun-eq3.log 2>&1')

    ## @TODO remove this in Campain 3 - this is for campain 1*, to fix the HVR allow it to fix the helix - also edit command below
    ## This run shold allow the HVR to move - fixing the few residues that need to form a helix
    #utils.sys_call(gromacs+' grompp -c lipids-water-eq3.gro -f '+particlesim.dirMartini +
    #               '/martini_v2.x_new-rf-eq4hvr.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-eq3.gro')
    ## To troubleshot topolgy add -pp temp-pp.top
    #utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq4hvr.gro >> mdrun-eq4hvr.log 2>&1')

    # Added a x2 pre eq4 setp eq4p - 1/10 the number of time step using 5fs and 1/5 the number of time steps, using 10fs
    utils.sys_call(gromacs+' grompp -c lipids-water-eq3.gro -f '+particlesim.dirMartini +
                   '/martini_v2.x_new-rf-eq4p.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-eq3.gro')
    utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq4p.gro >> mdrun-eq4p.log 2>&1')
    #utils.sys_call(mdrun_mpi+' -nt 24 -rdd 2.0 -ntomp 4 -c lipids-water-eq4p.gro >> mdrun-eq4p.log 2>&1')
    utils.sys_call(gromacs+' grompp -c lipids-water-eq4p.gro -f '+particlesim.dirMartini +
                   '/martini_v2.x_new-rf-eq4pp.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-eq4p.gro')
    utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq4pp.gro >> mdrun-eq4pp.log 2>&1')


    # Run final equilibration at 20 fs time setp (same as pruduction - test stability + allow the box are to equilibrate
    # @TESTING for now do only 2 ns but should move this to 20-100 ns (for the 3000 lipid systesm)
    # WARNING for PW use martini_v2.x_new-rf-eq4-PW.mdp
    #utils.sys_call(gromacs+' grompp -c lipids-water-eq3.gro -f '+particlesim.dirMartini +
    #               '/martini_v2.x_new-rf-eq4.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-eq3.gro')
    #utils.sys_call(gromacs+' grompp -c lipids-water-eq4hvr.gro -f '+particlesim.dirMartini +
    #               '/martini_v2.x_new-rf-eq4.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-eq4hvr.gro')
    utils.sys_call(gromacs+' grompp -c lipids-water-eq4pp.gro -f '+particlesim.dirMartini +
                   '/martini_v2.x_new-rf-eq4.mdp -p system.top -o topol.tpr -maxwarn 5 -n index.ndx -r lipids-water-eq4pp.gro')
    # To troubleshot topolgy add -pp temp-pp.top
    utils.sys_call(mdrun_mpi+' '+mdrun_o+' '+ntString+' -c lipids-water-eq4.gro >> mdrun-eq4.log 2>&1')

    # Cleanup temp files
    utils.sys_call('rm \#* step*.pdb || true')
    utils.sys_call('rm system.tope system-pull.tope state.cpt state_prev.cpt traj_comp.xtc traj.trr mdout.mdp md.log ener.edr index-old.ndx temp-index-end.txt || true')

    # Compress all cretsims files into createsims.tar.gz storage that are not needed for next steps and remove
    utils.sys_call('tar -czvf createsims.tar.gz psimSubPatch_lower.pat psimSubPatch_upper.pat protein-*-rotate.gro protein-ALL.gro ' +
                   'index-selection.txt lipids-water-em.gro lipids-water-eq1.gro lipids-water-eq2.gro lipids-water-eq3.gro ' +
                   'lipids-water.gro martini_v2.x_new-rf-eq3-PULL.mdp mdrun-em.log mdrun-eq1.log mdrun-eq2.log mdrun-eq3.log ' +
                   'mdrun-eq4.log pullf.xvg pullx.xvg system-pull.top lipid_counts.npz mdrun-eq4p.log lipids-water-eq4p.gro ' +
                   'mdrun-eq4pp.log lipids-water-eq4pp.gro || true')
    utils.sys_call('rm psimSubPatch_lower.pat psimSubPatch_upper.pat protein-*-rotate.gro protein-ALL.gro ' +
                   'index-selection.txt lipids-water-em.gro lipids-water-eq1.gro lipids-water-eq2.gro lipids-water-eq3.gro ' +
                   'lipids-water.gro martini_v2.x_new-rf-eq3-PULL.mdp mdrun-em.log mdrun-eq1.log mdrun-eq2.log mdrun-eq3.log ' +
                   'mdrun-eq4.log pullf.xvg pullx.xvg system-pull.top lipid_counts.npz mdrun-eq4p.log lipids-water-eq4p.gro ' +
                   'mdrun-eq4pp.log lipids-water-eq4pp.gro || true')


def run_multible_cg_sim(pSim, simpath, simrange, genvel):
    """Run multible production CG simulation."""
    for i in simrange:
        run_cg_sim(pSim, simpath, rName="-r"+str(i), genvel=genvel)


def prime_cg_sim(pSim, simpath, rName="", genvel=False):
    """Run production CG simulation."""
    LOGGER.info("---> Setup CG production sim in dir {}".format(simpath))
    os.chdir(simpath)

    # Copy over last eq frame and topology
    cgLastFrame = "lipids-water-eq4.gro"
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
        utils.sys_call(gromacs+' grompp -c '+cgLastFrame+' -r '+cgLastFrame+' -f '+particlesim.dirMartini +
                       '/martini_v2.x_new-rf-prod-vel.mdp -p '+cgTopFile+' -o topol.tpr -maxwarn 5 -n '+cgIndexFile)
    else:
        utils.sys_call(gromacs+' grompp -c '+cgLastFrame+' -r '+cgLastFrame+' -f '+particlesim.dirMartini +
                       '/martini_v2.x_new-rf-prod.mdp -p '+cgTopFile+' -o topol.tpr -maxwarn 5 -n '+cgIndexFile)

    # Get SLURM submission file
    #utils.get_file(particlesim.dirMartini+'/fancy-slurm-syrah.sh', 'fancy-slurm-syrah.sh')
    #utils.sys_call("sed -ie 's/XXX/"+pSim.descriptiveName+rName+"/g' fancy-slurm-syrah.sh")

    # Get quartz submission file
    #utils.get_file(particlesim.dirMartini+'/fancy.sh', 'fancy.sh')
    #utils.sys_call("sed -ie 's/XXX/"+pSim.descriptiveName+rName+"/g' fancy.sh")

    # Make .pdb file for use as input in ddcMD
    pdbFileName = cgLastFrame[0:-3]+"pdb"
    utils.sys_call(gromacs+' editconf -f '+cgLastFrame+' -o '+pdbFileName)

    ###############
    # Get all files needed for ddcMD Martini run

    utils.sys_call("mkdir snapshot.mem")
    #utils.sys_call("mkdir snapshot.000000000000")

    # Fix .pdb file to include protein seperation section
    pdbFileFileName = cgLastFrame[0:-4]+"-fix.pdb"
    utils.sys_call("cp "+pdbFileName+" "+pdbFileFileName)
    utils.sys_call("sed -ie 's/MG+ ION/MG   MG/g' "+pdbFileFileName)
    utils.sys_call("sed -ie 's/NA+ ION/NA   NA/g' "+pdbFileFileName)
    utils.sys_call("sed -ie 's/CL- ION/CL   CL/g' "+pdbFileFileName)
    utils.sys_call("sed -ie 's/ZN+ ION/ZN   ZN/g' "+pdbFileFileName)
    #utils.sys_call("sed -ie '/BB  THR     1/i TER' "+pdbFileFileName)
    #utils.sys_call("sed -ie '/BB  THR     1/i PROTEIN        1     KRAS' "+pdbFileFileName)
    #utils.sys_call("sed -ie '/MG  MG/a ENDPROTEIN' "+pdbFileFileName)

    # Fix names in .top file if needed
    utils.sys_call("sed -ie 's/NA+ /NA  /g' "+cgTopFile)
    utils.sys_call("sed -ie 's/CL- /CL  /g' "+cgTopFile)

    # Copy in needed ddcMD setup files
    utils.sys_call("cp "+particlesim.dirddcMD+"/object.data .")
    utils.sys_call("cp "+particlesim.dirddcMD+"/resItpList .")
    utils.sys_call("cp "+particlesim.dirddcMD+"/POPX_Martini_v2.0_lipid.itp .")

    # Copy in ddcMD Martini toloplogy files - note this can be default form mummi_resources or feedback from specific feedback dir.
    utils.sys_call("cp "+pSim.ddcMDParams+"/martini.data .")
    utils.sys_call("cp "+pSim.ddcMDParams+"/molecule.data .")
    utils.sys_call("cp "+pSim.ddcMDParams+"/ConsAtom.data .")

    # Convert final pdb to ddcMD atom format + make retart file and put in correct places
    #   This is the new way of doing this move back to that when ddcMD new vesion is working
    ## @TODO @WARNING this is a hack to get proItpList to work; Xiaohua, needs to change - maybe we change from system call to pyton call?
    #utils.sys_call("echo  " + particlesim.dirddcMD + "/KRAS-GTP-04-HVR-best-guess-M22-CYFpos.itp > proItpList")
    #utils.sys_call('pdbmartini2obj -p '+pdbFileFileName+' -t proItpList -f '+particlesim.dirddcMD+'/ConsAtom.data -o snapshot.mem/atoms#000000')
    ##   This is the old way of calling the ddcMDConvertor - used now for testing
    #utils.sys_call('pdbmartini2obj -p '+pdbFileFileName+' -o snapshot.mem/atoms#000000')
    # This is the new way convert 1.0.4
    utils.sys_call('pdbmartini2obj -p '+pdbFileFileName+' -t '+cgTopFile+' -m martini.data -f ConsAtom.data -o snapshot.mem/atoms#000000')

    utils.sys_call("mv restart snapshot.mem/")
    utils.sys_call("ln -s snapshot.mem/restart")
    #utils.sys_call('pdbmartini2obj -p '+pdbFileFileName+' -o snapshot.000000000000/atoms#000000')
    #utils.sys_call("mv restart snapshot.000000000000/")
    #utils.sys_call("ln -s snapshot.000000000000/restart")

    # Generate ddcMD restraint file
    #@TODO resItpList is hardcoded also the file it refearnces - needs to change more to resources - and/or compleatly change
    #@WARNING works on lassen - will not work on summit
    #utils.sys_call("restraint -i snapshot.mem/atoms#000000 -o restraint.data -p /usr/gapps/kras/templates/ddcMD/resItpList")
    #utils.sys_call("restraint -i snapshot.000000000000/atoms#000000 -o restraint.data -p /usr/gapps/kras/templates/ddcMD/resItpList")
    utils.sys_call("restraint -i snapshot.mem/atoms#000000 -o restraint.data -p resItpList")

    # Cleanup temp files
    utils.sys_call('rm \#* || true')
    utils.sys_call('rm fancy-slurm-syrah.she fancy.she lipids-water-eq4-fix.pdbe mdout.mdp system.tope || true')


def run_cg_sim_GROMACS(pSim, simpath, inpath, intype):
    """Run GROMACS in a alreaddy setup directory."""
    # @WARNING this on was used for testing and is currently not working @TODO fix
    LOGGER.info("---> Run CG sim in dir {}".format(simpath))
    os.chdir(simpath)

    # To troubleshot topolgy add -pp temp-pp.top
    utils.sys_call(mdrun_mpi+' '+mdrun_o+' -maxh 0.3 -c final.pdb >> mdrun-prod.log 2>&1')

    # Dump pdb file from Sim
    utils.sys_call('echo 0 | ' + gromacs+' trjconv -f traj_comp.xtc -o '+pSim.name+'_cg000_f.pdb -sep')

    '''
    # Load pdb fiels to DataBroker
    if (intype == "db"):
        # Save .pdb files to DB
        db = DataBroker()
        db.connect(inpath)
        for pdbFile in os.listdir('.'):
            if os.path.isfile(pdbFile) and pdbFile[0:7] == "pfpatch":
                LOGGER.info("PDB_File: " + pdbFile)
                db.put(pdbFile, particlesim.add_simpath(pdbFile))
        db.disconnect()
    '''

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

    # initram uses GROMACS and has to have it in path - now source, on the server we might want to change initram
    # @TODO this dose not work for some reason ??? currently need to source gromacs before running script,
    #  change to pass gmx path to initram.sh
    # utils.sys_call('source '+gromacs+'GMXRC')

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


def creat_cg_patch(specifications, simnum):
    """Creates a default mix8 one KRAS patch, then modifies based on specifications."""

    # Average Mix8b - no sub-composision
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
        cProteinDict["structure"] = "/g/g90/helgi/mummi_resources/martini/test_RAS4A/KRAS4A_RAF_shift.gro"
        cProteinDict["charge"] = -1
    elif "RAS4A" in simnum:
        # RAS4A
        cProteinDict["type"] = "RAS4A-ONLY"
        cProteinDict["beads"] = -1
        cProteinDict["state"] = "N/A"
        cProteinDict["position"] =  [15,15]
        cProteinDict["structure_pre_assigned"] = True
        cProteinDict["structure"] = "/g/g90/helgi/mummi_resources/martini/test_RAS4A/KRAS4A_shift.gro"
        cProteinDict["charge"] = -2
    elif "RAS_RAF" in simnum:
        # RAS-RAF (with regular KRAS4b)
        cProteinDict["type"] = "RAS-RAF"
        cProteinDict["angle"] = 0
        cProteinDict["beads"] = -1
        cProteinDict["state"] = "N/A"
        cProteinDict["position"] =  [15,15]
        cProteinDict["structure_pre_assigned"] = True
        cProteinDict["structure"] = "/g/g90/helgi/mummi_resources/martini/test_RAS4A/RAS-RAF_state_ma_1.gro"
        cProteinDict["charge"] = -1
    else:
        # RAS-only (regular KRAS4b)
        cProteinDict["type"] = "RAS-ONLY"
        cProteinDict["beads"] = -1
        cProteinDict["state"] = "N/A"
        cProteinDict["position"] =  [15,15]
        cProteinDict["structure_pre_assigned"] = True
        cProteinDict["structure"] = "/g/g90/helgi/mummi_resources/martini/test_RAS4A/RAS-ONLY_state_a_100.gro"
        cProteinDict["charge"] = -2

    proteinUnits.append(cProteinDict)

    patchSim = particlesim.ParticleSim(lipidTypes, lipidUpperRatioArray, lipidLowerRatioArray, lipidSuPatches=[1, 1], proteinUnits=proteinUnits, asym=-118, simName=simnum, lipidFullAsignment=False)
    #patchSim.cgnt = 8 not supportet if you also specified in commond input
    # @TODO generalize RAS-RAF temp test
    #krasStruct = particlesim.dirKRASStructures + "/KRAS-04-protein-cg-M22-em-shift.gro"
    #krasStruct = "/g/g90/helgi/mummi_resources/kras-structures/RAS-RAF-01-fix-shift.gro"
    #patchSim.rasList = [krasStruct]
    #patchSim.rasPos = [[15,15]]

    # Alter default based on specifications
    for cAttr in specifications.split():
        [key,value] = cAttr.split('=')
        setattr(patchSim, key, value)

    return patchSim

def convert_gc_patch_to_cg_patch(p, simnum):
    """Converts a GC patch to a CG patch."""

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
    lipidTypes = ["POPX", "POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"]

    lInner = lconcentrations[:, :, 0:8]
    lOuter = lconcentrations[:, :, 8:14]
    sumOuter = np.sum(lOuter) / (lconcentrations.shape[0]*lconcentrations.shape[1])
    sumInner = np.sum(lInner) / (lconcentrations.shape[0]*lconcentrations.shape[1])
    asym = np.rint(1600 * (1 - (sumInner / sumOuter))).astype(int)
    # For testing of lipid consentrations
    #np.save("concentrations_asym.npy", asym)

    # Save full lipid type asignmet patch for insane to read
    lipid_types_inner = ["POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"]
    lipid_types_outer = ["POPX", "PAPC", "POPE", "DIPE", "DPSM", "CHOL"]
    #p.concentrations.shape # (37, 37, 14)  (for new patches)
    #lconcentrations_full = pdcglob.reinterp(37, 37, 14, p.concentrations, 40, 40)
    # To fit protein placement we need to use the transpose of x,y columns
    lconcentrations_full = pdcglob.reinterp(72, 72, 14, p.concentrations.transpose((1,0,2)), 40, 40)
    lInner_full = lconcentrations_full[:, :, 0:8]
    lOuter_full = lconcentrations_full[:, :, 8:14]
    sumOuter_full = np.sum(lOuter_full) / (lconcentrations_full.shape[0]*lconcentrations_full.shape[1])
    sumInner_full = np.sum(lInner_full) / (lconcentrations_full.shape[0]*lconcentrations_full.shape[1])

    # @WARNING should not be asym_full always 0
    asym_full = np.rint(1600 * (1 - (sumInner_full / sumOuter_full))).astype(int)
    LOGGER.info("Both asym(s) shoudl be the same asym {} and asym_full {} outer {} inner {}".format(asym, asym_full, sumOuter_full, sumInner_full))
    lipid_counts_inner = pdcglob.pdcglob(40, 40, 8, lInner_full)
    lipid_counts_outer = pdcglob.pdcglob(40, 40, 6, lOuter_full)
    lipid_counts_file = "lipid_counts.npz"
    np.savez_compressed(lipid_counts_file, lipid_counts_inner = lipid_counts_inner, lipid_counts_outer = lipid_counts_outer, asym = asym_full, lipid_types_inner = lipid_types_inner, lipid_types_outer = lipid_types_outer)
    LOGGER.info("Full patch lipid type assignment saved to {}, asym {}, shape {}".format(lipid_counts_file, asym_full, lconcentrations_full.shape))


    if (lInner.min() < 0 or lOuter.min() < 0):
        error = "---> Negative lipid consentrations found "
        LOGGER.error(error)
        raise ValueError(error)

    if (asym > 500 or asym < -500):
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

    # Convert protein bead list to protein "unit" list - where each unit is inserted a one piece later in the pipeline
    # For each protein unit - get inital structure, insertion properties (location, angle etc), and charge

    LOGGER.info("Protein bead ids    = {}".format(p.protein_ids))
    LOGGER.info("Protein bead states = {}".format(p.protein_states))
    LOGGER.info("Protein bead pos.   = {}".format(p.protein_positions))
    LOGGER.info("Protein complex ids = {}".format(p.complex_ids))
    #LOGGER.info("Protein comple dis. = {}".format(p.complex_dists)) # removed from storage!

    # Initial shift not needed anymore as first beads is at center pixel (not 0,0 but within 1/2 pixel of that)
    # posArray = p.protein_positions - p.protein_positions[0]  # Center at 0,0

    '''
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
    '''
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
        if (len(p.complex_ids) > 0 and any(remainig_ids[0] in string for string in p.complex_ids)):
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
            deltaPost = correctPBCjump(deltaPost, [30, 30]) # in current GCMM schema beads are allways within PBC here we don't want that but the RAF that is closest to RAS
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
                #cProteinDict["structure"] = "/g/g90/helgi/mummi_resources/martini/test_RAS4A/KRAS4A_RAF_shift_updateitp.gro"
                #cProteinDict["structure"] = "/g/g90/helgi/mummi_resources/martini/test_RAS4A/KRAS4A_RAF_shift.gro"
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
                #cProteinDict["structure"] = "/p/gpfs1/helgi/mummi_root_20210513/sims-cg/RAS4A_states/temp-box-shift-em.gro"
                #cProteinDict["structure_pre_assigned"] = True
                #cProteinDict["structure"] = "/g/g90/helgi/mummi_resources/martini/test_RAS4A/KRAS4A_shift.gro"
                #cProteinDict["structure"] = "/g/g90/helgi/mummi_resources/martini/test_RAS4A/KRAS4A_shift_updateitp.gro"
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
    lipid_counts_file = "lipid_counts.npz"
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
                #cProteinDict["structure"] = "/g/g90/helgi/mummi_resources/martini/test_RAS4A/KRAS4A_RAF_shift_updateitp.gro"
                #cProteinDict["structure"] = "/g/g90/helgi/mummi_resources/martini/test_RAS4A/KRAS4A_RAF_shift.gro"
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
                #cProteinDict["structure"] = "/g/g90/helgi/mummi_resources/martini/test_RAS4A/KRAS4A_shift.gro"
                #cProteinDict["structure"] = "/g/g90/helgi/mummi_resources/martini/test_RAS4A/KRAS4A_shift_updateitp.gro"
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

                ''' @TODO, this
                # Temp add for RAS-ONLY seperate build using old patches
                cProteinDict["beads"] = 1
                cProteinDict["state"] = "n/a"
                cProteinDict["position"] =  posArray[0]
                assert cProteinDict["position"].shape[0] == 2
                #/p/gpfs1/helgi/init_structures_ras/chris_fixed_10/sr2_pfpatch_000000343624_1000ns.gro
                cStructureName = "/p/gpfs1/helgi/init_structures_ras/chris_fixed_10/{}".format(simnum)
                LOGGER.info("Get {} structure: state {}, filename {}".format(cProteinDict["type"], cProteinDict["state"], cStructureName))
                if not os.path.isfile(cStructureName):
                    error = "---> Protein structure file {} not found.".format(cStructureName)
                    LOGGER.error(error)
                    raise ValueError(error)
                cProteinDict["structure"] = cStructureName
                # Temp end add for RAS-ONLY seperate build using old patches
                '''
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

    ''' @TODO fix this and make it general for RAS and RAF
    # Select KRAS structures acording to states - pfPatch state to saved files + random value (e.g. 1-1000)
    # "KRAS-s"+str(cS)+"-r"+str(randint(1, librarySize))+".gro"
    #librarySize = 1000  # numer of diffrent configurations saved for each state
    #krasStructs = []
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

    # @TODO remove RAS-RAF temp test - NOW this is RAS only campain 1* structs
    krasStructs = []
    #krasStruct = particlesim.dirKRASStructures + "/KRAS-04-protein-cg-M22-em-shift.gro"
    #krasStruct = "/g/g90/helgi/mummi_resources/kras-structures/RAS-RAF-01-fix-shift.gro"
    structNumber = simnum.split("_")[2] # get XXX number only RAS_RAF_XXX
    #structPath = "/p/gpfs1/helgi/init_structures_ras_craf/crafForHelgi_2020June20/sel_1000_weight_random"
    structPath = "/p/gpfs1/helgi/init_structures_ras/converterd"
    fileMatch = fnmatch.filter(os.listdir(structPath), '{}_*.gro'.format(structNumber))
    krasStruct = "{}/{}".format(structPath, fileMatch[0])
    for i in range(len(p.protein_states)):
        krasStructs.append(krasStruct)
    #print("s num {} file {}".format(structNumber, krasStruct))
    ## @TODO remove this - it'shere Make this a none RAS patch:
    #krasStructs = []
    #krasPosArray = []
    #exit()
    '''

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

    return parser


def store_files(outpath: str, sim_path: str, failure: bool):
    if sim_path != outpath and os.path.exists(sim_path):
        dest_path = f'{outpath}/{utils.add_timestamp("failed")}' if failure else outpath
        LOGGER.debug("---> Copy local %s to global folder %s", sim_path, dest_path)
        os.makedirs(dest_path, exist_ok=True)
        utils.sys_call("cp -r {}/* {}/".format(sim_path, dest_path))


def cleanup_files(path: str):
    if path and os.path.exists(path):
        return # @TODO put this back in commented out for local debug
        #utils.sys_call(f"rm -r {path}")


def sigterm_wrapper(outpath: str, local_path: str):

    def sigterm_handler(signum, frame):
        LOGGER.info(f"Caught {signum} -- Aborting.")
        store_files(outpath, local_path, True)
        cleanup_files(local_path)
        # terminate all the way so not exicute finally block
        os._exit(1)

    return sigterm_handler


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

    global gromacs
    global mdrun_o
    global mdrun_mpi
    gromacs = args.gromacs
    mdrun_o = args.mdrunopt
    mdrun_mpi = args.mpi

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
    outpath = args.outpath
    localpath = args.outlocal

    # Show options
    LOGGER.info("Running createsims, options are:")
    LOGGER.info(args)
    LOGGER.info("GROMACS       = {}".format(gromacs))
    LOGGER.info("GROMACS single precision path " + shutil.which(gromacs))
    LOGGER.info("mdrun_o       = {}".format(mdrun_o))
    LOGGER.info("mdrun_mpi     = {}".format(mdrun_mpi))
    LOGGER.info("init outpath  = {}".format(outpath))
    LOGGER.info("init outlocal = {}".format(localpath))

    # This is the working dir for the setup
    simpath = outpath if not localpath else localpath
    # Set up the SIGTERM handler
    signal.signal(signal.SIGTERM, sigterm_wrapper(outpath, simpath))

    # @TODO This statement is to kill existing GROMACS runs on the node. This is
    # a temporary fix for flux overscheduling to some ranks.
    utils.sys_call('pkill gmx_mpi*', return_codes=[0,1])

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
            patchSim = creat_cg_patch(patchName[6:].strip(), args.simnum)
        elif (patchName[:8] == "test_raw"):
            LOGGER.info("Testing raw patch '%s'", patchName)
            ppath = "/p/gpfs1/bhatia4/pilot2/campaign1star/patches_37"
            pids = np.array([1062730])
            patch_concentrations, patch_rstates, patch_rpositions = read_patches_from_combined(ppath, pids)
            #patchSim = convert_macro_to_cg_patch("rawpatch_1062730", patch_concentrations[0], patch_rstates[0], patch_rpositions[0])
            patchSim = convert_macro_to_cg_patch("rawpatch_1062730", patch_concentrations[0], patch_rstates, patch_rpositions)
        else:
            error = "createsims was called with -p {} - value not supported".format(patchName)
            LOGGER.error(error)
            raise ValueError(error)

        # If outpath no set (default in MuMMI) get full paths
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
                #iointerface.send_signal(outpath, flag_success)  # no need for this!
                exit()
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
            LOGGER.warning("---> Sim {} folder didn't excist and was created - this is not expected within a MuMMI run".format(outpath))
            utils.sys_call('mkdir -p ' + outpath, test=False)

        LOGGER.info("Creating PG Particle simulation...")
        patchSim.create(simpath) # Creates folder and sub-grid concentration files
        patchSim.save() # Saves patchSim object
        if p != None:
            # This is saved locally (local node) in current woking dir and then coped to gloabl file system at the end
            pfPatch.save_npz("ppatch.npz", p)
        setup_cg_sim(patchSim, simpath)

        LOGGER.info("Running conversion from GROMACS to ddcMD...")
        prime_cg_sim(patchSim, simpath)

        # Setup successfull
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
            LOGGER.error("---> In creatsim finally but din't reach success - most likely a gromacs error/exit.")
            iointerface.send_signal(outpath, flag_failure)
        # Copy back from local and cleanup
        try:
            LOGGER.info("---> Copied form local to global - {}/* {}/".format(simpath, outpath))
            store_files(outpath, simpath, marked_failure)
            cleanup_files(localpath)
        except Exception as error:
            # @TODO - add better error handling and saving of creatsims_error file if always faital
            LOGGER.error("---> ERROR in copying form local to global. {}".format(error))
            # For now don't delete local
            #if localpath != "":
            #    utils.sys_call("rm -r {}".format(simpath))
            iointerface.send_signal(outpath, flag_failure)
            cleanup_files(localpath)
            raise error


if __name__ == "__main__":
    main()
