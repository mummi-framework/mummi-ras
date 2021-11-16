###############################################################################
# @todo add Pilot2-splash-app disclaimer
# @authors  Helgi I. Ingolfsson <ingolfsson1@llnl.gov>
###############################################################################

import sys
import math
import os
import random
import numpy as np
import numpy.linalg
import MDAnalysis as mda

from mummi_core.utils import Naming
from logging import getLogger
LOGGER = getLogger(__name__)

def kras_Rotate_Shift(krasFile, krasPos, angle, krasType):
    """For x1 KRAS rotate by a specific angle around z axes and shit to given x,y pos
       Currently supports two kRAS 4A and 4B - if 4A in type we rotate around COM of CYF and CYP
       else we assume it's kRAS 4B and use CYF only.
    """

    kras = mda.Universe(krasFile)
    com = []
    if "4A" in krasType:
        # Assuming protein is KRAS 4A and center/roate about COM Farnistalation and Palmitilization
        com = kras.select_atoms("(resname CYF and name SC1) or (resname CYP and name SC1)").center_of_mass()
    else:
        # Assuming protein is KRAS 4B and center/roate about Farnistalation
        com = kras.select_atoms("resname CYF and name SC1").center_of_mass()
    com_z = 5.455 # Set all KRAS FAR group at z = 5.455 (in the default simulations this is a little below the bilayer)
    com[2] = com[2] - com_z
    kras.atoms.translate(-1 * com) # Center x,y based on Farnistalation
    #angle = random.randint(1, 360)
    kras.atoms.rotateby(angle, [0,0,1], point=[0,0,com_z])
    shiftBy = np.zeros(3)
    shiftBy[0:2] = 10 * np.array(krasPos)
    kras.atoms.translate(shiftBy)
    LOGGER.debug("Rotate {}, shift {}".format(angle, shiftBy))

    return kras

def kras_Rotate_Shift_old(krasFile, krasPos, angle):
    """For x1 KRAS rotate by a specific angle around z axes and shit to given x,y pos - Works for 4B only CYF based"""

    kras = mda.Universe(krasFile)
    farCOM = kras.select_atoms("resname CYF and name SC1").center_of_mass()
    farCOM_z = 5.455 # Set all KRAS FAR group at z = 5.455 (in the default simulations this is a little below the bilayer)
    farCOM[2] = farCOM[2] - farCOM_z
    kras.atoms.translate(-1 * farCOM) # Center x,y based on Farnistalation
    #angle = random.randint(1, 360)
    kras.atoms.rotateby(angle, [0,0,1], point=[0,0,farCOM_z])
    shiftBy = np.zeros(3)
    shiftBy[0:2] = 10 * np.array(krasPos)
    kras.atoms.translate(shiftBy)
    LOGGER.debug("Rotate {}, shift {}".format(angle, shiftBy))

    return kras

def protein_align(proteinUnits):
    """ For a list protein molecules places them together.
            RAS-ONLY randomly rotate around z, palce farnesyl at x,y
            RAS-RAF roate to a specified angle, palce farnesyl at x,y
        check if overlapp and save down roated structures,
        if no solution repeat with all now randoms
        if those fail max times shift proteins out

        If no placment found - move new protein out by "incShiftprotien"
        If "maxShiftprotein" reach with no placement throw error
    """

    # If no protein return
    if len(proteinUnits) == 0:
        return []

    # $TODO check for protein type - if RAS-RAF use angle in file
    maxTryPerprotein = 50 # Max tries for diffrent random angles
    incShiftprotein = 0.5 # Dist to shift if no placment found in nm
    maxShiftprotein = 4   # Max distance to try before throwing error in nm
    maxStructureTry = 10  # Max number of times new structures are sampled per shift interval

    currentShiftprotein = 0
    currentStructureTryCount = 0
    currentAngleTryCount = 0

    proteinNewNames = []
    proteinPos = []

    # loop over shift of all proteins "away" from other proteins
    while proteinNewNames == [] and currentShiftprotein <= maxShiftprotein:
        currentStructureTryCount = 0

        # loop over diffrent new input structurs
        while proteinNewNames == [] and currentStructureTryCount <= maxStructureTry:
            currentAngleTryCount = 0

            # Get random new input structures for all proteins - changes protein dictionaries in palce
            getRandomInputStructurs(proteinUnits)

            # loop over diffrent random inital angle tries
            while proteinNewNames == [] and currentAngleTryCount < maxTryPerprotein:
                currentAngleTryCount += 1
                LOGGER.info("--> protein_align angle_t {} of max {}, struct_t {} of max {}, cshift {} of max {}".format(currentAngleTryCount, maxTryPerprotein, currentStructureTryCount, maxShiftprotein, currentShiftprotein, maxShiftprotein))
                proteinNewNames, proteinPos = protein_align_one_try(proteinUnits, currentShiftprotein)

            currentStructureTryCount += 1

        currentShiftprotein += incShiftprotein

    # No structure found
    if proteinNewNames == []:
        errorMessage = "Failed to insert protein molecules - max shift {} reached for all angle tries {}.".format(currentShiftprotein, maxTryPerprotein)
        LOGGER.error(errorMessage)
        raise ValueError(errorMessage)

    return proteinNewNames, proteinPos


def getRandomInputStructurs(proteinUnits):
    """Assign random input structures for all protein units form library"""

    # Lookup table to support diffrent library sizes for diffrent structure/state input
    # @WARNING library format has to be exact - see naming
    # @WARNING these 4A numbers are for C4 healing run - update for C4
    strictire_library_size = { "RAS-ONLY":  {  "a":5000 ,  "b":5000,  "c":5000 } \
                             , "RAS-RAF":   { "ma":5000 , "mb":5000, "za":5000 } \
                             , "RAS4A-ONLY":{  "a":5000 ,  "b":5000,  "c":5000 } \
                             , "RAS4A-RAF": { "ma":5000 , "mb":5000, "za":1400 }  }

    for i in range(len(proteinUnits)):
        cProteinDict = proteinUnits[i]
        if ("structure_pre_assigned" in cProteinDict and cProteinDict["structure_pre_assigned"] == True):
            # Structure has already been assigned don't do anything use assigned structure.
            LOGGER.debug("Protein id {} structure alreaddy asign to {} use that one".format(cProteinDict["id"], cProteinDict["structure"]))
        else:
            # This is for both RAS-only (4B) and RAS-RAF (4B) and all states are equal
            # Structure not set get a random configuration
            # Numer of diffrent configurations saved for each RAS-ONLY states as well as RAS-RAF states (not all need to thave the same count)
            # @WARNING hardcoded  _0, _1, ... _librarySize-1
            librarySize = strictire_library_size[cProteinDict["type"]][cProteinDict["state"]]
            cInt = random.randint(0, librarySize - 1)
            cStructureName = "{}/{}".format(Naming.dir_res('structures'), Naming.protein_structure(cProteinDict["type"], cProteinDict["state"], cInt))
            LOGGER.debug("Get {} structure: state {}, randint {}, filename {}, librarySize {}".format(cProteinDict["type"], cProteinDict["state"], cInt, cStructureName, librarySize))
            if not os.path.isfile(cStructureName):
                error = "---> Protein structure file {} not found.".format(cStructureName)
                LOGGER.error(error)
                raise ValueError(error)
            cProteinDict["structure"] = cStructureName


def getCurrentAngle(proteinUnit):
    """Get random 0-360 angle except if RAS-RAF (not z sate) then return angel that will give macro model angle alignment for current structure with a small random addition"""

    if (proteinUnit["type"] == "RAS4A-RAF" or proteinUnit["type"] == "RAS-RAF") and proteinUnit["state"] != "za":
        macroAngle = proteinUnit["angle"]
        ras_raf = mda.Universe(proteinUnit["structure"])
        farCOM = ras_raf.select_atoms("resname CYF and name SC1").center_of_mass()
        # from get tilt z state mulitple - and now (resname ACE1 and name BB) is always = 1 - first residue in the file
        start_of_g = 1
        zshifting=-2
        crdCOM = ras_raf.select_atoms('resid '+str(start_of_g+zshifting+291)+':'+str(start_of_g+zshifting+294)+' ' +\
                                             str(start_of_g+zshifting+278)+':'+str(start_of_g+zshifting+281)+' ' +\
                                             ' and (name CA or name BB)').center_of_mass()
        #crdCOM = ras_raf.select_atoms('resid '+str(start_of_g+291)+':'+str(start_of_g+294)+' ' +\
        #                                       str(start_of_g+278)+':'+str(start_of_g+281)+' ' +\
        #                                      ' and (name CA or name BB)').center_of_mass()
        deltaCOM = 0.1 * (crdCOM - farCOM) # Convert from A to nm
        cgAngle = math.degrees(math.atan2(deltaCOM[1], deltaCOM[0])) % 360

        # Then do a small random displacment - +15 - -15 (addjust this angle)
        cgRandom = (random.random()*30) - 15.0
        cgAngle += cgRandom
        LOGGER.debug("Get RAS-RAF new angle - farCOM {}, crdCOM {}, deltaCOM {} / nm, angle {}, macroAngle {}, deltaAngle {}, cgRandom {}".format(farCOM, crdCOM, deltaCOM, cgAngle, macroAngle, macroAngle - cgAngle, cgRandom))

        return macroAngle - cgAngle
    else:
        return random.random()*360


def protein_align_one_try(proteinUnits, currentShiftprotein):
    """Try alignment at specific shift - return [], [] if unsucessfull"""

    # Place first protein with rotation if needed
    # Iterativly go through other protein molecules
    # Place protein - rotate if needed and check if overlapp with placed atom within "withinCutoff"
    # If overlapp return fail ([], [])

    proteinNewNames = []
    proteinPos = []

    # First
    angle = getCurrentAngle(proteinUnits[0])
    protein = kras_Rotate_Shift(proteinUnits[0]["structure"], proteinUnits[0]["position"], angle, proteinUnits[0]["type"])
    proteinName = "protein-{}-rotate.gro".format(proteinUnits[0]["id"])
    protein.atoms.write(proteinName) # Save file to dist
    proteinNewNames.append(proteinName)
    proteinPos.append(proteinUnits[0]["position"])

    # All others
    protein_len = len(protein.atoms)
    allprotein = protein
    withinCutoff = 8  # This is in Armstrong check values 5 - 10

    for i in range(1, len(proteinUnits)):
        cProtein = []
        j = 0
        # @TODO clean print stuff
        avgPost = np.average(proteinPos, axis=0) # Average protein posision up to now
        #print(avgPost)
        #print(proteinUnits[i]["position"])
        cVect = proteinUnits[i]["position"] - avgPost
        #print(cVect)
        awayUnitVector = cVect / np.linalg.norm(cVect)
        #print(awayUnitVector)

        '''
        while(True):
            cProteinPos = proteinUnits[i]["position"] + currentShiftprotein * awayUnitVector
            angle = random.randint(1, 360)
            cProtein = kras_Rotate_Shift(proteinUnits[i]["structure"], cProteinPos, angle)
            sim_merge = mda.Merge(cProtein.atoms, allprotein.atoms)
            stringSel = 'bynum {}:{} and around {} (bynum {}:{})'.format(1, protein_len, withinCutoff, len(cProtein.atoms)+1, len(sim_merge.atoms))

            if len(sim_merge.select_atoms(stringSel)) == 0:
                allprotein = sim_merge
                proteinName = "protein-{}-rotate.gro".format(proteinUnits[i]["id"])
                cProtein.atoms.write(proteinName)
                proteinNewNames.append(proteinName)
                proteinPos.append(cProteinPos)  # Note protein posision might have been sifted from orginal
                break
            elif (j > maxTryPerprotein):
                if (currentShiftprotein > maxShiftprotein):
                    LOGGER.error("Failed to insert protein molecules - max shift reached for all angle tries.")
                    sys.exit(-1)
                else:
                    currentShiftprotein = currentShiftprotein + incShiftprotein
                    j = 0
                    LOGGER.info("--> Max try for roated reached, now shift pos by {} and try again".format(currentShiftprotein))
            else:
                LOGGER.info("--> Protein overlapped for protein# {} with rest, angle try #{} of {} with extra shift of {} - retry".format(i, j, maxTryPerprotein, currentShiftprotein))
            j = j + 1
        '''
        cProteinPos = proteinUnits[i]["position"] + currentShiftprotein * awayUnitVector
        angle = getCurrentAngle(proteinUnits[i])
        cProtein = kras_Rotate_Shift(proteinUnits[i]["structure"], cProteinPos, angle, proteinUnits[i]["type"])
        sim_merge = mda.Merge(cProtein.atoms, allprotein.atoms)
        stringSel = 'bynum {}:{} and around {} (bynum {}:{})'.format(1, protein_len, withinCutoff, len(cProtein.atoms)+1, len(sim_merge.atoms))

        if len(sim_merge.select_atoms(stringSel)) == 0:
            allprotein = sim_merge
            proteinName = "protein-{}-rotate.gro".format(proteinUnits[i]["id"])
            cProtein.atoms.write(proteinName)
            proteinNewNames.append(proteinName)
            proteinPos.append(cProteinPos)  # Note protein posision might have been sifted from orginal
        else:
            LOGGER.debug("--> Protein {} overlapped with rest, with pos {}, with angle {}, and extra shift {}".format(i, cProteinPos, angle, currentShiftprotein))
            return [], []

    # Test location of all proteins
    allprotein.atoms.write("protein-ALL.gro")

    return proteinNewNames, proteinPos

'''
# Test with a number of imputs - assuming this is run in
#  pilot2-splash-app/pysplash/pysplash
#  python createsim/placeproteins.py
if __name__ == "__main__":

    # Get array of KRAS structures and array of KRAS posision - test sparsh x4
    krasStruct = "../../resources/KRAS-structures/KRAS-04-protein-cg-M22-em-shift.gro"
    krasStructs = [krasStruct, krasStruct, krasStruct, krasStruct]
    krasPos = [[3, 3], [3, 16], [20, 3], [20, 16]]
    krasNewName = kras_align(krasStructs, krasPos)
    print(krasNewName)

    # Get array of KRAS structures and array of KRAS posision - test dense x8
    krasStructs = [krasStruct, krasStruct, krasStruct, krasStruct, krasStruct, krasStruct, krasStruct, krasStruct]
    krasPos = [[2, 2], [2, 8], [8, 2], [5, 16], [20, 3], [20, 16], [10, 20], [20, 20]]
    krasNewName = kras_align(krasStructs, krasPos)
    print(krasNewName)
'''
