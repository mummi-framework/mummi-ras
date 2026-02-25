from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.lib import distances
import numpy as np
import time
import math


def findFirstProteinResidue(pdbPatch, firstProteinResname, numberOfMarkedResiduesPerProtein):
    # find the first protein residue
    # NOTES: - The use of .resindex rather than .resid makes this safe even if the topology is a misnumbered .pdb file
    #        - must add error checking in case no THR residue is found or if len(allMarkedResidues) % numberOfMarkedResiduePerProtein != 0
    # allMarkedResidues = pdbPatch.select_atoms("resname %s" % firstProteinResname)
    # if len(allMarkedResidues) > 0:
    #     firstProteinStartResidue = allMarkedResidues[0].resindex
    #     numProteins = int(len(allMarkedResidues)/numberOfMarkedResiduesPerProtein)
    #     return firstProteinStartResidue, numProteins
    # else:
    #     # This happens if there is no protein 
    #     return 0, 0
    return 0, 0


#def getSel_findProtLipidContacts(pdbPatch, numProteins, firstProteinStartResidue, numResPerProteinGTPMg):
#    #ContactProt = list(pdbPatch.select_atoms("resid %d:%d" % (firstProteinStartResidue, firstProteinStartResidue + numProteins*numResPerProteinGTPMg - 1)))
#    ContactProt = []
#    for p in range(numProteins):
#        ContactProt.append(pdbPatch.select_atoms("resid %d:%d" % (firstProteinStartResidue + p*numResPerProteinGTPMg, firstProteinStartResidue + (p+1)*numResPerProteinGTPMg - 1)))
#    plContactLipid = pdbPatch.select_atoms("name C1A D1A T1A R1")
#    return ContactProt, plContactLipid

def getSel_findProtLipidContacts(pdbPatch, numProteins, firstProteinStartResidue, numResPerProteinGTPMg, numResPerProtein, numResGdomain):
    ContactProt = []
    ContactProtBB = []
    ContactG = []  # contactG[] does not contain the GTP or the Mg2+
    ContactHVR = []
    # for p in range(numProteins):
    #     ContactProt.append(pdbPatch.select_atoms("resid %d:%d" % (firstProteinStartResidue + p*numResPerProteinGTPMg, firstProteinStartResidue + (p+1)*numResPerProteinGTPMg - 1)))
    #     ContactProtBB.append(pdbPatch.select_atoms("resid %d:%d and name BB" % (firstProteinStartResidue + p*numResPerProteinGTPMg, firstProteinStartResidue + (p+1)*numResPerProteinGTPMg - 1)))
    #     ContactG.append(pdbPatch.select_atoms("resid %d:%d" % (firstProteinStartResidue + p*numResPerProteinGTPMg, firstProteinStartResidue + p*numResPerProteinGTPMg + numResGdomain - 1)))
    #     ContactHVR.append(pdbPatch.select_atoms("resid %d:%d" % (firstProteinStartResidue + p*numResPerProteinGTPMg + numResGdomain, firstProteinStartResidue + p*numResPerProteinGTPMg + numResPerProtein - 1)))
    plContactLipid = pdbPatch.select_atoms("name C1A D1A T1A R1")
    return ContactProt, ContactG, ContactHVR, ContactProtBB, plContactLipid,


def getSel_allProteinFromAllRas(pdbPatch, numProteins, firstProteinStartResidue, numResPerProteinGTPMg):
    # we will use this in fast selection routines such as findProtLipidContactsQuick()
    # could be incorporated into getSel_findProtLipidContacts() but Chris' previous branch not merged to master yet
    ContactAllProt = pdbPatch.select_atoms("resid %d:%d" % (firstProteinStartResidue, firstProteinStartResidue + (numProteins)*numResPerProteinGTPMg - 1))
    return ContactAllProt


def findProtLipidContactsQuick(pdbPatch, groupContactAllProt, groupPlContactLipid, numProteins, firstProteinStartResidue, numResPerProteinGTPMg, contactcutoff, distanceBackendType, outputFileName):
    # find protein-lipid contacts -- and do it quicker than the original findProtLipidContacts() routine
    # This mechanism could be updated into findProtLipidContacts() but we will not likely use that routine at all anymore so not updating it
    # output file has one row for each protein-lipid contact (file may be empty)
    # format is resA(protA)_lipidType(lipidResID) -- have included record of which protein it is so that state-dependent contacts can later be analyzed if desired
    # WARNING NOTES: - this function fails if the loaded topology was a PDB (not a .tpr) and the residue numbers are not perfectly ordered
    #                  the failure is likely due to the pdbPatch.select_atoms("resid %d:%d" ... line and no way to use "resindex" here
    # plContactSet=set() # this will not double count residue-residue contacts
    # plContactList=[] # this will keep track of multiple events

    # # make subset lists so that we don't need to compute a larger than necessary distance matrix
    # ns_groupCloseLipids = NeighborSearch.AtomNeighborSearch(groupPlContactLipid, box=pdbPatch.dimensions)
    # groupCloseLipids = ns_groupCloseLipids.search(groupContactAllProt,contactcutoff)

    # # compute the distance matrix
    # plContactCmap = distances.distance_array(groupContactAllProt.atoms.positions, groupCloseLipids.atoms.positions, box=pdbPatch.dimensions, result=None, backend=distanceBackendType)

    # # process the distance matrix into useable output
    # foundContacts = np.where(plContactCmap < contactcutoff)

    # for fc in range(len(foundContacts[0])):
    #   tapi = (1 + groupContactAllProt.atoms.resnums[foundContacts[0][fc]] - firstProteinStartResidue) // numResPerProteinGTPMg
    #   plca = 1 + groupContactAllProt.atoms.resnums[foundContacts[0][fc]] - firstProteinStartResidue - tapi*numResPerProteinGTPMg
    #   plcb = 1 + groupCloseLipids.atoms.resnums[foundContacts[1][fc]]
    #   talj = groupCloseLipids.atoms.resnames[foundContacts[1][fc]]
    #   plcelement="%d(%d)_%s(%d)" %(plca, tapi, talj, plcb)
    #   plContactSet.add(plcelement)
    # plContactList += list(plContactSet)
    # return np.array(plContactList)
    return np.array([])


def findProtLipidContacts(pdbPatch, groupContactProt, groupPlContactLipid, numProteins, firstProteinStartResidue, numResPerProteinGTPMg, contactcutoff, distanceBackendType):
    # THIS FUNCTION IS OUTDATED. USE findProtLipidContactsQuick() instead
    # find protein-lipid contacts
    # output file has one row for each protein-lipid contact (file may be empty)
    # format is resA(protA)_lipidType(lipidResID) -- have included record of which protein it is so that state-dependent contacts can later be analyzed if desired
    # WARNING NOTES: - this function fails if the loaded topology was a PDB (not a .tpr) and the residue numbers are not perfectly ordered
    #                  the failure is likely due to the pdbPatch.select_atoms("resid %d:%d" ... line and no way to use "resindex" here
    # plContactSet=set() # this will not double count residue-residue contacts
    # plContactList=[] # this will keep track of multiple events
    # # NOTE: this could be sped up further by doing a single call to distances.distance_array() and processing it differently
    # for tapi in range(numProteins):
    #     # -1 below means I am evaluating GTP and Mg2+ too
    #     plContactSet.clear()  # new set for each case
    #     plContactCmap = distances.distance_array(groupContactProt[tapi].atoms.positions, groupPlContactLipid.residues.atoms.positions, box=pdbPatch.dimensions, result=None, backend=distanceBackendType)
    #     foundContacts = np.where(plContactCmap < contactcutoff)
    #     for fc in range(len(foundContacts[0])):
    #       plca = 1 + groupContactProt[tapi].atoms.resnums[foundContacts[0][fc]] - firstProteinStartResidue - tapi*numResPerProteinGTPMg
    #       plcb = 1 + groupPlContactLipid.residues.atoms.resnums[foundContacts[1][fc]]
    #       talj =  groupPlContactLipid.residues.atoms.resnames[foundContacts[1][fc]]
    #       plcelement="%d(%d)_%s(%d)" %(plca, tapi, talj, plcb)
    #       plContactSet.add(plcelement)
    #     # want to count duplicates for different pairs of proteins (would get it anyway now that I am numbering the protein, but leave this for robustness
    #     plContactList += list(plContactSet)

    # #outputFileHandle = open(outputFileName, "w")
    # #for vari in plContactList:
    # #    outputFileHandle.write("%s\n" % vari)
    # #outputFileHandle.close()
    # return np.array(plContactList)
    return np.array([])


def findProtProtContacts(pdbPatch, groupContactProt, numProteins, firstProteinStartResidue, numResPerProteinGTPMg, contactcutoff, distanceBackendType):
    # find protein-protein intermolecular contacts
    # output file has one row for each intermolecular protein-protein residue-residue contact (file may be empty)
    # format is resA(protA)_resB(protB) -- have included record of which protein it is so that state-dependent contacts can later be analyzed if desired
    # WARNING NOTES: - this function fails if the loaded topology was a PDB (not a .tpr) and the residue numbers are not perfectly ordered
    #                  the failure is likely due to the pdbPatch.select_atoms("resid %d:%d" ... lines and no way to use "resindex" here
    # ppContactSet=set() # this will not double count residue-residue contacts (unless they are present in both combinations like symmetry)
    # ppContactList=[] # this will keep track of multiple events
    # # NOTE: this could be sped up further by doing a single call to distances.distance_array() and processing it differently
    # for tapi in range(numProteins):
    #     # -1 below means I am evaluating GTP and Mg2+ too
    #     for tapj in range(tapi+1, numProteins):
    #         ppContactSet.clear()  # new set for each case
    #         ppContactCmap = distances.distance_array(groupContactProt[tapi].atoms.positions, groupContactProt[tapj].atoms.positions, box=pdbPatch.dimensions, result=None, backend=distanceBackendType)
    #         foundContacts = np.where(ppContactCmap < contactcutoff)
    #         for fc in range(len(foundContacts[0])):
    #             ppca = 1 + groupContactProt[tapi].atoms.resnums[foundContacts[0][fc]] - firstProteinStartResidue - tapi*numResPerProteinGTPMg
    #             ppcb = 1 + groupContactProt[tapj].atoms.resnums[foundContacts[1][fc]] - firstProteinStartResidue - tapj*numResPerProteinGTPMg
    #             ppcelement="%d(%d)_%d(%d)" %(ppca, tapi, ppcb, tapj)
    #             ppContactSet.add(ppcelement)
    #         # want to count duplicates for different pairs of proteins (would get it anyway now that I am numbering the protein, but leave this for robustness
    #         ppContactList += list(ppContactSet)
    # # output
    # #outputFileHandle = open(outputFileName, "w")
    # #for vari in ppContactList:
    # #    outputFileHandle.write("%s\n" % vari)
    # #outputFileHandle.close()
    # return np.array(ppContactList)
    return np.array([])


def findGHVRContacts(pdbPatch, groupContactG, groupContactHVR, numProteins, firstProteinStartResidue, numResPerProteinGTPMg, contactcutoff, distanceBackendType):
    # find G domain - HVR intermolecular contacts
    # output file has one row for each intramolecular protein-protein residue-residue contact (file may be empty)
    # format is resA(protA)_resB(protB) -- have included record of which protein it is so that state-dependent contacts can later be analyzed if desired
    # WARNING NOTES: - this function fails if the loaded topology was a PDB (not a .tpr) and the residue numbers are not perfectly ordered
    #                  the failure is likely due to the pdbPatch.select_atoms("resid %d:%d" ... lines and no way to use "resindex" here
    # ppContactSet=set() # this will not double count residue-residue contacts (unless they are present in both combinations like symmetry)
    # ppContactList=[] # this will keep track of multiple events
    # # NOTE: this could be sped up further by doing a single call to distances.distance_array() and processing it differently
    # for tapi in range(numProteins):
    #     ppContactSet.clear()  # new set for each case
    #     ppContactCmap = distances.distance_array(groupContactG[tapi].atoms.positions, groupContactHVR[tapi].atoms.positions, box=pdbPatch.dimensions, result=None, backend=distanceBackendType)
    #     foundContacts = np.where(ppContactCmap < contactcutoff)
    #     for fc in range(len(foundContacts[0])):
    #         ppca = 1 + groupContactG[tapi].atoms.resnums[foundContacts[0][fc]] - firstProteinStartResidue - tapi*numResPerProteinGTPMg
    #         ppcb = 1 + groupContactHVR[tapi].atoms.resnums[foundContacts[1][fc]] - firstProteinStartResidue - tapi*numResPerProteinGTPMg
    #         ppcelement="%d(%d)_%d(%d)" %(ppca, tapi, ppcb, tapi)
    #         ppContactSet.add(ppcelement)
    #     ppContactList += list(ppContactSet)
    # return np.array(ppContactList)
    return np.array([])


def getSel_findProtRigidBB(pdbPatch, numProteins, firstProteinStartResidue, numResPerProteinGTPMg):
    #Ras starts at THR2, so use: 1-25, 39-55, 68-165
    ContactProt = []
    # for p in range(numProteins):
    #     proteinStartResidue = firstProteinStartResidue + p*numResPerProteinGTPMg
    #     ContactProt.append(pdbPatch.select_atoms("resid %d:%d or resid %d:%d or resid %d:%d and name BB" % (proteinStartResidue, proteinStartResidue + 24,proteinStartResidue + 38, proteinStartResidue + 54,proteinStartResidue + 67, proteinStartResidue + 164)))
    return ContactProt


def findGdomainRMSD(groupProtBB, groupRefBB, numProteins):
    # find RMSD of G domain to reference structure
    # output file has one column per Ras, with the RMSD to the reference for selected residues (units=A)
    # WARNING NOTES: - this function fails if the loaded topology was a PDB (not a .tpr) and the residue numbers are not perfectly ordered
    #                  the failure is likely due to the pdbPatch.select_atoms("resid %d:%d" ... lines and no way to use "resindex" here
    #outputFileHandle = open(outputFileName, "w")
    outList = []
    # for p in range (numProteins):
    #     R = RMSD(groupProtBB[p], groupRefBB[0], select="name BB", filename="null.not.written", weights=None, tol_mass=999.99)
    #     R.run()
    #     outList.append(R.rmsd[0][2])
        #outputFileHandle.write("%f\t" % R.rmsd[0][2])
    #outputFileHandle.write("\n")
    #outputFileHandle.close()
    return np.array(outList)


def findGdomainRMSDsurrogate(groupProtBB, groupRefBB, numProteins, distanceBackendType):
    # The RMSD calculation takes 1.6 seconds, so use this as something a lot faster but still a sanity check on the structure
    #outputFileHandle = open(outputFileName, "w")
    outList = []
    # refSelfDistMap = distances.self_distance_array(groupRefBB[0].atoms.positions, box=groupRefBB[0].dimensions, result=None, backend=distanceBackendType)
    # nelem =  len(refSelfDistMap)
    # for p in range (numProteins):
    #     thisSelfDistMap = distances.self_distance_array(groupProtBB[p].atoms.positions, box=groupProtBB[p].dimensions, result=None, backend=distanceBackendType)
    #     euclidDist = np.linalg.norm(thisSelfDistMap - refSelfDistMap)
    #     norm = math.sqrt((euclidDist*euclidDist) / nelem)
    #     outList.append(norm)
        #outputFileHandle.write("%f\t" % norm)
    #        for rli in range(len(thisSelfDistMap)):
    #            outputFileHandle.write("%d\t%d\t%f\t%f\n" % (p, rli, refSelfDistMap[rli], thisSelfDistMap[rli]))
    #outputFileHandle.write("\n")
    #outputFileHandle.close()
    return np.array(outList)


def findGdomainRgyr(groupProtBB, numProteins):
    # find Rgyr of G domain (residues 1-165)
    # output file has one column per Ras, with the Rg of residue 1-165 BB beads
    # WARNING NOTES: - this function fails if the loaded topology was a PDB (not a .tpr) and the residue numbers are not perfectly ordered
    #                  the failure is likely due to the pdbPatch.select_atoms("resid %d:%d" ... lines and no way to use "resindex" here
    #outputFileHandle = open(outputFileName, "w")
    outList = []
    # for p in range (numProteins):
    #     thisRgyr=groupProtBB[p].radius_of_gyration()
    #     outList.append(thisRgyr)
        #outputFileHandle.write("%f\t" % thisRgyr)
    #outputFileHandle.write("\n")
    #outputFileHandle.close()
    return np.array(outList)


def getSel_findGdomainBB(pdbPatch, numProteins, firstProteinStartResidue, numResPerProteinGTPMg, numResGdomain):
    #ContactProt = list(pdbPatch.select_atoms("resid %d:%d" % (firstProteinStartResidue, firstProteinStartResidue + numProteins*numResPerProteinGTPMg - 1)))
    ContactProt = []
    # for p in range(numProteins):
    #     proteinStartResidue = firstProteinStartResidue + p*numResPerProteinGTPMg
    #     ContactProt.append(pdbPatch.select_atoms("resid %d:%d and name BB" % (proteinStartResidue, proteinStartResidue + numResGdomain - 1)))
    return ContactProt


def findGdomainCOM(groupProtBB, numProteins):
    # find center of mass of G domain X,Y,Z (res Thr1 to His165)
    # output file has three columns per Ras, with the X,Y,Z com of residue 1-165 BB beads
    # WARNING NOTES: - this function fails if the loaded topology was a PDB (not a .tpr) and the residue numbers are not perfectly ordered
    #                  the failure is likely due to the pdbPatch.select_atoms("resid %d:%d" ... lines and no way to use "resindex" here
    #outputFileHandle = open(outputFileName, "w")
    outList = []
    # for p in range (numProteins):
    #     thisCOM=groupProtBB[p].center_of_mass()
    #     outList.append(thisCOM)
        #outputFileHandle.write("%f\t%f\t%f\t" % (thisCOM[0], thisCOM[1], thisCOM[2]))
    #outputFileHandle.write("\n")
    #outputFileHandle.close()
    return np.array(outList)


def getSel_findLipidLeaflet(pdbPatch):
    # This list of lipid types and special beads is hard-coded here
    identList = 'POPX_PO4_C4A_C4B PAPC_PO4_C5A_C4B POPE_PO4_C4A_C4B DIPE_PO4_C4A_C4B DPSM_PO4_C3A_C4B PAPS_PO4_C5A_C4B PAP6_PO4_C5A_C4B POPC_PO4_C4A_C4B'

    groupLipidHeadSelection = []
    groupLipidTailASelection = []
    groupLipidTailBSelection = []

    # for thisToken in identList.split():
    #     lipidIdent,headIdent,tail1Ident,tail2Ident = thisToken.split("_")
    #     groupLipidHeadSelection.append(pdbPatch.select_atoms("resname %s and name %s" % (lipidIdent, headIdent)))
    #     groupLipidTailASelection.append(pdbPatch.select_atoms("resname %s and name %s" % (lipidIdent, tail1Ident)))
    #     groupLipidTailBSelection.append(pdbPatch.select_atoms("resname %s and name %s" % (lipidIdent, tail2Ident)))

    return groupLipidHeadSelection, groupLipidTailASelection, groupLipidTailBSelection


def findLipidLeaflet(groupLipidHeadSelection, groupLipidTailASelection, groupLipidTailBSelection):
    # find out what leaflet each lipid is in. This is not expected to be exact frame-by-frame, but a smoothed version should show if there are flip-flops
    # output file has one column for each lipid with the z-displacement between head and tail
    # orientation vector is defined so that positive means the lipid is in the upper leaflet
    # NOTES: - This doesn't select by residue number so it is safe even if the topology is a misnumbered .pdb file

    #outputFileHandle = open(outputFileName, "w")
    outArray = np.empty(len(groupLipidHeadSelection), dtype=object) # @WARNING will have variable size objects inside so not a matrix
    # for nsel in range(len(groupLipidHeadSelection)):
    #     cArray = np.zeros(len(groupLipidHeadSelection[nsel]))
    #     for lipi in range(len(groupLipidHeadSelection[nsel])):
    #         zdist = groupLipidHeadSelection[nsel][lipi].position[2] - (groupLipidTailASelection[nsel][lipi].position[2] + groupLipidTailBSelection[nsel][lipi].position[2])/2.0
    #         cArray[lipi] = zdist
    #         #outputFileHandle.write("%f\t" % zdist)
    #     outArray[nsel] = cArray
    #outputFileHandle.write("\n")
    #outputFileHandle.close()
    return outArray


def getSel_writeCompressedSelection(pdbPatch):
    groupCompressedTrajSelection = pdbPatch.select_atoms("resname GTP or name BB SC1 SC2 SC3 F1 F2 F3 F4 MG+ C1A D1A T1A R1")
    return groupCompressedTrajSelection


def writeCompressedSelection(groupCompressedTrajSelection, outputFileName):
    # write compressed trajectory with coordinates of one atom per lipid and BB atoms of protein
    # BB is protein and F1 is farnesyl tail bead for analysis (but also keep F2, F3, F4 for visualization)
    # C1A bead is in DIPE, POPC, POPE, POPX; D1A bead is in PAP6, PAPC, PAPS; T1A bead is in DPSM; R1 bead is in CHOL
    # MG+ bead is magnesium; resname GTP
    # NOTES: - This doesn't select by residue number so it is safe even if the topology is a misnumbered .pdb file
    groupCompressedTrajSelection.write(outputFileName)


def writeBoxDimensions(pdbPatch):
    # write box dimensions (pdbPatch.dimensions has the same values as cg.syst.dimensions)
    # NOTES: - This doesn't select by residue number so it is safe even if the topology is a misnumbered .pdb file
    #outputFileHandle = open(outputFileName, "w")
    #outputFileHandle.write("%f\t%f\t%f\n" % (pdbPatch.dimensions[0], pdbPatch.dimensions[1], pdbPatch.dimensions[2]))
    #outputFileHandle.close()
    return pdbPatch.dimensions[0:3]


def writeNumRas(numKRAS):
    # write the number of Ras
    # NOTES: - This doesn't select by residue number so it is safe even if the topology is a misnumbered .pdb file
    #outputFileHandle = open(outputFileName, "w")
    #outputFileHandle.write("%f\h" % numKRAS)
    #outputFileHandle.close()
    return np.array([numKRAS])
