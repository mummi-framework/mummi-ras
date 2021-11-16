###############################################################################
# @todo add Pilot2-splash-app disclaimer
###############################################################################

""" Get's KRAS states """

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.lib.mdamath import make_whole
import os
import numpy as np
import math

############## Below section needs to be uncommented ############
import mummi_core
import mummi_ras
from mummi_core.utils import Naming

# Logger has to be initialized the first thing in the script
from logging import getLogger
LOGGER = getLogger(__name__)

# Innitilize MuMMI if it has not been done before
MUMMI_ROOT = mummi_core.init(True)

# This is needed so the Naming works below
#@TODO fix this so we don't have these on import make them as an init
# mummi.init()

dirKRASStates = Naming.dir_res('states')
dirKRASStructures = Naming.dir_res('structures')
#RAS_ONLY_macrostate = np.loadtxt(os.path.join(dirKRASStates, "RAS-ONLY.microstates.txt"))
RAS_ONLY_macrostate = np.loadtxt(os.path.join(dirKRASStates, "ras-states.txt"),comments='#')
#RAS_RAF_macrostate = np.loadtxt(os.path.join(dirKRASStates, "RAS-RAF.microstates.txt"))
RAS_RAF_macrostate = np.loadtxt(os.path.join(dirKRASStates, "ras-raf-states.txt"),comments='#')  # Note diffrent number of columns so index change below
RAS4A_ONLY_macrostate = np.loadtxt(os.path.join(dirKRASStates, "ras4a-states.txt"),comments='#') # New RAS4A states created by Tomas from C4 healing run
# @WARNING RAS4A_RAF_macrostate not defined, in C4 we use RAS4B-RAF macro states also for RAS4A-RAF (distribtion from healing runs looked very similar)
############## above section needs to be uncommented ############

kras_ref_universe = mda.Universe(os.path.join(dirKRASStructures, "RAS-ONLY-reference-structure.gro"))

######### Below hard codes the number of residues within RAS-only and RAS-RAF ##########
RAS_only_num_res = 184               ###Updated for C4###
RAS_RAF_num_res = 320                ###Updated for C4###
RAS4A_only_num_res = 185               ###Updated for C4###
RAS4A_RAF_num_res = 321                ###Updated for C4###
######### Above hard codes the number of residues within RAS-only and RAS-RAF ##########

####### This can be removed
# def get_kras(syst, kras_start):
#     """Gets all atoms for a KRAS protein starting at 'kras_start'."""
#     return syst.atoms[kras_start:kras_start+428]
####### This can be removed

def get_segids(u):
    """Identifies the list of segments within the system. Only needs to be called x1 time"""
    segs = u.segments
    segs = segs.segids
    ras_segids = []
    rasraf_segids = []
    ras4a_segids = []
    ras4araf_segids = []
    for i in range(len(segs)):
    #     print(segs[i])
        val = -1
        for f in range(0, 2):
            val = segs[i].find('_', val + 1)
        if segs[i][val+1:] == 'RAS':               ###Updated for C4###
            ras_segids.append(segs[i])
        if segs[i][val+1:] == 'RAS_RAF':               ###Updated for C4###
            rasraf_segids.append(segs[i])
        if segs[i][val+1:] == 'RAS4A':               ###Updated for C4###
            ras4a_segids.append(segs[i])
        if segs[i][val+1:] == 'RAS4A_RAF':               ###Updated for C4###
            ras4araf_segids.append(segs[i])
    return ras_segids, rasraf_segids, ras4a_segids, ras4araf_segids              ###Updated for C4###

def get_protein_info(u,tag):
    """Uses the segments identified in get_segids to make a list of all proteins in the systems.\
    Outputs a list of the first residue number of the protein, and whether it is 'RAS-ONLY', or 'RAS-RAF'.\
    The 'tag' input defines what is used to identify the first residue of the protein. i.e. 'resname ACE1 and name BB'.\
    Only needs to be called x1 time"""
    ras_segids, rasraf_segids, ras4a_segids, ras4araf_segids = get_segids(u)

    if len(ras_segids) > 0:
        RAS = u.select_atoms('segid '+ras_segids[0]+' and '+str(tag))
    else:
        RAS = []

    if len(rasraf_segids) > 0:
        RAF = u.select_atoms('segid '+rasraf_segids[0]+' and '+str(tag))
    else:
        RAF = []

    if len(ras4a_segids) > 0:
        RAS4A = u.select_atoms('segid '+ras4a_segids[0]+' and '+str(tag))       ###Updated for C4###
    else:
        RAS4A = []

    if len(ras4araf_segids) > 0:
        RAF4A = u.select_atoms('segid '+ras4araf_segids[0]+' and '+str(tag))        ###Updated for C4###
    else:
        RAF4A = []

    protein_info = []#np.empty([len(RAS)+len(RAF),2])
    for i in range(len(RAS)):
        protein_info.append((RAS[i].resid,'RAS-ONLY'))              ###Updated for C4###
    for i in range(len(RAF)):
        protein_info.append((RAF[i].resid,'RAS-RAF'))               ###Updated for C4###
    for i in range(len(RAS4A)):
        protein_info.append((RAS4A[i].resid,'RAS4A-ONLY'))              ###Updated for C4###
    for i in range(len(RAF4A)):
        protein_info.append((RAF4A[i].resid,'RAS4A-RAF'))               ###Updated for C4###

    ######## sort protein info
    protein_info = sorted(protein_info)
    ######## sort protein info
    return protein_info

def get_ref_kras():
    """Gets the reference KRAS struct. Only called x1 time when class is loaded"""
    start_of_g_ref = kras_ref_universe.residues[0].resid
    ref_selection = 'resid '+str(start_of_g_ref)+':'+str(start_of_g_ref+24)+' ' +\
                             str(start_of_g_ref+38)+':'+str(start_of_g_ref+54)+' ' +\
                             str(start_of_g_ref+67)+':'+str(start_of_g_ref+164)+' ' +\
                             'and (name CA or name BB)'

    r2_26r40_56r69_166_ref = kras_ref_universe.select_atoms(str(ref_selection))
    return kras_ref_universe.select_atoms(str(ref_selection)).positions - kras_ref_universe.select_atoms(str(ref_selection)).center_of_mass()


# Load inital ref frames (only need to do this once)
ref0 = get_ref_kras()

def getKRASstates(u,kras_indices):
    """Gets states for all KRAS proteins in path."""

#     res_shift = 8

#     all_glycine = u.select_atoms("resname GLY")

#     kras_indices = []
#     for i in range(0, len(all_glycine), 26):
#         kras_indices.append(all_glycine[i].index)

########## Below is taken out of the function so it is only done once #########
#     kras_indices = get_protein_info(u,'resname ACE1 and name BB')
########## Above is taken out of the function so it is only done once #########


    ALLOUT = []
    for k in range(len(kras_indices)):
        start_of_g = kras_indices[k][0]
        protein_type = str(kras_indices[k][1])
########## BELOW SECTION TO DETERMINE WHICH RESIDUES ARE PART OF THE PROTEIN GROUP - NEEDED FOR PBC REMOVAL ##############
########## POTENTIALLY REDO WITH A 'HARD-CODED' NUMBER OF RESIDUES PER PROTEIN GROUP (WHETHER RAS-ONLY OR RAS-RAF) #######
########## HAS BEEN REDONE WITH A 'HARD-CODED' NUMBER OF RESIDUES PER PROTEIN GROUP (WHETHER RAS-ONLY OR RAS-RAF) ########
#         if len(kras_indices) == 1:
#             krases0_BB = u.select_atoms('resid '+str(start_of_g)+':'+str(len(u.residues))+' and name BB') ####### HAS TO BE FIXED FOR BACKBONE ATOMS FOR SPECIFIC PROTEIN
#         elif len(kras_indices) > 1:
#             if k == len(kras_indices)-1:
#                 krases0_BB = u.select_atoms('resid '+str(start_of_g)+':'+str(len(u.residues))+' and name BB')
#             else:
#                 krases0_BB = u.select_atoms('resid '+str(start_of_g)+':'+str(kras_indices[k+1][0])+' and name BB')
########## ABOVE SECTION TO DETERMINE WHICH RESIDUES ARE PART OF THE PROTEIN GROUP - NEEDED FOR PBC REMOVAL ##############

########## Below hard codes the number of residues/beads in the RAS-ONLY and RAS-RAF simulations #########################
        if protein_type == 'RAS-ONLY':          ###Updated for C4###
            num_res = RAS_only_num_res
        elif protein_type == 'RAS-RAF':         ###Updated for C4###
            num_res = RAS_RAF_num_res
        elif protein_type == 'RAS4A-ONLY':          ###Updated for C4###
            num_res = RAS4A_only_num_res
        elif protein_type == 'RAS4A-RAF':         ###Updated for C4###
            num_res = RAS4A_RAF_num_res
########## Above hard codes the number of residues/beads in the RAS-ONLY and RAS-RAF simulations #########################

        krases0_BB = u.select_atoms('resid '+str(start_of_g)+':'+str(start_of_g+num_res)+' and name BB')

        r2_26r40_56r69_166 = u.select_atoms('resid '+str(start_of_g)+':'+str(start_of_g+24)+' ' +\
                                                  str(start_of_g+38)+':'+str(start_of_g+54)+' ' +\
                                                  str(start_of_g+67)+':'+str(start_of_g+164)+\
                                                  ' and (name CA or name BB)')

        u_selection = \
            'resid '+str(start_of_g)+':'+str(start_of_g+24)+' '+str(start_of_g+38)+':'+str(start_of_g+54)+' ' +\
            str(start_of_g+67)+':'+str(start_of_g+164)+' and (name CA or name BB)'
        mobile0 = u.select_atoms(str(u_selection)).positions - u.select_atoms(str(u_selection)).center_of_mass()

        R, RMSD_junk = align.rotation_matrix(mobile0, ref0)

        lipids = u.select_atoms('resname POPX POPC PAPC POPE DIPE DPSM PAPS PAP6 CHOL')

        coords = ref0
        RotMat = []
        OS = []
        r152_165 = krases0_BB.select_atoms('resid '+str(start_of_g+150)+':'+str(start_of_g+163)+' and (name CA or name BB)')
        r65_74 = krases0_BB.select_atoms('resid '+str(start_of_g+63)+':'+str(start_of_g+72)+' and (name CA or name BB)')
        timeframes = []

        make_whole(krases0_BB)

        j, rmsd_junk = mda.analysis.align.rotation_matrix((r2_26r40_56r69_166.positions-r2_26r40_56r69_166.center_of_mass()), coords)
        RotMat.append(j)
        OS.append(r65_74.center_of_mass()-r152_165.center_of_mass())
        timeframes.append(u.trajectory.time)

        if protein_type == 'RAS-RAF' or protein_type == 'RAS4A-RAF':           ###Updated for C4###
            z_pos = []
############### NEED TO CONFIRM THE SELECTION OF THE RAF LOOP RESIDUES BELOW ####################
            if protein_type == 'RAS-RAF':
                zshifting=-2         ###Updated for C4###
            elif protein_type == 'RAS4A-RAF':
                zshifting=-1        ###Updated for C4### ##CHECK THIS!!!!!!!!!!!!!!!!!!!!
            raf_loops_selection = u.select_atoms('resid '+str(start_of_g+zshifting+291)+':'+str(start_of_g+zshifting+294)+' ' +\
                                                 str(start_of_g+zshifting+278)+':'+str(start_of_g+zshifting+281)+' ' +\
                                                 ' and (name CA or name BB)')
############### NEED TO CONFIRM THE SELECTION OF THE RAF LOOP RESIDUES ABOVE ####################
            diff = (lipids.center_of_mass()[2]-raf_loops_selection.center_of_mass(unwrap=True)[2])/10
            if diff < 0:
                diff = diff+(u.dimensions[2]/10)
            z_pos.append(diff)
            z_pos = np.array(z_pos)

        RotMatNP = np.array(RotMat)

        OS = np.array(OS)
        OA = RotMatNP[:, 2, :]/(((RotMatNP[:, 2, 0]**2)+(RotMatNP[:, 2, 1]**2)+(RotMatNP[:, 2, 2]**2))**0.5)[:, None]
        OWAS = np.arccos(RotMatNP[:, 2, 2])*180/math.pi
        OC_temp = np.concatenate((OA, OS), axis=1)
        t = ((OC_temp[:, 0]*OC_temp[:, 3])+(OC_temp[:, 1]*OC_temp[:, 4]) +
             (OC_temp[:, 2]*OC_temp[:, 5]))/((OC_temp[:, 0]**2)+(OC_temp[:, 1]**2)+(OC_temp[:, 2]**2))
        OC = OA*t[:, None]
        ORS_tp = np.concatenate((OC, OS), axis=1)
        ORS_norm = (((ORS_tp[:, 3]-ORS_tp[:, 0])**2)+((ORS_tp[:, 4]-ORS_tp[:, 1])**2)+((ORS_tp[:, 5]-ORS_tp[:, 2])**2))**0.5
        ORS = (OS - OC)/ORS_norm[:, None]
        OACRS = np.cross(OA, ORS)
        OZCA = OA * OA[:, 2][:, None]
        Z_unit = np.full([len(OZCA), 3], 1)
        Z_adjust = np.array([0, 0, 1])
        Z_unit = Z_unit*Z_adjust
        Z_OZCA = Z_unit-OZCA
        OZPACB = Z_OZCA/((Z_OZCA[:, 0]**2+Z_OZCA[:, 1]**2+Z_OZCA[:, 2]**2)**0.5)[:, None]
        OROTNOTSIGNED = np.zeros([len(ORS)])
        for i in range(len(ORS)):
            OROTNOTSIGNED[i] = np.arccos(np.dot(OZPACB[i, :], ORS[i, :]) /
                                         (np.sqrt(np.dot(OZPACB[i, :], OZPACB[i, :]))) *
                                         (np.sqrt(np.dot(ORS[i, :], ORS[i, :]))))*180/math.pi

        OZPACBCRS_cross = np.cross(OZPACB, ORS)
        OZPACBCRS = OZPACBCRS_cross/((OZPACBCRS_cross[:, 0]**2+OZPACBCRS_cross[:, 1]**2+OZPACBCRS_cross[:, 2]**2)**0.5)[:, None]
        OFORSIGN_temp = (OA - OZPACBCRS)**2
        OFORSIGN = OFORSIGN_temp[:, 0]+OFORSIGN_temp[:, 1]+OFORSIGN_temp[:, 2]
        OROT = OROTNOTSIGNED

        for i in range(len(OROT)):
            if OROT[i] < 0:
                OROT[i] = -(OROT[i])
        for i in range(len(OROT)):
            if OFORSIGN[i] < 0.25:
                OROT[i] = -(OROT[i])

###### Below introduces new shift to account for upper vs. lower leaflet #####
        for i in range(len(OWAS)):
            OWAS[i] = abs(-(OWAS[i])+180) # made this an absolute value so that the tilt remains positive

        for i in range(len(OROT)):
            if OROT[i] < 0:
                OROT[i] = OROT[i]+180
            elif OROT[i] > 0:
                OROT[i] = OROT[i]-180
###### Above introduces new shift to account for upper vs. lower leaflet #####



###### Below might have to be updated to take into account the periodic nature of the rotation ######
        if protein_type == 'RAS-ONLY' or protein_type == 'RAS4A-ONLY':
            c_macrostate = RAS_ONLY_macrostate
            if protein_type == 'RAS4A-ONLY':
                c_macrostate = RAS4A_ONLY_macrostate
            states = np.zeros(len(OROT))
            for j in range(len(OROT)):
                diff0 = []
                for i in range(len(c_macrostate)):
                    #diff0.append([((RAS_ONLY_macrostate[i,0]-OWAS[j])**2+(RAS_ONLY_macrostate[i,1]-OROT[j])**2)**0.5, RAS_ONLY_macrostate[i,6]])
                    diff0.append([((c_macrostate[i,1]-OWAS[j])**2+(c_macrostate[i,0]-OROT[j])**2)**0.5, c_macrostate[i,5]])
                diff0.sort()
                states[j] = diff0[0][1]
        elif protein_type == 'RAS-RAF' or protein_type == 'RAS4A-RAF':
            states = np.zeros(len(OROT))
            for j in range(len(OROT)):
### below: adding in the requirements for the 'high-z' state ###
                if (OROT[j] < -45 or OROT[j] > 140) and z_pos[j] > 4.8:
                    states[j] = 3
                else:
### above: adding in the requirements for the 'high-z' state ###
                    diff0 = []
                    for i in range(len(RAS_RAF_macrostate)):
                        #diff0.append([((RAS_RAF_macrostate[i,0]-OWAS[j])**2+(RAS_RAF_macrostate[i,1]-OROT[j])**2)**0.5, RAS_RAF_macrostate[i,6]])
                        diff0.append([((RAS_RAF_macrostate[i,1]-OWAS[j])**2+(RAS_RAF_macrostate[i,0]-OROT[j])**2)**0.5, RAS_RAF_macrostate[i,4]])
                    diff0.sort()
                    states[j] = diff0[0][1]
###### Above might have to be updated to take into account the periodic nature of the rotation ######



###### Assume we want to remove this? Where is the code that reads this information? i.e. will there be knock-on effects? ######
###### If feedback code needs index 5 (two_states) from the output, deleting this four_states will shift that to index 4 #######
#         four_states = np.zeros(len(OROT))
#         for j in range(len(OROT)):
#             diff0 = []
#             for i in range(len(macrostate4)):
#                 diff0.append([((macrostate4[i,0]-OWAS[j])**2+(macrostate4[i,1]-OROT[j])**2)**0.5, macrostate4[i,6]])
#             diff0.sort()
#             four_states[j] = diff0[0][1]+1


###### below: old output details.... ######################################
###### Updated - RAS-only to NOT HAVE the Z-distance ######################
###### Updated - Added in the protein 'tag', i.e. RAS-ONLY or RAS-RAF #####
#         OUTPUT = np.zeros([len(OROT), 6])

#         for i in range(len(OROT)):
#             OUTPUT[i] = timeframes[i], OWAS[i], OROT[i], z_pos[i], four_states[i], two_states[i]
###### above: old output details.... ######################################

###### below: NEW output details.... ######################################
        if protein_type == 'RAS-ONLY' or protein_type == 'RAS4A-ONLY':
            OUTPUT = np.zeros([len(OROT), 6]).astype(object)

            for i in range(len(OROT)):
                OUTPUT[i] = str(protein_type), timeframes[i], OWAS[i], OROT[i], 'n/a', int(states[i])
        elif protein_type == 'RAS-RAF' or protein_type == 'RAS4A-RAF':
            OUTPUT = np.zeros([len(OROT), 6]).astype(object)

            for i in range(len(OROT)):
                OUTPUT[i] = str(protein_type), timeframes[i], OWAS[i], OROT[i], z_pos[i], int(states[i])




        ALLOUT.append(OUTPUT)
    return np.asarray(ALLOUT)
    #np.savetxt(str(tpr)+"_tilt_rot_z_state.KRAS_"+str(k+1)+".txt", OUTPUT, fmt=['%i','%10.3f','%10.3f','%10.3f','%i','%i'], delimiter=' ')
