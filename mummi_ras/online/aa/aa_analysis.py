"""
Functions for online analysis of all-atom (AA) simulations.
"""

import os
import time
import argparse
import traceback
import shutil
import io
import tarfile
from time import sleep
import numpy as np

import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSF
# from MDAnalysis.analysis import distances
import MDAnalysis.transformations as trans
from mummi_core.utils.logger import init_logger

import mummi_core
from mummi_core.utils import Naming
import mummi_ras
from mummi_ras import get_io
from .aa_simulation import AAsimulation
from .mock_process import MockProcess, MockAAsimulation
from . import aa_get_tiltrot_z_state

from logging import getLogger
LOGGER = getLogger(__name__)


'''
This module analyzes AA data.

Default is for AMBER files (using .gro/.mdcrd as topo/traj for mdanalysis loading).
Other file types can be handled.

For an example run script, see: mummi/tests/analysis_aa.sh

'''


# ------------------------------------------------------------------------------
def setup_paths(args):
    """Extract localpath and outpath from args"""
    locpath = args.locpath
    outpath = args.outpath

    # default feedback path is the full path given by args.fbpath
    #fbpath = args.fbpath

    assert os.path.isdir(outpath)

    #@TODO this is not right - fix
    #if locpath == "":
    if True:
        locpath = outpath
    else:
        assert os.path.isdir(locpath)

    # from Harsh: switch fbpath to namespace if using databroker
    fbpath = ''

    return outpath, locpath, fbpath


def setup_filetypes(args):
    """Determine topology and trajectory file types for loading"""

    # Trajectory filename
    traj = args.trajname
    topo = args.toponame

    # Check for file extensions.
    # Assume AMBER (our gro/mdcrd combo for mdanalysis) if none given.
    trajtype = os.path.splitext(traj)[1]
    topotype = os.path.splitext(topo)[1]

    if trajtype == '':
        LOGGER.warning('No file extension given for trajname. Assuming .mdcrd')
        traj = '%s.mdcrd' % (traj)

    if topotype == '':
        LOGGER.warning('No file extension given for toponame. Assuming .gro')
        topo = '%s.gro' % (topo)

    return topo, traj, trajtype


# ------------------------------------------------------------------------------
def aa_load(topo, traj, trajtype, trajlist, aasimulation):
    """Load AA trajectory into Universe"""

    # first check if maxsimtime is reached
    limit_reached = \
        bool(aasimulation.get_current_simtime() >= aasimulation.get_maxsimtime())

    # if maxsimtime is reached, don't load any new files (e.g. md.10.mdcrd)
    # and don't checkstate
    if limit_reached:
        LOGGER.info("maxsimtime is reached; not loading any new trajectory files")
    # create a list from sequence of traj files
    # and check the state of mdcrd trajectory
    elif traj not in trajlist:
        trajlist.append(traj)
        aasimulation.mdcrd_checkstate()
    # check the state of mdcrd trajectory
    else:
        aasimulation.mdcrd_checkstate()

    # first try to load AMBER gro/mdcrd files
    try:
        LOGGER.info("topo = {}".format(topo))
        LOGGER.info("trajlist = {}".format(trajlist))
        # need to specify dt to prevent a get_dt error when 1 frame
        # hardcode dt=100 for now; will update as variable from aa_simulation
        # aa = mda.Universe(topo, traj, format='NCDF', dt=100)  # load single traj
        aa = mda.Universe(topo, trajlist, format='NCDF', dt=100)  # load list of traj
    except Exception as e:

        LOGGER.error(e)
        # aa = mda.Universe(topo, traj)  # load single traj
        #TODO: it is a hack for now to process unfinish frames
        # Chris: Xiaohua please check new logic above, if we need the [:-1] here anymore?
        aa = mda.Universe(topo, trajlist[:-1], format='NCDF', dt=100)
        #aa = mda.Universe(topo, trajlist)  # load list of traj

    return aa


def make_selections(aa):
    protein = aa.select_atoms('protein or resname ACE CYP CYF NMA')

    # individual selections to reorder same as CG (for feedback)
    gtp = aa.select_atoms('resname GTP')
    mg = aa.select_atoms('resname MG')
    zn = aa.select_atoms('resname ZN2')
    wat = aa.select_atoms('resname CGW')

    # gtp_ions_wat = gtp + mg + zn + wat

    # combined = protein + gtp_ions_wat
    combined = protein + gtp + mg + zn + wat

    # select everything except protein and water (SOL) and free ions (NA, CL)
    # needed for wrap PBC calc
    not_protein = aa.select_atoms('not protein and not (resname SOL NA CL)')

    # use a list of dict items to select based on n_residues:
    protein_systems = [{'186':'ras4a', '324':'ras4araf', '185':'ras', '323':'rasraf'}]  # ras=kras4b

    # count protein residues to determine which of 4 protein_systems
    protein_n_residues = protein.n_residues
    protein_type = [item[str(protein_n_residues)] for item in protein_systems][0]
    LOGGER.info('protein_type: {}'.format(protein_type))

    # make selections based on protein_type
    if protein_n_residues == 186:  # kras4a
        hvr = combined.residues[168:186]  # 18 residues
        crd = combined.residues[0:0]  # (empty residue group)

    elif protein_n_residues == 324:  # kras4a-rbdcrd
        hvr = combined.residues[168:186]  # 18 residues
        crd = combined.residues[268:322]  # 54 residues

    elif protein_n_residues == 185:  # kras4b
        hvr = combined.residues[168:185]  # 17 residues
        crd = combined.residues[0:0]  # (empty residue group)

    elif protein_n_residues == 323:  # kras4b-rbdcrd
        hvr = combined.residues[168:185]  # 17 residues
        crd = combined.residues[267:321]  # 54 residues

    else:
        LOGGER.error('Check sequence - number of residues does not match known seqeuences.')
        # raise Exception('Error: unknown sequence in KRas structure')

    # now make hvr + crd selection
    hvr_crd = hvr + crd

    '''
    kras4a:
    hvr.resnames
    ['LYS', 'LYS', 'ILE', 'SER', 'LYS', 'GLU', 'GLU', 'LYS', 'THR',
       'PRO', 'GLY', 'CYP', 'VAL', 'LYS', 'ILE', 'LYS', 'LYS', 'CYF']

    kras4b:
    hvr.resnames
    ['LYS', 'MET', 'SER', 'LYS', 'ASP', 'GLY', 'LYS', 'LYS', 'LYS',
       'LYS', 'LYS', 'LYS', 'SER', 'LYS', 'THR', 'LYS', 'CYF']

    crd.resnames
    ['VAL', 'PRO', 'LEU', 'THR', 'THR', 'HIS', 'ASN', 'PHE', 'ALA',
       'ARG', 'LYS', 'THR', 'PHE', 'LEU', 'LYS', 'LEU', 'ALA', 'PHE',
       'CYM', 'ASP', 'ILE', 'CYM', 'GLN', 'LYS', 'PHE', 'LEU', 'LEU',
       'ASN', 'GLY', 'PHE', 'ARG', 'CYM', 'GLN', 'THR', 'CYM', 'GLY',
       'TYR', 'LYS', 'PHE', 'HIS', 'GLU', 'HIS', 'CYM', 'SER', 'THR',
       'LYS', 'VAL', 'PRO', 'THR', 'MET', 'CYM', 'VAL', 'ASP', 'TRP']
    '''

    # Guess bonds needed for PBC calc below
    # takes a while to run; quicker to guess bonds on selections than all Universe
    # can pass for some traj that don't need and complain
    try:
        protein.guess_bonds()
        # not_protein.guess_bonds()
    except:
        pass

    # TODO: c4 update, return protein_type
    # return protein, combined, not_protein, hvr, crd, hvr_crd
    return protein_type, protein, combined, not_protein, hvr, crd, hvr_crd


# ------------------------------------------------------------------------------
# Save a checkpoint file that can be used to restore/restart
def checkpoint(frame, step, simname, toponame, trajname, trajtype, trajlist,
               outpath, fbpath):

    # TODO: this function assumes checkpointing using npz and assumes simple io
    LOGGER.info('Checkpointing: frame = {}, step = {}'.format(frame, step))

    ddict = {'fcount': frame,
             'step': step,
             'simname': simname,
             'toponame': toponame,
             'trajname': trajname,
             'trajtype': trajtype,
             'trajlist': trajlist,
             'outpath': outpath,
             # 'locpath': locpath,
             #'fbpath': fbpath,
             }

    # TODO: CS, may want better file naming
    filename = 'aa_checkpoint.npz'

    # check for existing checkpoint file and make a backup before saving new
    file_path = os.path.join(outpath, filename)
    if os.path.isfile(file_path):
        backup = os.path.join(outpath, filename[:-4]+'_backup.npz')
        shutil.copy(file_path, backup)

    get_io('simple').save_npz(outpath, filename, ddict)
    LOGGER.info('Checkpointed to ({})'.format(filename))


# Restore from a previous checkpoint file
def restore(outpath, args):

    filename = os.path.join(outpath, 'aa_checkpoint.npz')
    LOGGER.info('Restoring from ({})'.format(filename))

    success = True
    try:
        args_restore = np.load(filename)

        args.simname = args_restore['simname'].item(0)
        args.fcount = args_restore['fcount'].item(0)
        args.step = args_restore['step'].item(0)

        args.outpath = args_restore['outpath'].item(0)
        # args.fbpath = args_restore['fbpath'].item(0)
        # args.locpath = args_restore['locpath'].item(0)

        # now advance the start frame by step interval
        args.fcount += args.step

        toponame = args_restore['toponame'].item(0)
        trajname = args_restore['trajname'].item(0)
        trajtype = args_restore['trajtype'].item(0)
        trajlist = args_restore['trajlist'].tolist()

    except Exception as e:
        LOGGER.error(e)
        success = False
        toponame = ''
        trajname = ''
        trajtype = ''
        trajlist = []

    return success, args, toponame, trajname, trajtype, trajlist


# ------------------------------------------------------------------------------
# Individual analysis functions

# Correct for periodic boundary conditions (PBC)
# Transform can only be applied once on a new Universe
def transform_PBC(aa, protein, not_protein):
    farnesyl = protein.select_atoms('resname CYF')
    all_sel = protein + not_protein
    ag = all_sel.atoms

    transforms = (trans.unwrap(ag),
                  # trans.center_in_box(protein, center='mass'),
                  trans.center_in_box(farnesyl, center='mass'),
                  trans.translate([0,0,20]),
                  trans.wrap(ag, compound='fragments'))

    aa.trajectory.add_transformations(*transforms)

    return aa


# Save protein xyz coordinates, in array and pdb file
def calc_coord(simname, protein, frame, outpath):

    traj_coord = []

    traj_coord.append(protein.positions)

    # Turning off pdb write here; using taridx save_pdb_data_aa()
    # filename = os.path.join(outpath, Naming.aaframe(simname, frame)+'.pdb')

    # LOGGER.debug('writing ({})'.format(filename))
    # protein.write(filename)

    # convert list to np.array
    traj_coord = np.asarray(traj_coord)

    return traj_coord


# RMSD calculation
def calc_rmsd(ag, ref, select):

    # aa.trajectory[0] # first frame
    # traj_ref = ag.select_atoms(select).positions

    # aa.trajectory[frame] # current frame
    traj_frame = np.asarray(ag.select_atoms(select).positions, dtype=np.float64)

    traj_rmsd = rms.rmsd(traj_frame, ref, superposition=True)

    return traj_rmsd


def calc_protalign(aa):

    # TODO: pass protein selection in from first definition
    # protein = aa.select_atoms('protein or resname ACE CYF NMA')
    protein = aa.select_atoms('protein or resname ACE CYP CYF NMA')

    # select = 'protein and name CA'
    # select = '(protein or resname ACE CYF NMA) and name CA'
    select = '(protein or resname ACE CYP CYF NMA) and name CA'
    # select = 'protein and backbone'

    # superimpose the frames based on selected region of protein,
    # which will be used to calculate rmsf
    prealigner = align.AlignTraj(aa, aa,
                                 select=select,
                                 in_memory=True).run()

    ref_coordinates = aa.trajectory.timeseries(asel=protein).mean(axis=1)

    ref = mda.Merge(protein).load_new(ref_coordinates[:, None, :], order="afc")

    aligner = align.AlignTraj(aa, ref,
                              select=select,
                              in_memory=True).run()
    return aa, protein


# RMSF calculation
def calc_rmsf(aa, frame_start, frame_stop, step):
    """ note that .run() can take the following:
        run(start=None, stop=None, step=None, verbose=None)
    """

    # pass through align
    aa, protein = calc_protalign(aa)

    # TODO: CS, hardcode to start/test; can pull this back out later as an option
    protein_sel = 'name CA'
    # protein_sel = 'backbone'

    # get RMSF per residue of selected protein region
    # based on the average fluctuations
    protein_region = protein.select_atoms(protein_sel)

    # Calculate RMSF over cumulative frames (and step = 1)
    rmsfer = RMSF(protein_region).run()
    # traj_rmsf = np.array([[[0, frame_stop, 1], [rmsfer.rmsf]]])
    traj_rmsf = np.array([[rmsfer.rmsf]])

    # Calculate RMSF over specified frames only
    # rmsfer = RMSF(protein_region).run(frame_start, frame_stop, step)
    # traj_rmsf = np.array([[[frame_start, frame_stop+1, step], [rmsfer.rmsf]]])

    LOGGER.debug('Computed rmsf: {}'.format(rmsfer.rmsf))

    return traj_rmsf


# Calculates HVR end-to-end distance
def calc_hvr_dist(hvr):

    hvr_dist = []

    # TODO: CS, make selections more general
    # select N-term and C-term of HVR region
    # K-Ras-4B HVR: res 169-185 (KMSKDGKKKKKKSKTKC)
    # ref: A. Stephen, Cancer Cell (2014)

    hvr_nterm = hvr.atoms.select_atoms('name N')[0]
    hvr_cterm = hvr.atoms.select_atoms('name C')[-1]

    r = hvr_cterm.position - hvr_nterm.position  # end-to-end vector from atom positions
    d = np.linalg.norm(r)  # end-to-end distance

    hvr_dist.append(d)

    # convert to np.array
    hvr_dist = np.asarray(hvr_dist)

    return hvr_dist


# tilt/rotation/depth analysis
# adapted from cganalysis function: processFrame()
def calc_tiltrot_z(frame, aa, kras_indices):
    """Run all analysis returns a dictionary of all results."""

    return_dic = {"frameNum": frame}

    # bbox = [0, cg.syst.dimensions[0], 0, cg.syst.dimensions[1]]

    # ####################################################
    # A2:    KRAS states
    # For each KRAS return (KRAS number, t_1 angle, t_2 angle, z dist from
    # bilayer, full state[1-4], main state[1-2])
    # output file is saved as numpy matix (number of KRAS x 6)

    # TODO: CS, my edit testing
    # statesList = \
    #     mummi_ras.online.cg.get_tiltrot_z_state_multiple.getKRASstates(cg.syst, options.cgSel.kras_indices)
    statesList = aa_get_tiltrot_z_state.getKRASstates(aa, kras_indices)


    # Update state number to state name - maybe shoudl be part of the files ???
    #print("HII-before " + str(statesList))
    ''' Refactor - use code form Tomas
    statemap = {  "RAS-ONLY" : { 1:"a"  , 2:"b"  , 3:"c"  } \
                  "RAS-RAF"  : { 1:"ma" , 2:"mb" , 3:"za" } }
    for pState in statesList:
        for innerState in pState:
            try:
                innerState[5] = statemap[innerState[0]][innerState[5]]
            except:
                LOGGER.error("-get, States, protein type {} state {} not supported".format(innerSate[0], innerSate[5]))
                raise Exception("Error in protein type definition {}, state {} not supported".format(innerSate[0], innerSate[5]))
    return_dic["states"] = statesList
    '''
    for pState in statesList:
        for innerSate in pState:
            if innerSate[0] == "RAS-ONLY":
                if innerSate[5] == 1:
                    innerSate[5] = "a"
                elif innerSate[5] == 2:
                    innerSate[5] = "b"
                elif innerSate[5] == 3:
                    innerSate[5] = "c"
                else:
                    LOGGER.error("-get, States, protein type {} state {} not supported".format(innerSate[0], innerSate[5]))
                    raise Exception("Error in protein type definition {}, state {} not supported".format(innerSate[0], innerSate[5]))
            elif innerSate[0] == "RAS-RAF":
                if innerSate[5] == 1:
                    innerSate[5] = "ma"
                elif innerSate[5] == 2:
                    innerSate[5] = "mb"
                elif innerSate[5] == 3:
                    innerSate[5] = "za"
                else:
                    LOGGER.error("-get, States, protein type {} state {} not supported".format(innerSate[0], innerSate[5]))
                    raise Exception("Error in protein type definition {}, state {} not supported".format(innerSate[0], innerSate[5]))
            else:
                LOGGER.error("-get, States, protein type {} not supported".format(innerSate[0]))
                raise Exception("Error in protein type definition, type {} not supported".format(innerSate[0]))
    return_dic["states"] = statesList
    #print("HII-after " + str(statesList))

    if len(statesList) > 0:
        LOGGER.info(
            "-get, protein states for shape {}, first is {}".format(statesList[0].shape, statesList[0][0]))
    else:
        # This happens if there are no proteins so a lipid only patch
        LOGGER.info("-get, we have a lipid only patch")

    # TODO: CS, to remove; for original hardcode testing for fb save criteria
    # current_crd_depth = 1
    # return current_crd_depth

    # Return all data as a dictonary
    return return_dic


# OPTIONAL ANALYSIS:

# Calculates protein-lipid contacts based on a cutoff distance
def calc_pl_contacts_simple(aa):
    """PLACEHOLDER

        1. finds all protein atoms within x (cutoff) of any lipid atom
           make selection as those protein residues

        2. selects all lipid atoms within x of those protein atoms
           make selection as those lipids
    """

    # arr_pl_contacts = np.empty((0,2))
    arr_pl_contacts = []

    # TODO: CS, turn into param to set; right now hardcoded
    contact_cutoff = 12.0  # (Angstroms) cutoff distance to use as 'interaction'

    # loops over traj subset
    # for ts in aa.trajectory[frame_start:frame_stop+1:step]:

    # select protein within cutoff of lipid
    protein = aa.select_atoms('protein and around ' + str(contact_cutoff) + ' (resname POPC POPC)')

    # select lipids within cutoff of protein
    lipid = aa.select_atoms('(resname POPC POPC) and (around ' + str(contact_cutoff) + ' protein)')

    # just capture unique resnums as arrays
    prot_residues = protein.residues.ix_array
    lip_residues = lipid.residues.ix_array

    # item = [[prot_residues, lip_residues]]
    item = (prot_residues, lip_residues)

    # arr_pl_contacts = np.append(arr_pl_contacts, item, axis=0)
    arr_pl_contacts.append(item)

    # convert to np.array
    arr_pl_contacts = np.asarray(arr_pl_contacts)

    return arr_pl_contacts


# Calculates protein-lipid COM distance
def calc_pl_comdist(aa):

    protein = aa.select_atoms('protein')
    lipid = aa.select_atoms('resname POPS POPC')

    # arr_pl_comdist = np.empty((0,1))
    arr_pl_comdist = []

    # for ts in aa.trajectory[frame_start:frame_stop+1:step]:

    protein_com = protein.center_of_mass()
    lipid_com = lipid.center_of_mass()

    prot_lip_dist = np.linalg.norm(protein_com - lipid_com)

    # arr_pl_comdist = np.append(arr_pl_comdist, [[prot_lip_dist]], axis=0).astype(object)
    arr_pl_comdist.append(prot_lip_dist)

    # convert to np.array
    arr_pl_comdist = np.asarray(arr_pl_comdist)

    return arr_pl_comdist


# ------------------------------------------------------------------------------
# save all analysis
def npz_save(protein_type, simname, frame,
             outpath, fbpath,
             ioi, iofb,
             fb_start, fb_crd_depth, current_crd_depth,
             traj_coord,
             rmsd_hvr,
             rmsd_crd,
             rmsd_hvr_crd,
             traj_rmsf,
             traj_hvr_dist,
             traj_tiltrot):
             # arr_pl_contacts,
             # arr_pl_comdist):

    # TODO: CS, ** other way to attach name is update naming.py, e.g. 'aaframe': {rastype}_{simname}_aaf{aaf:08d}
    # aa_filename = Naming.aaframe(simname, frame)
    basename = protein_type + '_' + simname
    aa_filename = Naming.aaframe(basename, frame)

    # create a dict of all data needed to be stored
    ddict = {'frameNum': frame,
              'pdb_coord': traj_coord,
              # rmsd=rmsd,
              'rmsd_hvr': rmsd_hvr,
              'rmsd_crd': rmsd_crd,
              'rmsd_hvr_crd': rmsd_hvr_crd,
              # rmsf=rmsfer.rmsf,
              'rmsf': traj_rmsf,
              'hvr_dist': traj_hvr_dist,
              'tiltrot_z': traj_tiltrot
              }
              # 'pl_contacts': arr_pl_contacts,
              # 'pl_comdist': arr_pl_comdist
              # }

    # save actual arrays
    # turning off outpath save; instead using tar save to outpath below
    # ioi.save_npz(outpath, aa_filename, ddict, interfaces=['simple'])

    # check criteria for feedback save
    if (frame >= fb_start) and (current_crd_depth < fb_crd_depth):
        # fdinatal -- removal of dbr due to hang
        # iofb.save_npz(fbpath, aa_filename, ddict, interfaces=['simple', 'dbr'])
        iofb.save_npz(fbpath, aa_filename, ddict, interfaces=['simple'])
        LOGGER.info('feedback criteria met; saving frame {} for feedback'.format(frame))
    else:
        LOGGER.info('feedback criteria was not met for frame {}'.format(frame))
        pass

    # also create a list of dict to pass into tar file save
    llist = [ddict]

    try:
        # TODO: CS, check to update naming for x4 types here (also pdb save)
        get_io('mummi').save_analysis_data(llist, outpath)

        # count files in tar file (to determine success)
        with tarfile.open(os.path.join(outpath, 'analysis.tar')) as archive:
            count_npz_tar = sum(1 for member in archive if member.isreg())

        ddict_analysis_keys = len(ddict.keys()) - 1  # subtract 1 for frameNum
        if ddict_analysis_keys != 0:
            count_frames_tar = count_npz_tar // ddict_analysis_keys

    except Exception as e:
        LOGGER.error(e)

    return count_frames_tar


# ------------------------------------------------------------------------------
def setup_argparser():
    """Set up the program's argument parser."""

    parser = argparse.ArgumentParser(description='Analyze AA simulations individually.')

    parser.add_argument('--simname',
                        help='Name of the simulation')

    parser.add_argument('--trajname',
                        help='Name of the trajectory file')

    parser.add_argument('--toponame', default=None,
                        help='Name of the topology file (default: same as trajname)')

    parser.add_argument("--step", type=int, default=1,
                       help="Frame step to use for analysis (default: 1)")

    parser.add_argument("--fcount", type=int, default=0,
                       help="Frame number for inital frame (default: 0)")

    parser.add_argument("--outpath", required=True, type=str,
                       help="Path for input and output")

    parser.add_argument("--locpath", type=str, default="",
                       help="If set, use local path, and periodically copy files to outpath")

    parser.add_argument('--fstype', type=str, default="mummi",
                        choices=mummi_ras.get_interfaces(),
                        help='Choose how to write files/data (default: mummi)')

    parser.add_argument('--fbio', type=str, default="mummi",
                        choices=mummi_ras.get_interfaces(),
                        help='Choose how to write files/data for feedback (default: mummi)')

    parser.add_argument("-d", "--loglevel", type=int, default=2,
                        help="Level of logging messages to be output:\n"
                             "5 - Critical\n"
                             "4 - Error\n"
                             "3 - Warning\n"
                             "2 - Info (Default)\n"
                             "1 - Debug")

    parser.add_argument("-l", "--logpath", type=str, default="./",
                        help="Path for log information.")

    parser.add_argument("-c", "--logstdout", action="store_true",
                        help="Log to stdout in addition to a file.")

    parser.add_argument("--noninteractive", action="store_true",
                       help="Perform analysis on static MD files - no simulation running")

    parser.add_argument("--test", action="store_true",
                        help="Analyze a single structure file and exit.")

    parser.add_argument("--logfile", type=str, default="aa_analysis",
                        help="Name of aa analysis logger")

    parser.add_argument("--backuprate", type=int, default=5,
                        help="Name of aa analysis logger")

    parser.add_argument("--maxsimtime", type=int, default=12500000,
                        help="Stop simulation and analysis after this time step is reached")
    return parser


# ------------------------------------------------------------------------------
def run_aa_analysis(simname, toponame, trajname, trajtype, trajlist,
                    frame_start, frame_step,
                    outpath, fbpath,
                    ioi, iofb,
                    process, aasimulation):
    """Perform online analysis of an AA simulation

        Current usage and defaults:
        step = 1        # frame step to take for analysis
        time.sleep(10)  # sleep for 10 sec before checking for new frames
    """


    # load trajectory
    aa = aa_load(toponame, trajname, trajtype, trajlist, aasimulation)
    LOGGER.info(aa)

    # get number of frames to process
    frame_stop = len(aa.trajectory)
    frames_total = (frame_stop - frame_start)//frame_step

    # make all needed selections from universe
    # also apply protein.guess_bonds() needed for PBC
    # protein, combined, not_protein, hvr, crd, hvr_crd = make_selections(aa)

    # --------- Online loop
    # keep looking for new frames to process and stop when no new frames

    # Initialize Loop Booleans
    guarantee_count=0
    sim_running = True if process and process.poll() is None else False
    frames_left = True if frame_start < frame_stop else False

    frame0_initialized = False
    hvr_frame0 = None
    crd_frame0 = None
    hvr_crd_frame0 = None

    LOGGER.info(
        f'BEFORE LOOP: found {frame_stop} frames in trajectory. '
        f'will process total {frames_total}')
    LOGGER.debug("BEFORE LOOP: MD Running? %s", bool(process and process.poll() is None))

    while sim_running or frames_left:

        LOGGER.info(
            f'LOOP START: found {frame_stop} frames in trajectory. '
            f'will process total {frames_total}')
        LOGGER.debug("IN LOOP: MD Running? %s", bool(process and process.poll() is None))

        # make all needed selections from universe; keep in loop to update with each aa_load()
        # also apply protein.guess_bonds() needed for PBC
        # TODO: c4 update, new var
        # protein, combined, not_protein, hvr, crd, hvr_crd = make_selections(aa)
        protein_type, protein, combined, not_protein, hvr, crd, hvr_crd = make_selections(aa)
        # Apply PBC transform
        # (within loop to update for new frames found)
        aa = transform_PBC(aa, protein, not_protein)

        # (these selections work for x4 types in C4)
        if not frame0_initialized:
            hvr_frame0 = np.asarray(hvr.atoms.select_atoms('backbone').positions.copy(), dtype=np.float64)
            crd_frame0 = np.asarray(crd.atoms.select_atoms('backbone').positions.copy(), dtype=np.float64)
            hvr_crd_frame0 = np.asarray(hvr_crd.atoms.select_atoms('backbone').positions.copy(), dtype=np.float64)
            frame0_initialized = True

        # TODO: CS, old selections method; keep as placeholder while testing because:
        # (fragments need PBC calc above to be run first)
        # hvr = protein.fragments[0].residues[-17:]
        # crd = protein.fragments[1].residues[-56:-2]

        # capture frame 0 coords *after* PBC and *before* analysis (ts) loop
        # Important: these are needed in calc_rmsd()

        # Runs each analysis per frame

        # 1) RMSF of protein
        # note: RMSF needs to be calculated across multi-frames so setting above loop
        # traj_rmsf = calc_rmsf(aa, frame_start, frame_stop, frame_step)
        # LOGGER.info('Calculated traj_rmsf')

        # TODO: CS, revisit - applying a placeholder RMSF calc here (with NaNs array)
        # 1/11/2021 - current issue is calc_protalign() takes too long for online
        traj_rmsf = np.empty((1,1,319))
        traj_rmsf[:] = np.NaN
        LOGGER.info('Passed traj_rmsf')

        # Perform analysis for all new frames found
        for ts in aa.trajectory[frame_start:frame_stop+1:frame_step]:

            frame = ts.frame

            LOGGER.info(f'starting frame {frame}')

            # calculate protein Rg and check if reasonable
            # i.e. check PBC transform applied is reasonable
            protein_rg = protein.radius_of_gyration()
            if protein_rg > aa.dimensions[2]*0.4:
                LOGGER.warning('frame {} has protein Rg = {:.2f} A; check if PBC is correct'.format(frame, protein_rg))
            else:
                pass

            # TODO: CS, check to update naming for x4 types
            # save frame pdb into tar file
            try:
                LOGGER.info(f'saving frame {frame} to pdb.tar')
                get_io('mummi').save_pdb_data_aa(frame, combined, outpath)
                LOGGER.info('successfully saved frame {} to pdb.tar'.format(frame))
            except Exception as e:
                LOGGER.error(e)

            # 2) Trajectory coordinates per frame
            # traj_coord = calc_coord(simname, protein, frame, path_analysis)
            traj_coord = calc_coord(simname, combined, frame, outpath)
            LOGGER.info('Calculated traj_coord')

            # 3a) RMSD of all protein
            # rmsd = calc_rmsd(aa, protein.atoms, 'backbone', frame)

            # 3b) RMSD of HVR
            try:
                rmsd_hvr = calc_rmsd(hvr.atoms, hvr_frame0, 'backbone')
                LOGGER.info('Calculated rmsd_hvr')
            except Exception as e:
                rmsd_hvr = np.NaN
                LOGGER.error(e)

            # wrap these 2 crd analyses in IF statement for x4 cases
            if (protein_type == 'ras4araf') or (protein_type == 'rasraf'):
                # 3c) RMSD of CRD
                try:
                    rmsd_crd = calc_rmsd(crd.atoms, crd_frame0, 'backbone')
                    LOGGER.info('Calculated rmsd_crd')
                except Exception as e:
                    rmsd_crd = np.NaN
                    LOGGER.error(e)

                # 3d) RMSD of HVR + CRD
                try:
                    rmsd_hvr_crd = calc_rmsd(hvr_crd.atoms, hvr_crd_frame0, 'backbone')
                    LOGGER.info('Calculated rmsd_hvr_crd')
                except Exception as e:
                    rmsd_hvr_crd = np.NaN
                    LOGGER.error(e)
            else:
                # no raf present in structure
                rmsd_crd = np.NaN
                rmsd_hvr_crd = np.NaN
                LOGGER.info('protein_type has no RAF; skipping CRD calculations')

            # 4) End-to-end distance of HVR region
            try:
                traj_hvr_dist = calc_hvr_dist(hvr)
                LOGGER.info('Calculated traj_hvr_dist')
            except Exception as e:
                traj_hvr_dist = np.asarray([np.NaN])
                LOGGER.error(e)

            # 5) Tilt/rotation/depth
            # pass protein_type to tiltrot_z and parse there into x4 options
            ''' note: 'RAS-ONLY' does not return a z-depth '''
            # kras_indices = [(2,'RAS-RAF')] # TODO: TSC, changed residue 1 to residue 2 due to ACE cap
            kras_indices = [(2, protein_type)]
            traj_tiltrot = calc_tiltrot_z(frame, aa, kras_indices)
            LOGGER.info('Calculated traj_tiltrot')

            # TODO: CS, hardcoding 2 variables for feedback save criteria
            #           also, add these to ckpt/restore
            # set criteria for saving to feedback
            fb_start = 0  # frame to reach before saving for feedback
            fb_crd_depth = 3  # crd threshold for saving for feedback

            # TODO: c4 update, acct for 'RAS-ONLY' no fb_crd_depth
            if (protein_type == 'ras4araf') or (protein_type == 'rasraf'):
            # if protein_n_residues >= 323:  # for 'kras4b-rbdcrd' or 'kras4a-rbdcrd'
                current_crd_depth = traj_tiltrot['states'][0][0][4]  # access z depth for frame
                # current_crd_depth = traj_tiltrot
            else:
                current_crd_depth = fb_crd_depth//2  # set to < fb_crd_depth to always save to fb

            # based on Rg check above, don't send frame to feedback if PBC may be incorrect
            if protein_rg > aa.dimensions[2]*0.4:
                current_crd_depth = fb_crd_depth*2  # set to > fb_crd_depth to not save to fb
                LOGGER.warning('frame {} will not be sent to feedback; check if PBC is correct'.format(frame))
            else:
                pass

            # OPTIONAL ANALYSIS:
            # turning these off for now

            # 6) Protein-lipid contacts
            # arr_pl_contacts = calc_pl_contacts_simple(aa)
            # LOGGER.info('Calculated arr_pl_contacts')

            # 7) Protein-lipid center-of-mass (COM) distance
            # arr_pl_comdist = calc_pl_comdist(aa)
            # LOGGER.info('Calculated arr_pl_comdist')

            # TODO: CS, local test; including protein_type to use in aa_filename
            # save all analysis as per-frame npz
            count_frames_tar = \
                npz_save(
                    protein_type,
                    simname, frame, outpath, fbpath, ioi, iofb, fb_start, fb_crd_depth, current_crd_depth,
                    traj_coord, rmsd_hvr, rmsd_crd, rmsd_hvr_crd, traj_rmsf, traj_hvr_dist, traj_tiltrot) # arr_pl_contacts, # arr_pl_comdist)

            # save a checkpoint upon completing analysis per frame
            checkpoint(frame, frame_step,
                        simname, toponame, trajname, trajtype, trajlist,
                        outpath, fbpath)
            LOGGER.debug("Frames in analysis.tar = %s", count_frames_tar)

        # wait before checking for new frames
        # current Amber frame rate is ~5 frames/hr
        LOGGER.info("AA analysis is waiting to check for new frames.")
        time.sleep(120)

        # TODO: CS, reset traj before re-loading below
        aa.trajectory.close()

        # update traj and look for new frames
        aa = aa_load(toponame, trajname, trajtype, trajlist, aasimulation)
        frame_start = frame_stop
        frame_stop = len(aa.trajectory)
        frames_total = (frame_stop - frame_start)//frame_step

        # UPDATE BOOLEANS FOR CONDITIONAL
        if process and process.poll() is None and guarantee_count < 1:
            sim_running = True
        else:
            guarantee_count = guarantee_count + 1
            sim_running = False

        frames_left = True if frame_start < frame_stop else False
        LOGGER.debug(
            "END LOOP: Simulation running is %s...\n"
            "END LOOP: MD Running? %s\n"
            "END LOOP: guarantee_count = %s\n",
            "END LOOP: Frames Left [start=%s, end=%s]? %s",
            sim_running, bool(process and process.poll() is None), guarantee_count,
            frame_start, frame_stop, frames_left
        )


    # 11/18:    Chris, please take a look at this one!
    # TODO: we should check if things finished successfully
    # did we analyze *all* frames or not

    # 11/30:    My criteria check is implemented below

    flag_success, flag_failure = Naming.status_flags('aa')

    # Method 1: can compare npz files in analysis.tar to total sim frames
    # get total frames generated from AA simulation
    # Note: ntwx=25000 is hardcoded in aa_siminputs
    # total_sim_frames = args.maxsimtime / 25000

    # if count_frames_tar == total_sim_frames:

    # Method 2: compare npz files in analysis.tar to total frames processed
    # *important* also ensure simulation reached completion by checking simtime
    limit_reached = \
        bool(aasimulation.get_current_simtime() >= aasimulation.get_maxsimtime())
    LOGGER.debug(
        "AFTER LOOP: Sim time reached? %s\n"
        "AFTER LOOP: Frame pointers -- start = %s, stop = %s",
    )
    if  limit_reached and frame_start >= frame_stop:
        ioi.send_signal(outpath, flag_success)
    else:
        ioi.send_signal(outpath, flag_failure)

    # Backup the file
    if aasimulation.get_islocal():
        aasimulation.backup()


# ------------------------------------------------------------------------------
def aa_analysis(simname, toponame, trajname, trajtype, trajlist,
                frame_start, frame_step,
                outpath, locpath, fbpath, ioi, iofb,
                process, aasimulation):
    #frame_start = int(args.fcount)  # default: 0
    #frame_step = int(args.step)  # default: 1

    LOGGER.info("trajname = {}".format(trajname))
    LOGGER.info("toponame = {}".format(toponame))
    LOGGER.info("trajtype = {}".format(trajtype))
    LOGGER.info("outpath = {}".format(outpath))
    LOGGER.info("locpath = {}".format(locpath))
    #LOGGER.info("fbpath = {}".format(fbpath))
    LOGGER.info("frame: start = {}, step = {}".format(frame_start, frame_step))

    # Run online analysis
    try:
        run_aa_analysis(simname, toponame, trajname, trajtype, trajlist,
                        frame_start, frame_step,
                        outpath, fbpath,
                        ioi, iofb,
                        process, aasimulation)

    except Exception as e:
        LOGGER.error(e)
        traceback.print_exc()
        raise e


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def main():
    aasimulation = None
    try:
        # Set up the necessary base data structures to begin study set up.
        parser = setup_argparser()
        args = parser.parse_args()

        # Do something with args, for testing
        print(vars(args))

        # Sets up MuMMI and the logger
        mummi_core.init()
        mummi_core.init_logger(argparser=args)

        # parse the arguments to get what is needed
        outpath, locpath, fbpath = setup_paths(args)

        # Initilize output DBR/files etc
        ioi = get_io(args.fstype)  # io interface # mummi
        iofb = get_io(args.fbio)   # io feedback  # mummi

        LOGGER.info('Using IO_Interface ({}) for data output'.format(ioi))
        LOGGER.info('Using IO_Interface ({}) for feedback output'.format(iofb))

        # HB added this to override argument!
        fbpath = Naming.dir_root('feedback-aa')

        # --------------------------------------------------------------------------
        # Option for entering non-interactive mode
        if args.noninteractive:

            # Setup a mock process
            process = MockProcess()

            # First try to restore from a saved checkpoint
            success, args, toponame, trajname, trajtype, trajlist = restore(outpath, args)

            toponame, trajname, trajtype = setup_filetypes(args)

            # Setup a mock simulation class
            aasimulation = MockAAsimulation()
            aasimulation.get_curIdx(trajname)

        else:  # all the orig stuff

            # --------------------------------------------------------------------------
            ## setup the simulation class
            aasimulation = AAsimulation(outpath, locpath)
            aasimulation.set_backuprate(args.backuprate)
            aasimulation.set_maxsimtime(args.maxsimtime)
            aasimulation.mdrun()

            # process = test_run(aasimulation.get_process)
            process = aasimulation.get_process()

            # --------------------------------------------------------------------------
            # First try to restore from a saved checkpoint
            success, args, toponame, trajname, trajtype, trajlist = restore(outpath, args)

            toponame = aasimulation.get_toponame()
            trajname = aasimulation.get_trajname()
            trajtype = aasimulation.get_trajtype()

        # --------------------------------------------------------------------------
        LOGGER.debug("Before Process Process Running? %s", process and process.poll() is None)

        aa_analysis(args.simname, toponame, trajname, trajtype, trajlist,
                    int(args.fcount), int(args.step),
                    outpath, locpath, fbpath, ioi, iofb,
                    process, aasimulation)

        # --------------------------------------------------------------------------

    except Exception as e:
        LOGGER.error(e)
        _, fail_flag = Naming.status_flags('aa')
        ioi.send_signal(outpath, fail_flag)
        raise e

    finally:
        LOGGER.info("When going down in a finally stop all simulations")
        if aasimulation is None:
            aasimulation.stop()

    # --------------------------------------------------------------------------
    '''
    LOGGER.debug("Before Process Process Running? %s", process.poll())

    p_backup = Process(target=aasimulation.backup)
    p_analysis = Process(target=aa_analysis, args=(args.simname, toponame, trajname,
                                                    int(args.fcount), int(args.step),
                                                    trajtype, outpath, locpath, fbpath, process, aasimulation))

    p_backup.start()
    p_analysis.start()

    LOGGER.debug("After start Process Running? %s", process.poll())
    p_analysis.join()
    p_backup.join()

    '''


# ------------------------------------------------------------------------------
if __name__ == '__main__':
    main()

# ------------------------------------------------------------------------------
