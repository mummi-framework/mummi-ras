# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
import os
import glob
import re
import yaml
import shutil
import time, datetime
import secrets
from collections import deque
from collections import Counter
import multiprocessing
import numpy as np
from zipfile import BadZipFile

import mummi_core
from mummi_core.utils import Naming
from mummi_core.utils.utilities import add_timestamp, remove_parallel, mkdir, sig_ign_and_rename_proc
from logging import getLogger
LOGGER = getLogger(__name__)

from mummi_core.workflow.feedback_manager import FeedbackManager, FeedbackManagerType


# ------------------------------------------------------------------------------
class FeedbackManager_AA2CG(FeedbackManager):

    # TODO: remove all arguments on 2nd and 3rd line
    def __init__(self, type, workspace, outpath, fbpath, iointerface, fbinterface, hvr_th, crd_th, nframes_fb_increment):

        assert type == FeedbackManagerType.Manager
        super().__init__(type, 'AA2CGFeedback_Manager')

        # this code has been written with this assumption
        assert workspace == outpath
        self.naggregates = 0            # number of aggregations
        self.nreports    = 0            # number of reported aggregations
        self.naggs_at_last_report = 0   # number of aggregations in last report

        self.workspace = workspace

        # outpath for the feedback (where final files will be saved for macro)
        self.iointerface = iointerface
        self.outpath = outpath

        # fbpath where communication between worker and manager will happen
        self.fbinterface = fbinterface
        self.fbpath = fbpath

        # give me a full path as the fb path
        # and i will figure out if i need to trim down to the dbr namespace
        assert len(os.path.dirname(self.fbpath)) > 0
        assert self.fbpath == self.outpath

        ## TODO: dont copy if they already exist!
        res_path = os.path.join(Naming.MUMMI_RESOURCES, 'kras-aa/feedback-aa2cg')

        # copy the resources
        files_to_copy = [
            'martinize_2.6_3_modified.py',
            'RAS_RAF_TERNARY_scfix-CYFpos.itp',
            'RAS_scfix-CYFpos.itp',
            'RAS4A_RAF_scfix-CYFpos.itp',
            'RAS4A_scfix-CYFpos.itp']

        for f in files_to_copy:
            if not os.path.exists(os.path.join(self.outpath, f)):
                shutil.copy(os.path.join(res_path, f), self.outpath)

        pdbs_to_tar = ['protein_chainA_ras4braf', 'protein_chainB_ras4braf', 'protein_chainAB_full_ras4braf', 'protein_chainA_ras4b', 'protein_chainA_full_ras4b',
                       'protein_chainA_ras4araf', 'protein_chainB_ras4araf', 'protein_chainAB_full_ras4araf', 'protein_chainA_ras4a', 'protein_chainA_full_ras4a', 'protein_chainA_ras4a_elastic']
        for f in pdbs_to_tar:
            _fpath = os.path.join(self.outpath, f"{f}.pdb")
            _tpath = os.path.join(self.outpath, f"{f}.tar")
            _rpath = os.path.join(res_path, f"{f}.pdb")

            if not os.path.exists(_fpath):
                shutil.copy(_rpath, self.outpath)
                with open(_fpath, 'r') as fp:
                    _d = fp.read()
                self.iointerface.save_files(_tpath, f"{f}.pdb", _d)

        # Tar archives (HB moved these here since they are supposed to be the same files)
        # ras4b-raf
        self.pdb_hvr_tar_ras4braf = os.path.join(self.outpath, 'protein_chainA_ras4braf.tar')
        self.pdb_crd_tar_ras4braf = os.path.join(self.outpath, 'protein_chainB_ras4braf.tar')
        self.pdb_full_tar_ras4braf = os.path.join(self.outpath, 'protein_chainAB_full_ras4braf.tar')
        # ras4a-raf
        self.pdb_hvr_tar_ras4araf = os.path.join(self.outpath, 'protein_chainA_ras4araf.tar')
        self.pdb_crd_tar_ras4araf = os.path.join(self.outpath, 'protein_chainB_ras4araf.tar')
        self.pdb_full_tar_ras4araf = os.path.join(self.outpath, 'protein_chainAB_full_ras4araf.tar')
        # ras4b
        self.pdb_hvr_tar_ras4b = os.path.join(self.outpath, 'protein_chainA_ras4b.tar')
        self.pdb_full_tar_ras4b = os.path.join(self.outpath, 'protein_chainA_full_ras4b.tar')
        # ras4a
        self.pdb_hvr_tar_ras4a = os.path.join(self.outpath, 'protein_chainA_ras4a.tar')
        self.pdb_full_tar_ras4a = os.path.join(self.outpath, 'protein_chainA_full_ras4a.tar')    
        # ras4a_elastic
        self.pdb_hvr_tar_ras4a_elastic = os.path.join(self.outpath, 'protein_chainA_ras4a_elastic.tar')

        # from FeedbackAA2CG
        self.revised_itp = ' '
        self.respath     = res_path
        
        # ras4b-raf 
        self.main_itp_ras4braf    = os.path.join(self.respath,'RAS_RAF_TERNARY_scfix-CYFpos.itp')
        self.itpkey_ras4braf      = os.path.join(self.outpath, 'RAS_RAF_TERNARY_scfix-CYFpos.itp')
        self.main_itp_latest_ras4braf    = os.path.join(self.outpath,'RAS_RAF_TERNARY_scfix-CYFpos.itp')
        # ras4b
        self.main_itp_ras4b = os.path.join(self.respath,'RAS_scfix-CYFpos.itp')
        self.itpkey_ras4b = os.path.join(self.outpath, 'RAS_scfix-CYFpos.itp')
        self.main_itp_latest_ras4b = os.path.join(self.outpath,'RAS_scfix-CYFpos.itp')
        # ras4a-raf 
        self.main_itp_ras4araf    = os.path.join(self.respath,'RAS4A_RAF_scfix-CYFpos.itp')
        self.itpkey_ras4araf      = os.path.join(self.outpath, 'RAS4A_RAF_scfix-CYFpos.itp')
        self.main_itp_latest_ras4araf  = os.path.join(self.outpath,'RAS4A_RAF_scfix-CYFpos.itp')
        # ras4a
        self.main_itp_ras4a = os.path.join(self.respath,'RAS4A_scfix-CYFpos.itp')
        self.itpkey_ras4a = os.path.join(self.outpath, 'RAS4A_scfix-CYFpos.itp')
        self.main_itp_latest_ras4a = os.path.join(self.outpath,'RAS4A_scfix-CYFpos.itp')
    
        # ras4b-raf 
        self.pdbkey_hvr_ras4braf  = 'protein_chainA_ras4braf.pdb'
        self.pdbkey_crd_ras4braf  = 'protein_chainB_ras4braf.pdb'
        self.pdbkey_full_ras4braf = 'protein_chainAB_full_ras4braf.pdb'
        self.main_pdb_hvr_ras4braf = os.path.join(self.outpath, self.pdbkey_hvr_ras4braf)
        self.main_pdb_crd_ras4braf = os.path.join(self.outpath, self.pdbkey_crd_ras4braf)
        self.main_pdb_full_ras4braf = os.path.join(self.outpath, self.pdbkey_full_ras4braf)
        
        # ras4b    
        self.pdbkey_hvr_ras4b  = 'protein_chainA_ras4b.pdb'
        self.pdbkey_full_ras4b = 'protein_chainA_full_ras4b.pdb'
        self.main_pdb_hvr_ras4b = os.path.join(self.outpath, self.pdbkey_hvr_ras4b)
        self.main_pdb_full_ras4b = os.path.join(self.outpath, self.pdbkey_full_ras4b)
        
        # ras4a-raf 
        self.pdbkey_hvr_ras4araf  = 'protein_chainA_ras4araf.pdb'
        self.pdbkey_crd_ras4araf  = 'protein_chainB_ras4araf.pdb'
        self.pdbkey_full_ras4araf = 'protein_chainAB_full_ras4araf.pdb'
        self.main_pdb_hvr_ras4araf = os.path.join(self.outpath, self.pdbkey_hvr_ras4araf)
        self.main_pdb_crd_ras4araf = os.path.join(self.outpath, self.pdbkey_crd_ras4araf)
        self.main_pdb_full_ras4araf = os.path.join(self.outpath, self.pdbkey_full_ras4araf)
        
        # ras4a    
        self.pdbkey_hvr_ras4a  = 'protein_chainA_ras4a.pdb'
        self.pdbkey_full_ras4a = 'protein_chainA_full_ras4a.pdb'
        self.main_pdb_hvr_ras4a = os.path.join(self.outpath, self.pdbkey_hvr_ras4a)
        self.main_pdb_full_ras4a = os.path.join(self.outpath, self.pdbkey_full_ras4a)

        # ras4a for elastic network                                                                             
        self.pdbkey_hvr_ras4a_elastic = 'protein_chainA_ras4a_elastic.pdb'
        self.main_pdb_hvr_ras4a_elastic = os.path.join(self.outpath, self.pdbkey_hvr_ras4a_elastic)

        # feedback status
        self.ras4b_feedback_done = False
        self.ras4braf_feedback_done = False
        self.ras4a_feedback_done = False
        self.ras4araf_feedback_done = False

        # tresholds for hvr and crd update
        # ras4b-raf
        self.hvr_th_ras4braf = hvr_th
        self.crd_th_ras4braf = crd_th
        # ras4b
        self.hvr_th_ras4b = hvr_th
        # ras4a-raf
        self.hvr_th_ras4araf = hvr_th
        self.crd_th_ras4araf = crd_th
        # ras4a
        self.hvr_th_ras4a = hvr_th
        
        # Aggregates
        # ras4b-raf
        self.aggregate_rmsf_ras4braf = []
        self.aggregate_coord_ras4braf = []
        # ras4b
        self.aggregate_rmsf_ras4b = []
        self.aggregate_coord_ras4b = []
        # ras4a-raf
        self.aggregate_rmsf_ras4araf = []
        self.aggregate_coord_ras4araf = []
        # ras4a
        self.aggregate_rmsf_ras4a = []
        self.aggregate_coord_ras4a = []

        # observed secondary structures
        # ras4b-raf
        self.ss_hvr_all_ras4braf = []
        self.ss_crd_all_ras4braf = []
        # ras4b
        self.ss_hvr_all_ras4b = []
        # ras4a-raf
        self.ss_hvr_all_ras4araf = []
        self.ss_crd_all_ras4araf = []
        # ras4a
        self.ss_hvr_all_ras4a = []
        
        # store pdb files associated with secondary structures
        # ras4b-raf
        self.dict_hvr_ss_pdb_ras4braf = {}
        self.dict_crd_ss_pdb_ras4braf = {}
        self.dict_full_ss_pdb_ras4braf = {}
        # ras4b
        self.dict_hvr_ss_pdb_ras4b = {}
        self.dict_full_ss_pdb_ras4b = {}
        # ras4a-raf
        self.dict_hvr_ss_pdb_ras4araf = {}
        self.dict_crd_ss_pdb_ras4araf = {}
        self.dict_full_ss_pdb_ras4araf = {}
        # ras4a
        self.dict_hvr_ss_pdb_ras4a = {}
        self.dict_full_ss_pdb_ras4a = {}

        # To get the percentage of most occurred secondary structure
        # ras4b-raf
        self.tot_num_hvr_ras4braf = 0
        self.perc_hvr_ras4braf = 0.0
        self.most_hvr_ras4braf = ' '
        self.tot_num_crd_ras4braf = 0
        self.perc_crd_ras4braf = 0.0
        self.most_crd_ras4braf = ' '
        # ras4b
        self.tot_num_hvr_ras4b = 0
        self.perc_hvr_ras4b = 0.0
        self.most_hvr_ras4b = ' '
        # ras4a-raf
        self.tot_num_hvr_ras4araf = 0
        self.perc_hvr_ras4araf = 0.0
        self.most_hvr_ras4araf = ' '
        self.tot_num_crd_ras4araf = 0
        self.perc_crd_ras4araf = 0.0
        self.most_crd_ras4araf = ' '
        # ras4a
        self.tot_num_hvr_ras4a = 0
        self.perc_hvr_ras4a = 0.0
        self.most_hvr_ras4a = ' '

        # initial secondary structures of feedback regions
        # ras4b-raf
        self.prev_hvr_ras4braf = 'HHHH2222CCCCCCCCCCCCC'  # residues 164-184
        self.prev_crd_ras4braf = 'EEEEEEECTTTCSEEEE'    # residues 279-295
        # ras4b
        self.prev_hvr_ras4b = 'HHHH2222CCCCCCCCCCCCC'  # residues 164-184
        # ras4a-raf
        self.prev_hvr_ras4araf = 'HHHH2222CCCCCCCCCCCCCC'  # residues 164-185
        self.prev_crd_ras4araf = 'EEEEEEECTTTCSEEEE'    # residues 280-296
        # ras4a
        self.prev_hvr_ras4a = 'HHHH2222CCCCCCCCCCCCCC'  # residues 164-185

        # residue number for the end of helix
        #ras4b-raf
        self.hvr_end_main_ras4braf = 171
        # ras4b
        self.hvr_end_main_ras4b = 171
        # ras4a-raf
        self.hvr_end_main_ras4araf = 171
        # ras4a
        self.hvr_end_main_ras4a = 171

        # pdb files used in the feedback
        # ras4b-raf
        self.pdb_full_crd_ras4braf = ' '
        self.pdb_full_hvr_ras4braf = ' '
        # ras4b
        self.pdb_full_hvr_ras4b = ' '
        # ras4a-raf
        self.pdb_full_crd_ras4araf = ' '
        self.pdb_full_hvr_ras4araf = ' '
        # ras4a
        self.pdb_full_hvr_ras4a = ' '

        # pdb files corresponding to the most occurred secondary structures
        # ras4b-raf
        self.pdb_most_crd_ss_ras4braf = ' '
        self.pdb_most_hvr_ss_ras4braf = ' '
        # ras4b
        self.pdb_most_hvr_ss_ras4b = ' '
        # ras4a-raf
        self.pdb_most_crd_ss_ras4araf = ' '
        self.pdb_most_hvr_ss_ras4araf = ' '
        # ras4a
        self.pdb_most_hvr_ss_ras4a = ' '


        # Feedback region info (RAS4B, RAS4B-RAF)
        # HVR
        self.main_hvr_index_start = 162   # resid - 2 (164-2=162)
        self.main_hvr_index_end_ras4b = 183    # resid - 1 (184-1=183)
        self.main_hvr_index_end_ras4a = 184    # resid - 1 (185-1=184)
        self.temp_hvr_index_start = 161   # resid - 3 (164-3=161)
        self.temp_hvr_index_end_ras4b = 182     # resid - 2 (184-2=182)
        self.temp_hvr_index_end_ras4a = 183     # resid - 2 (185-2=183)
        # CRD
        self.main_crd_index_start_ras4b = 276   # resid - 3 (279-3 = 276)
        self.main_crd_index_start_ras4a = 277   # resid - 3 (280-3 = 277)
        self.main_crd_index_end_ras4b = 293      # resid - 2 (295-2 = 293)
        self.main_crd_index_end_ras4a = 294      # resid - 2 (296-2 = 294)
        self.temp_crd_index_start_ras4b = 91    # resid - 188 (279-188 = 91)
        self.temp_crd_index_start_ras4a = 91    # resid - 188 (280-188 = 92)
        self.temp_crd_index_end_ras4b = 108      # resid - 187 (295-187 = 108)
        self.temp_crd_index_end_ras4a = 108      # resid - 187 (296-187 = 109)

        self.crd_4b_index_diff = 421
        self.crd_4a_index_diff = 425

        # Feedback frequency
        self.nframes_fb_threshold_ras4a = nframes_fb_increment
        self.nframes_fb_threshold_ras4b = nframes_fb_increment
        self.nframes_fb_threshold_ras4araf = nframes_fb_increment
        self.nframes_fb_threshold_ras4braf = nframes_fb_increment
        #self.nframes_fb_threshold = nframes_fb_increment
        
        self.nframes_fb_increment = nframes_fb_increment

        #self.restore()
        

    # --------------------------------------------------------------------------
    @staticmethod
    def load_npz(filename):
        data = np.load(filename, allow_pickle=True)
        return data['rmsf'][0], data['pdb_coord']

    # --------------------------------------------------------------------------


    # 2x types of martinize commands
    def _martinize_with_dssp(self, namespace, key, itp_name):
        os.chdir(self.workspace)
        filename = self._create_tmp(namespace, key)

        cmd = 'python martinize_2.6_3_modified.py -f {} -o temp-system.top' \
              ' -x cg.pdb -logfile martinize_log.log -p backbone -name {} -ff martini22 -nt -elastic -dssp mkdssp' \
              ' > /dev/null'
        cmd = cmd.format(filename, itp_name)
        LOGGER.debug("Executing: {}".format(cmd))
        os.system(cmd)
        if not os.path.isfile(itp_name + '.itp'):
            raise Exception('Martinize cannot create itp file: {}'.format(itp_name)) 
        os.remove(filename)

    def _martinize_with_ss(self, namespace, key, ss_key, itp_name):
        os.chdir(self.workspace)
        filename = self._create_tmp(namespace, key)

        cmd = 'python martinize_2.6_3_modified.py -f {} -o temp-system.top' \
              ' -x cg.pdb -logfile martinize_log.log -p backbone -name {} -ff martini22 -nt -elastic -ss "{}"' \
              ' > /dev/null'
        cmd = cmd.format(filename, itp_name, ss_key)
        LOGGER.debug("Executing: {}".format(cmd))
        os.system(cmd)
        if not os.path.isfile(itp_name + '.itp'):
            raise Exception('Martinize cannot create itp file: {}'.format(itp_name))
        os.remove(filename)

    def _create_tmp(self, namespace, key, outfile=''):
        data = self.iointerface.load_files(namespace, key)
        if data is None:
            err_msg = 'Failed to create tmp file: key = {}, namespace = {}'.format(key, namespace)
            LOGGER.error(err_msg)
            raise Exception(err_msg)

        if outfile == '':
            dirname, filename = os.path.split(key)
            # There shouldn't be any name collisions, but add random hash just in case
            outfile = os.path.join(dirname, "tmp_{}_{}".format(secrets.token_hex(4), filename))
        with open(outfile, 'w') as f:
            f.write(data.decode('utf-8'))
        return outfile

    # --------------------------------------------------------------------------
    @staticmethod
    def create_pdb(main_pdb, coord, region, ras_type):

        # LOGGER.info('Saving pdb file {}'.format(new_pdb))

        if region == 'hvr':
            counter = 0
        elif region == 'crd' and ras_type == 'ras4b':
            counter = 2996
        elif region == 'crd' and ras_type == 'ras4a':
            counter = 3071
        else:
            counter = 0

        data_string = ""

        if region == 'hvr' or region == 'crd':
            with open(main_pdb) as infile:

                for line in infile:
                    data = line.split()
                    if data[0] == "ATOM":
                        if len(data[2]) == 1:
                            data_string += str(data[0] + data[1].rjust(7) + "  " + data[2].rjust(1) + "   " + data[3].ljust(6) + data[4].rjust(3) + (str("{:12.3f}".format(coord[0][counter][0])).rjust(12)) + (str("{:8.3f}".format(coord[0][counter][1])).rjust(8)) + (str("{:8.3f}".format(coord[0][counter][2])).rjust(8)) + data[8].rjust(6) + data[9].rjust(6) + "\n")
                            counter = counter + 1

                        elif len(data[2]) == 2:
                            data_string += str(data[0] + data[1].rjust(7) + "  " + data[2].rjust(1) + "  " + data[3].ljust(5) + data[4].rjust(4) + (str("{:12.3f}".format(coord[0][counter][0])).rjust(12)) + (str("{:8.3f}".format(coord[0][counter][1])).rjust(8)) + (str("{:8.3f}".format(coord[0][counter][2])).rjust(8)) + data[8].rjust(6) + data[9].rjust(6) + "\n")
                            counter = counter + 1

                        elif len(data[2]) == 3:
                            data_string += str(data[0] + data[1].rjust(7) + "  " + data[2].rjust(1) + " " + data[3].ljust(4) + data[4].rjust(5) + (str("{:12.3f}".format(coord[0][counter][0])).rjust(12)) + (str("{:8.3f}".format(coord[0][counter][1])).rjust(8)) + (str("{:8.3f}".format(coord[0][counter][2])).rjust(8)) + data[8].rjust(6) + data[9].rjust(6) + "\n")
                            counter = counter + 1

                        else:
                            data_string += str(data[0] + data[1].rjust(7) + " " + data[2].ljust(2) + " " + data[3].ljust(4) + data[4].rjust(5) + (str("{:12.3f}".format(coord[0][counter][0])).rjust(12)) + (str("{:8.3f}".format(coord[0][counter][1])).rjust(8)) + (str("{:8.3f}".format(coord[0][counter][2])).rjust(8)) + data[8].rjust(6) + data[9].rjust(6) + "\n")
                            counter = counter + 1

                    else:
                        data_string += str(line)

            return data_string

        else:

            with open(main_pdb) as infile:

                for line in infile:
                    data = line.split()

                    if data[0] == "ATOM":
                        if len(data[2]) == 1:
                            data_string += str(data[0] + data[1].rjust(7) + "  " + data[2].rjust(1) + "   " + data[3].ljust(4) + (" ").ljust(0) + data[4].rjust(4) + (str("{:12.3f}".format(coord[0][counter][0])).rjust(12)) + (str("{:8.3f}".format(coord[0][counter][1])).rjust(8)) + (str("{:8.3f}".format(coord[0][counter][2])).rjust(8)) + data[8].rjust(6) + data[9].rjust(6) + (" ").rjust(12) + "\n")
                            counter = counter + 1

                        elif len(data[2]) == 2:
                            if data[2] == 'MG' or data[2] == 'ZN':
                                data_string += str(data[0] + data[1].rjust(7) + " "  + data[2] + "  " + data[3].rjust(4) + " " + (" ").ljust(0) + data[4].rjust(4) + (str("{:12.3f}".format(coord[0][counter][0])).rjust(12)) + (str("{:8.3f}".format(coord[0][counter][1])).rjust(8)) + (str("{:8.3f}".format(coord[0][counter][2])).rjust(8)) + data[8].rjust(6) + data[9].rjust(6) + (" ").rjust(12) + "\n")
                                counter = counter + 1
                            else:
                                data_string += str(data[0] + data[1].rjust(7) + "  " + data[2].rjust(1) + "  " + data[3].ljust(4) + (" ").ljust(0) + data[4].rjust(4) + (str("{:12.3f}".format(coord[0][counter][0])).rjust(12)) + (str("{:8.3f}".format(coord[0][counter][1])).rjust(8)) + (str("{:8.3f}".format(coord[0][counter][2])).rjust(8)) + data[8].rjust(6) + data[9].rjust(6) + (" ").rjust(12) + "\n")
                                counter = counter + 1

                        elif len(data[2]) == 3:
                            data_string += str(data[0] + data[1].rjust(7) + "  " + data[2].rjust(1) + " " + data[3].ljust(4) + (" ").ljust(0) + data[4].rjust(4) + (str("{:12.3f}".format(coord[0][counter][0])).rjust(12)) + (str("{:8.3f}".format(coord[0][counter][1])).rjust(8)) + (str("{:8.3f}".format(coord[0][counter][2])).rjust(8)) + data[8].rjust(6) + data[9].rjust(6) + (" ").rjust(12) + "\n")
                            counter = counter + 1

                        else:
                            data_string += str(data[0] + data[1].rjust(7) + " " + data[2].ljust(2) + " " + data[3].ljust(4) + (" ").ljust(0) + data[4].rjust(4) + (str("{:12.3f}".format(coord[0][counter][0])).rjust(12)) + (str("{:8.3f}".format(coord[0][counter][1])).rjust(8)) + (str("{:8.3f}".format(coord[0][counter][2])).rjust(8)) + data[8].rjust(6) + data[9].rjust(6) + (" ").rjust(12) + "\n")
                            counter = counter + 1

                    else:
                        data_string += str(line)

            return data_string

    # --------------------------------------------------------------------------
    def write_ss_ras4braf(self):

        ts = add_timestamp('')
        ss_output = os.path.join(self.outpath, 'ss_history_ras4braf' + ts + '.txt')

        with open(ss_output,"w") as outfile:

            outfile.write("HVR SS:" + '\n')
            for i in range(0, len(self.ss_hvr_all_ras4braf)):
                outfile.write(self.ss_hvr_all_ras4braf[i] + '\n')

            outfile.write("HVR SS counts (most common):" + '\n')
            c_hvr = Counter(self.ss_hvr_all_ras4braf)
            _most_com = c_hvr.most_common(15)
            for i in range(0, len(_most_com)):
                outfile.write(str(_most_com[i]) + '\n')

            outfile.write("CRD SS:" + '\n')
            for i in range(0, len(self.ss_crd_all_ras4braf)):
                outfile.write(self.ss_crd_all_ras4braf[i] + '\n')

            outfile.write("CRD SS counts (most common):" + '\n')
            c_crd = Counter(self.ss_crd_all_ras4braf)
            _most_com = c_crd.most_common(15)
            for i in range(0, len(_most_com)):
                outfile.write(str(_most_com[i]) + '\n')


        def _find_most_common(_list):
            _cntr = Counter(_list)
            _tot_num = sum(_cntr.values())
            _most_com = _cntr.most_common(1)[0]
            return _tot_num, _most_com[0], _most_com[1]

        self.tot_num_hvr_ras4braf, self.most_hvr_ras4braf, num_most_hvr = _find_most_common(self.ss_hvr_all_ras4braf)
        self.perc_hvr_ras4braf = num_most_hvr / self.tot_num_hvr_ras4braf
        LOGGER.info('Percentage of the most occurred HVR SS for ras4b-raf ({}) is {}'.format(self.most_hvr_ras4braf,self.perc_hvr_ras4braf))

        self.tot_num_crd_ras4braf, self.most_crd_ras4braf, num_most_crd = _find_most_common(self.ss_crd_all_ras4braf)
        self.perc_crd_ras4braf = num_most_crd / self.tot_num_crd_ras4braf
        LOGGER.info('Percentage of the most occurred CRD SS for ras4b-raf ({}) is {}'.format(self.most_crd_ras4braf,self.perc_crd_ras4braf))


        return ss_output

    # --------------------------------------------------------------------------
    
    def write_ss_ras4b(self):
        
        ts = add_timestamp('')
        ss_output = os.path.join(self.outpath, 'ss_history_ras4b' + ts + '.txt')

        with open(ss_output,"w") as outfile:

            outfile.write("HVR SS:" + '\n')
            for i in range(0, len(self.ss_hvr_all_ras4b)):
                outfile.write(self.ss_hvr_all_ras4b[i] + '\n')

            outfile.write("HVR SS counts (most common):" + '\n')
            c_hvr = Counter(self.ss_hvr_all_ras4b)
            _most_com = c_hvr.most_common(15)
            for i in range(0, len(_most_com)):
                outfile.write(str(_most_com[i]) + '\n')


        def _find_most_common(_list):
            _cntr = Counter(_list)
            _tot_num = sum(_cntr.values())
            _most_com = _cntr.most_common(1)[0]
            return _tot_num, _most_com[0], _most_com[1]


        self.tot_num_hvr_ras4b, self.most_hvr_ras4b, num_most_hvr_ras4b = _find_most_common(self.ss_hvr_all_ras4b)
        self.perc_hvr_ras4b = num_most_hvr_ras4b / self.tot_num_hvr_ras4b
        LOGGER.info('Percentage of the most occurred HVR SS for ras4b ({}) is {}'.format(self.most_hvr_ras4b,self.perc_hvr_ras4b))

        return ss_output

    # --------------------------------------------------------------------------
    
    def write_ss_ras4araf(self):

        ts = add_timestamp('')
        ss_output = os.path.join(self.outpath, 'ss_history_ras4araf' + ts + '.txt')

        with open(ss_output,"w") as outfile:

            outfile.write("HVR SS:" + '\n')
            for i in range(0, len(self.ss_hvr_all_ras4araf)):
                outfile.write(self.ss_hvr_all_ras4araf[i] + '\n')

            outfile.write("HVR SS counts (most common):" + '\n')
            c_hvr = Counter(self.ss_hvr_all_ras4araf)
            _most_com = c_hvr.most_common(15)
            for i in range(0, len(_most_com)):
                outfile.write(str(_most_com[i]) + '\n')

            outfile.write("CRD SS:" + '\n')
            for i in range(0, len(self.ss_crd_all_ras4araf)):
                outfile.write(self.ss_crd_all_ras4araf[i] + '\n')

            outfile.write("CRD SS counts (most common):" + '\n')
            c_crd = Counter(self.ss_crd_all_ras4araf)
            _most_com = c_crd.most_common(15)
            for i in range(0, len(_most_com)):
                outfile.write(str(_most_com[i]) + '\n')


        def _find_most_common(_list):
            _cntr = Counter(_list)
            _tot_num = sum(_cntr.values())
            _most_com = _cntr.most_common(1)[0]
            return _tot_num, _most_com[0], _most_com[1]

        self.tot_num_hvr_ras4araf, self.most_hvr_ras4araf, num_most_hvr = _find_most_common(self.ss_hvr_all_ras4araf)
        self.perc_hvr_ras4araf = num_most_hvr / self.tot_num_hvr_ras4araf
        LOGGER.info('Percentage of the most occurred HVR SS for ras4a-raf ({}) is {}'.format(self.most_hvr_ras4araf,self.perc_hvr_ras4araf))

        self.tot_num_crd_ras4araf, self.most_crd_ras4araf, num_most_crd = _find_most_common(self.ss_crd_all_ras4araf)
        self.perc_crd_ras4araf = num_most_crd / self.tot_num_crd_ras4araf
        LOGGER.info('Percentage of the most occurred CRD SS for ras4a-raf ({}) is {}'.format(self.most_crd_ras4araf,self.perc_crd_ras4araf))


        return ss_output

    # --------------------------------------------------------------------------
    
    def write_ss_ras4a(self):
        
        ts = add_timestamp('')
        ss_output = os.path.join(self.outpath, 'ss_history_ras4a' + ts + '.txt')

        with open(ss_output,"w") as outfile:

            outfile.write("HVR SS:" + '\n')
            for i in range(0, len(self.ss_hvr_all_ras4a)):
                outfile.write(self.ss_hvr_all_ras4a[i] + '\n')

            outfile.write("HVR SS counts (most common):" + '\n')
            c_hvr = Counter(self.ss_hvr_all_ras4a)
            _most_com = c_hvr.most_common(15)
            for i in range(0, len(_most_com)):
                outfile.write(str(_most_com[i]) + '\n')


        def _find_most_common(_list):
            _cntr = Counter(_list)
            _tot_num = sum(_cntr.values())
            _most_com = _cntr.most_common(1)[0]
            return _tot_num, _most_com[0], _most_com[1]


        self.tot_num_hvr_ras4a, self.most_hvr_ras4a, num_most_hvr_ras4a = _find_most_common(self.ss_hvr_all_ras4a)
        self.perc_hvr_ras4a = num_most_hvr_ras4a / self.tot_num_hvr_ras4a
        LOGGER.info('Percentage of the most occurred HVR SS for ras4a ({}) is {}'.format(self.most_hvr_ras4a,self.perc_hvr_ras4a))

        return ss_output

    # --------------------------------------------------------------------------

    def update_itp_hvr(self):
        LOGGER.info('Starting to update HVR region')

         
        if self.ras4braf_feedback_done == True:
            itp_type = self.main_itp_ras4braf
            print(itp_type)
        elif self.ras4b_feedback_done == True:
            itp_type = self.main_itp_ras4b
            print(itp_type)
        elif self.ras4araf_feedback_done == True:
            itp_type = self.main_itp_ras4araf
            print(itp_type)
        elif self.ras4a_feedback_done == True:
            itp_type = self.main_itp_ras4a
            print(itp_type)


        with open(self.revised_itp) as fp:

            for line in fp:
                data = line.split()

                # Find the line corresponding to secondary structure in the temporary itp file
                if len(data)>1 and data[1] == "Secondary":

                    data_ss = next(fp).split()
                    ss_new = data_ss[1]
                    LOGGER.debug('New secondary structure from temporary itp: {}'.format(ss_new))
                    #print("New secondary structure from temporary itp: ")
                    #print(ss_new)
                    #print("\n")
                    break

        #self.ss_list.append(ss_new) if ss_new not in self.ss_list else self.ss_list

        with open(itp_type) as fp:

            for line in fp:
                data = line.split()

                # Find the line corresponding to secondary structure in the main itp file
                if len(data)>1 and data[1] == "Secondary":

                    data_ss = next(fp).split()
                    ss_main = data_ss[1]
                    LOGGER.debug('Old secondary structure from main itp: {}'.format(ss_main))
                    #print("Old secondary structure from main itp: ")
                    #print(ss_main)
                    #print("\n")
                    break

        #self.ss_list.append(ss_main) if ss_main not in self.ss_list else self.ss_list

        # create lists for storing name and position of changed secondary structure
        site_changed = []
        site_position_main = []
        site_position_temp = []

        ss_hvr_revised = []

        '''
        hvr_end_main = 172

        if self.ras4braf_feedback_done == True:
            ss_hvr_end_residue = 'LYS172'
            ss_hvr_end_residue_num = '172'
            ss_hvr_end_atom_index = '379'
        elif self.ras4b_feedback_done == True:
            ss_hvr_end_residue = 'LYS172'
            ss_hvr_end_residue_num = '172'
            ss_hvr_end_atom_index = '379'
        elif self.ras4araf_feedback_done == True:
            ss_hvr_end_residue = 'SER172'
            ss_hvr_end_residue_num = '172'
            ss_hvr_end_atom_index = '381'
        elif self.ras4a_feedback_done == True:
            ss_hvr_end_residue = 'SER172'
            ss_hvr_end_residue_num = '172'
            ss_hvr_end_atom_index = '381'
        '''


        #Count difference due to the missing residue in temporary itp
        
        if self.ras4braf_feedback_done == True or self.ras4b_feedback_done == True:
            ss_main_hvr = ss_main[self.main_hvr_index_start:self.main_hvr_index_end_ras4b]
            ss_new_hvr = ss_new[self.temp_hvr_index_start:self.temp_hvr_index_end_ras4b]
        elif self.ras4araf_feedback_done == True or self.ras4a_feedback_done == True:
            ss_main_hvr = ss_main[self.main_hvr_index_start:self.main_hvr_index_end_ras4a]
            ss_new_hvr = ss_new[self.temp_hvr_index_start:self.temp_hvr_index_end_ras4a]


        if self.ras4braf_feedback_done == True or self.ras4b_feedback_done == True:

            ss_hvr_residues = {}                                                                             
            ss_hvr_residues[0] = 'ARG164'
            ss_hvr_residues[1] = 'LYS165'
            ss_hvr_residues[2] = 'HIS166'                                                                    
            ss_hvr_residues[3] = 'LYS167'                                                                     
            ss_hvr_residues[4] = 'GLU168'                                                                    
            ss_hvr_residues[5] = 'LYS169'                                                                    
            ss_hvr_residues[6] = 'MET170'                                                                    
            ss_hvr_residues[7] = 'SER171'                                                                    
            ss_hvr_residues[8] = 'LYS172'                                                                    
            ss_hvr_residues[9] = 'ASP173'                                                                    
            ss_hvr_residues[10] = 'GLY174'                                                                   
            ss_hvr_residues[11] = 'LYS175'                                                                   
            ss_hvr_residues[12] = 'LYS176'                                                                   
            ss_hvr_residues[13] = 'LYS177'                                                                   
            ss_hvr_residues[14] = 'LYS178'                                                                   
            ss_hvr_residues[15] = 'LYS179'                                                                   
            ss_hvr_residues[16] = 'LYS180'                                                                  
            ss_hvr_residues[17] = 'SER181'                                                                   
            ss_hvr_residues[18] = 'LYS182'                                                                   
            ss_hvr_residues[19] = 'THR183'                                                                   
            ss_hvr_residues[20] = 'LYS184'
 
        
            if ss_new[163]=='2' and ss_new[164]=='C':
                ss_new_hvr = ss_new[160:182]
                ss_main_hvr = ss_main[161:183]

                ss_hvr_residues = {}     
                ss_hvr_residues[0] = 'ILE163'
                ss_hvr_residues[1] = 'ARG164'
                ss_hvr_residues[2] = 'LYS165'
                ss_hvr_residues[3] = 'HIS166'
                ss_hvr_residues[4] = 'LYS167'
                ss_hvr_residues[5] = 'GLU168'
                ss_hvr_residues[6] = 'LYS169'
                ss_hvr_residues[7] = 'MET170'
                ss_hvr_residues[8] = 'SER171'
                ss_hvr_residues[9] = 'LYS172'
                ss_hvr_residues[10] = 'ASP173'
                ss_hvr_residues[11] = 'GLY174'
                ss_hvr_residues[12] = 'LYS175'
                ss_hvr_residues[13] = 'LYS176'
                ss_hvr_residues[14] = 'LYS177'
                ss_hvr_residues[15] = 'LYS178'
                ss_hvr_residues[16] = 'LYS179'
                ss_hvr_residues[17] = 'LYS180'
                ss_hvr_residues[18] = 'SER181'
                ss_hvr_residues[19] = 'LYS182'
                ss_hvr_residues[20] = 'THR183'
                ss_hvr_residues[21] = 'LYS184'
            

            elif ss_new[162]=='2' and ss_new[163]=='C':
                ss_new_hvr = ss_new[159:182]
                ss_main_hvr = ss_main[160:183]

                ss_hvr_residues = {}
                ss_hvr_residues[0] = 'GLU162'
                ss_hvr_residues[1] = 'ILE163'
                ss_hvr_residues[2] = 'ARG164'
                ss_hvr_residues[3] = 'LYS165'
                ss_hvr_residues[4] = 'HIS166'
                ss_hvr_residues[5] = 'LYS167'
                ss_hvr_residues[6] = 'GLU168'
                ss_hvr_residues[7] = 'LYS169'
                ss_hvr_residues[8] = 'MET170'
                ss_hvr_residues[9] = 'SER171'
                ss_hvr_residues[10] = 'LYS172'
                ss_hvr_residues[11] = 'ASP173'
                ss_hvr_residues[12] = 'GLY174'
                ss_hvr_residues[13] = 'LYS175'
                ss_hvr_residues[14] = 'LYS176'
                ss_hvr_residues[15] = 'LYS177'
                ss_hvr_residues[16] = 'LYS178'
                ss_hvr_residues[17] = 'LYS179'
                ss_hvr_residues[18] = 'LYS180'
                ss_hvr_residues[19] = 'SER181'
                ss_hvr_residues[20] = 'LYS182'
                ss_hvr_residues[21] = 'THR183'
                ss_hvr_residues[22] = 'LYS184'

            elif ss_new[161]=='2' and ss_new[162]=='C':
                ss_new_hvr = ss_new[158:182]                                                              
                ss_main_hvr = ss_main[159:183]     

                ss_hvr_residues = {}
                ss_hvr_residues[0] = 'ARG161'
                ss_hvr_residues[1] = 'GLU162'
                ss_hvr_residues[2] = 'ILE163'
                ss_hvr_residues[3] = 'ARG164'
                ss_hvr_residues[4] = 'LYS165'
                ss_hvr_residues[5] = 'HIS166'
                ss_hvr_residues[6] = 'LYS167'
                ss_hvr_residues[7] = 'GLU168'
                ss_hvr_residues[8] = 'LYS169'
                ss_hvr_residues[9] = 'MET170'
                ss_hvr_residues[10] = 'SER171'
                ss_hvr_residues[11] = 'LYS172'
                ss_hvr_residues[12] = 'ASP173'
                ss_hvr_residues[13] = 'GLY174'
                ss_hvr_residues[14] = 'LYS175'
                ss_hvr_residues[15] = 'LYS176'
                ss_hvr_residues[16] = 'LYS177'
                ss_hvr_residues[17] = 'LYS178'
                ss_hvr_residues[18] = 'LYS179'
                ss_hvr_residues[19] = 'LYS180'
                ss_hvr_residues[20] = 'SER181'
                ss_hvr_residues[21] = 'LYS182'
                ss_hvr_residues[22] = 'THR183'
                ss_hvr_residues[23] = 'LYS184'

        
            for i in range(0,len(ss_new_hvr)-1):   
                if (ss_new_hvr[i] == '2' and ss_new_hvr[i+1] == 'C'):
                    ss_hvr_end_residue = ss_hvr_residues[i]
                    hvr_end_main = int(ss_hvr_residues[i][-3:])
                    ss_hvr_end_residue_num = str(hvr_end_main - 2)
                    ss_hvr_end_residue_num_onebefore = str(hvr_end_main - 3)
                    ss_hvr_end_residue_num_fourbefore = str(hvr_end_main - 6)

            for i in range(0,len(ss_main_hvr)-1):
                if (ss_main_hvr[i] == '2' and ss_main_hvr[i+1] == 'C'):
                    ss_hvr_end_residue_old = ss_hvr_residues[i]
                    hvr_end_main_old = int(ss_hvr_residues[i][-3:])
                    ss_hvr_end_residue_num_old = str(hvr_end_main_old - 3)
                    ss_hvr_end_residue_num_onebefore_old = str(hvr_end_main - 4)
                    
        

        if self.ras4araf_feedback_done == True or self.ras4a_feedback_done == True:

            ss_hvr_residues = {}
            ss_hvr_residues[0] = 'ARG164'
            ss_hvr_residues[1] = 'GLN165'
            ss_hvr_residues[2] = 'TYR166'
            ss_hvr_residues[3] = 'ARG167'
            ss_hvr_residues[4] = 'LEU168'
            ss_hvr_residues[5] = 'LYS169'
            ss_hvr_residues[6] = 'LYS170'
            ss_hvr_residues[7] = 'ILE171'
            ss_hvr_residues[8] = 'SER172'
            ss_hvr_residues[9] = 'LYS173'
            ss_hvr_residues[10] = 'GLU174'
            ss_hvr_residues[11] = 'GLU175'
            ss_hvr_residues[12] = 'LYS176'
            ss_hvr_residues[13] = 'THR177'
            ss_hvr_residues[14] = 'PRO178'
            ss_hvr_residues[15] = 'GLY179'
            ss_hvr_residues[16] = 'CYP180'
            ss_hvr_residues[17] = 'VAL181'
            ss_hvr_residues[18] = 'LYS182'
            ss_hvr_residues[19] = 'ILE183'
            ss_hvr_residues[20] = 'LYS184'
            ss_hvr_residues[21] = 'LYS185'


            if ss_new[163]=='2' and ss_new[164]=='C':
                ss_new_hvr = ss_new[160:182]
                ss_main_hvr = ss_main[161:183]

                ss_hvr_residues = {}
                ss_hvr_residues[0] = 'ILE163'
                ss_hvr_residues[1] = 'ARG164'
                ss_hvr_residues[2] = 'GLN165'
                ss_hvr_residues[3] = 'TYR166'
                ss_hvr_residues[4] = 'ARG167'
                ss_hvr_residues[5] = 'LEU168'
                ss_hvr_residues[6] = 'LYS169'
                ss_hvr_residues[7] = 'LYS170'
                ss_hvr_residues[8] = 'ILE171'
                ss_hvr_residues[9] = 'SER172'
                ss_hvr_residues[10] = 'LYS173'
                ss_hvr_residues[11] = 'GLU174'
                ss_hvr_residues[12] = 'GLU175'
                ss_hvr_residues[13] = 'LYS176'
                ss_hvr_residues[14] = 'THR177'
                ss_hvr_residues[15] = 'PRO178'
                ss_hvr_residues[16] = 'GLY179'
                ss_hvr_residues[17] = 'CYP180'
                ss_hvr_residues[18] = 'VAL181'
                ss_hvr_residues[19] = 'LYS182'
                ss_hvr_residues[20] = 'ILE183'
                ss_hvr_residues[21] = 'LYS184'
                ss_hvr_residues[22] = 'LYS185'


            elif ss_new[162]=='2' and ss_new[163]=='C':
                ss_new_hvr = ss_new[159:182]
                ss_main_hvr = ss_main[160:183]

                ss_hvr_residues = {}
                ss_hvr_residues[0] = 'GLU162'
                ss_hvr_residues[1] = 'ILE163'
                ss_hvr_residues[2] = 'ARG164'
                ss_hvr_residues[3] = 'GLN165'
                ss_hvr_residues[4] = 'TYR166'
                ss_hvr_residues[5] = 'ARG167'
                ss_hvr_residues[6] = 'LEU168'
                ss_hvr_residues[7] = 'LYS169'
                ss_hvr_residues[8] = 'LYS170'
                ss_hvr_residues[9] = 'ILE171'
                ss_hvr_residues[10] = 'SER172'
                ss_hvr_residues[11] = 'LYS173'
                ss_hvr_residues[12] = 'GLU174'
                ss_hvr_residues[13] = 'GLU175'
                ss_hvr_residues[14] = 'LYS176'
                ss_hvr_residues[15] = 'THR177'
                ss_hvr_residues[16] = 'PRO178'
                ss_hvr_residues[17] = 'GLY179'
                ss_hvr_residues[18] = 'CYP180'
                ss_hvr_residues[19] = 'VAL181'
                ss_hvr_residues[20] = 'LYS182'
                ss_hvr_residues[21] = 'ILE183'
                ss_hvr_residues[22] = 'LYS184'
                ss_hvr_residues[23] = 'LYS185'


            elif ss_new[161]=='2' and ss_new[162]=='C':
                ss_new_hvr = ss_new[158:182]
                ss_main_hvr = ss_main[159:183]

                ss_hvr_residues = {}
                ss_hvr_residues[0] = 'ARG161'
                ss_hvr_residues[1] = 'GLU162'
                ss_hvr_residues[2] = 'ILE163'
                ss_hvr_residues[3] = 'ARG164'
                ss_hvr_residues[4] = 'GLN165'
                ss_hvr_residues[5] = 'TYR166'
                ss_hvr_residues[6] = 'ARG167'
                ss_hvr_residues[7] = 'LEU168'
                ss_hvr_residues[8] = 'LYS169'
                ss_hvr_residues[9] = 'LYS170'
                ss_hvr_residues[10] = 'ILE171'
                ss_hvr_residues[11] = 'SER172'
                ss_hvr_residues[12] = 'LYS173'
                ss_hvr_residues[13] = 'GLU174'
                ss_hvr_residues[14] = 'GLU175'
                ss_hvr_residues[15] = 'LYS176'
                ss_hvr_residues[16] = 'THR177'
                ss_hvr_residues[17] = 'PRO178'
                ss_hvr_residues[18] = 'GLY179'
                ss_hvr_residues[19] = 'CYP180'
                ss_hvr_residues[20] = 'VAL181'
                ss_hvr_residues[21] = 'LYS182'
                ss_hvr_residues[22] = 'ILE183'
                ss_hvr_residues[23] = 'LYS184'
                ss_hvr_residues[24] = 'LYS185'


            for i in range(0,len(ss_new_hvr)-1):
                if (ss_new_hvr[i] == '2' and ss_new_hvr[i+1] == 'C'):
                    ss_hvr_end_residue = ss_hvr_residues[i]
                    hvr_end_main = int(ss_hvr_residues[i][-3:])
                    #ss_hvr_end_residue_num = str(hvr_end_main - 2)
                    #ss_hvr_end_residue_num_onebefore = str(hvr_end_main - 3)
                    #ss_hvr_end_residue_num_fourbefore = str(hvr_end_main - 6)

            for i in range(0,len(ss_main_hvr)-1):
                if (ss_main_hvr[i] == '2' and ss_main_hvr[i+1] == 'C'):
                    ss_hvr_end_residue_old = ss_hvr_residues[i]
                    hvr_end_main_old = int(ss_hvr_residues[i][-3:])
                    #ss_hvr_end_residue_num_old = str(hvr_end_main_old - 3)
                    #ss_hvr_end_residue_num_onebefore_old = str(hvr_end_main - 4)


        ss_hvr_end_residue_num = str(hvr_end_main - 2)
        ss_hvr_end_residue_num_onebefore = str(hvr_end_main - 3)
        ss_hvr_end_residue_num_fourbefore = str(hvr_end_main - 6)

        ss_hvr_end_residue_num_old = str(hvr_end_main_old - 2)
        ss_hvr_end_residue_num_onebefore_old = str(hvr_end_main_old - 3)
        
        #print(ss_hvr_end_residue)
        #print(ss_hvr_end_residue_num)

        LOGGER.debug('Comparison between the SS of hvr region from main and temporary itp: {} and {}'.format(ss_main_hvr, ss_new_hvr))
        #print("comparison between the SS of hvr region from main and temporary itp")
        #print("main")
        #print(ss_main_hvr)
        #print("temporary")
        #print(ss_new_hvr)

        # count difference comes from the difference in numbering and missing residue (count_temp = count_main -2)
        count_main = self.main_hvr_index_start + 2
        count_temp = self.temp_hvr_index_start + 1

        if ss_new[163]=='2' and ss_new[164]=='C':
            count_main = count_main - 1
            count_temp = count_temp - 1
        elif ss_new[162]=='2' and ss_new[163]=='C':
            count_main = count_main - 2
            count_temp = count_temp - 2
        elif ss_new[161]=='2' and ss_new[162]=='C':
            count_main = count_main - 3
            count_temp = count_temp - 3


        u = zip(ss_new_hvr,ss_main_hvr)

        # compare strings for secondary structures from temporary and main itps
        # and add them to the lists. Only change the residues we want to focus (HVR)
        for i,j in u:
            #if i!=j and (i=='C' or i=='H' or i=='2'):
            if i!=j:
                site_changed.append(i)
                site_position_main.append(count_main)
                site_position_temp.append(count_temp)
                ss_hvr_revised.append(i)
            else:
                ss_hvr_revised.append(j)
            count_main+=1
            count_temp+=1

        #print("\n")

        LOGGER.debug('Revised secondary structure of the residues from hvr region: {}'.format(site_changed))
        #print("Revised secondary structure of the residues from hvr region: ")
        #print(site_changed)
        #print("\n")

        LOGGER.debug('Positions of the revised residues of the hvr in the temporary itp: {}'.format(site_position_temp))
        #print("Positions of the revised residues of the hvr in the temporary itp: ")
        #print(site_position_temp)
        #print("\n")

        LOGGER.debug('Positions of the revised residues of the hvr in the main itp: {}'.format(site_position_main))
        #print("Positions of the revised residues of the hvr in the main itp: ")
        #print(site_position_main)
        #print("\n")

        ss_hvr_revised = ''.join(ss_hvr_revised)

        if self.ras4araf_feedback_done == True or self.ras4a_feedback_done == True:
            if ss_new[163]=='2' and ss_new[164]=='C':
                ss_main_new = ss_main[0:self.main_hvr_index_start-1] + ss_hvr_revised + ss_main[self.main_hvr_index_end_ras4a:-1] + ss_main[-1]
            elif ss_new[162]=='2' and ss_new[163]=='C':
                ss_main_new = ss_main[0:self.main_hvr_index_start-2] + ss_hvr_revised + ss_main[self.main_hvr_index_end_ras4a:-1] + ss_main[-1]
            elif ss_new[161]=='2' and ss_new[162]=='C':
                ss_main_new = ss_main[0:self.main_hvr_index_start-3] + ss_hvr_revised + ss_main[self.main_hvr_index_end_ras4a:-1] + ss_main[-1]
            else:
                ss_main_new = ss_main[0:self.main_hvr_index_start] + ss_hvr_revised + ss_main[self.main_hvr_index_end_ras4a:-1] + ss_main[-1]


        if self.ras4braf_feedback_done == True or self.ras4b_feedback_done == True:
            if ss_new[163]=='2' and ss_new[164]=='C':
                ss_main_new = ss_main[0:self.main_hvr_index_start-1] + ss_hvr_revised + ss_main[self.main_hvr_index_end_ras4b:-1] + ss_main[-1]
            elif ss_new[162]=='2' and ss_new[163]=='C':
                ss_main_new = ss_main[0:self.main_hvr_index_start-2] + ss_hvr_revised + ss_main[self.main_hvr_index_end_ras4b:-1] + ss_main[-1]
            elif ss_new[161]=='2' and ss_new[162]=='C':
                ss_main_new = ss_main[0:self.main_hvr_index_start-3] + ss_hvr_revised + ss_main[self.main_hvr_index_end_ras4b:-1] + ss_main[-1]
            else:
                ss_main_new = ss_main[0:self.main_hvr_index_start] + ss_hvr_revised + ss_main[self.main_hvr_index_end_ras4b:-1] + ss_main[-1]

        
        LOGGER.debug('New secondary structure with changes in the protein region of interest: {}'.format(ss_main_new))
        #print("New secondary structure with changes in the protein region of interest: ")
        #print(ss_main_new)
        #print("\n")

        # check if there is any change in the secondary structure,
        # if not exit the function
        #if len(site_changed) == 0:
        #    print("No change in the secondary structure, exiting update_itp")
        #    return
        #else:
        #    self.multiple_itps = True

        #create lists to keep track of changed atom, angle and dihedral lines
        # from temporary, revised itp
        bb_atom_lines = []
        bb_bond_lines = []
        bb_angle_lines = []
        bb_dihedral_lines = []

        constraint_lines = []
        rubber_band_lines = []
        short_elastic_bond_lines = []
        long_elastic_bond_lines = []

        atom_indexes = []
        atom_indexes_main = []
        #dict_atm_ss = {}
        atom_indexes_elastic = []
        atom_indexes_main_elastic = []

        #find the changed lines in the temporary, revised itp file by using
        # the SS positions stored in site_position list
        with open(self.revised_itp) as fp:

            for line in fp:
                data = line.split()

                for i in range(0,len(site_position_temp)):

                    if len(data)>7 and data[2] == str(site_position_temp[i]) and data[4]=="BB":
                        bb_atom_lines.append(line)
                        atom_indexes.append(data[0])
                        atom_indexes_elastic.append(data[0])
                        print(atom_indexes)
                        #dict_atm_ss[str(int(data[0])+2)] = ss_main[(site_position_temp[i]+2)-1]
                        break

                    if len(data)>7 and data[2] == ss_hvr_end_residue_num and data[4]=="BB":
                        ss_hvr_end_atom_index = data[0]
                        # break
                    if len(data)>7 and data[2] == ss_hvr_end_residue_num_onebefore and data[4]=="BB":
                        ss_hvr_end_atom_index_onebefore = data[0]
                    if len(data)>7 and data[2] == ss_hvr_end_residue_num_fourbefore and data[4]=="BB":
                        ss_hvr_end_atom_index_fourbefore = data[0]

                    if len(data)>7 and data[2] == ss_hvr_end_residue_num_old and data[4]=="BB":
                        ss_hvr_end_atom_index_old = data[0]
                    if len(data)>7 and data[2] == ss_hvr_end_residue_num_onebefore_old and data[4]=="BB":
                        ss_hvr_end_atom_index_onebefore_old = data[0]



                # STORE THE PARAMETERS FROM DIFFERENT SECTIONS

                # BONDS
                for i in range(0,len(atom_indexes)):
                    if (len(data)==7 and data[4]=='1250' and len(data[-1])==13 and (data[0]==str(atom_indexes[i]) or data[1]==str(atom_indexes[i]))):
                        bb_bond_lines.append(line)
                        break

                # ANGLES
                    elif (len(data)==8 and len(data[7])==20 and (data[0]==str(atom_indexes[i]) or data[1]==str(atom_indexes[i]) or data[2]==str(atom_indexes[i]))):
                        bb_angle_lines.append(line)
                        break

                # DIHEDRALS
                    elif (len(data)>7 and len(data[-1])==27 and (data[0]==str(atom_indexes[i]) or data[1]==str(atom_indexes[i]) or data[2]==str(atom_indexes[i]) or data[3]==str(atom_indexes[i]))):
                        bb_dihedral_lines.append(line)
                        break

                # SHORT ELASTIC BONDS
                for i in range(0,len(atom_indexes)):
                    if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.64000" and (data[0]==str(atom_indexes[i]) or data[1]==str(atom_indexes[i]))):
                        short_elastic_bond_lines.append(line)
                        break


                # LONG ELASTIC BONDS
                for i in range(0,len(atom_indexes)):
                    if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.97000" and (data[0]==str(atom_indexes[i]) or data[1]==str(atom_indexes[i]))):
                        long_elastic_bond_lines.append(line)
                        break


                # ELASTIC BANDS
                if (len(data)==5 and len(data[-1])==18):
  
                    atom_indexes_elastic.append(ss_hvr_end_atom_index_old) if ss_hvr_end_atom_index_old not in atom_indexes_elastic else atom_indexes_elastic
                    atom_indexes_elastic.append(ss_hvr_end_atom_index_onebefore_old) if ss_hvr_end_atom_index_onebefore_old not in atom_indexes_elastic else atom_indexes_elastic
                    
                    atom_indexes_elastic.append(ss_hvr_end_atom_index_onebefore) if ss_hvr_end_atom_index_onebefore not in atom_indexes_elastic else atom_indexes_elastic

                
                
                for i in range(0,len(atom_indexes_elastic)):
                    if (len(data)==5 and len(data[-1])==18 and (data[0]==str(atom_indexes_elastic[i]) or data[1]==str(atom_indexes_elastic[i]))):

                        if (int(data[0])>=int(ss_hvr_end_atom_index) or int(data[1])>=int(ss_hvr_end_atom_index)):
                            continue

                        if (int(data[0])!=int(ss_hvr_end_atom_index_fourbefore) and int(data[1])==int(ss_hvr_end_atom_index_onebefore)):
                           
                            continue


                        rubber_band_lines.append(line)
                        break


                # CONSTRAINTS                                                                   
                for i in range(0,len(atom_indexes)):
                    if (len(data)==6 and len(data[3])==7 and data[4]==";" and (data[0]==str(atom_indexes[i]) or data[1]==str(atom_indexes[i]))):
                        constraint_lines.append(line)
                        break




        d_bb_atom_lines = deque(bb_atom_lines)
        d_bb_bond_lines = deque(bb_bond_lines)
        d_bb_angle_lines = deque(bb_angle_lines)
        d_bb_dihedral_lines = deque(bb_dihedral_lines)
        d_short_elastic_bond_lines = deque(short_elastic_bond_lines)
        d_long_elastic_bond_lines = deque(long_elastic_bond_lines)
        d_rubber_band_lines = deque(rubber_band_lines)
        d_constraint_lines = deque(constraint_lines)

        #print(d_bb_bond_lines)
        #print("ATOM INDEXES REVISED")
        #print(atom_indexes)

        
        if self.ras4braf_feedback_done == True:
            main_itp_2 = Naming.fname_fb_itp(self.itpkey_ras4braf[:-4],'temp_hvr')
        elif self.ras4b_feedback_done == True:
            main_itp_2 = Naming.fname_fb_itp(self.itpkey_ras4b[:-4],'temp_hvr')
        elif self.ras4araf_feedback_done == True:
            main_itp_2 = Naming.fname_fb_itp(self.itpkey_ras4araf[:-4],'temp_hvr')
        elif self.ras4a_feedback_done == True:
            main_itp_2 = Naming.fname_fb_itp(self.itpkey_ras4a[:-4],'temp_hvr')
             

        # write and return new itp file
        with open(itp_type) as infile, open(main_itp_2,"w") as outfile:

            for line in infile:
                line_changed = False

                data = line.split()

                if len(data)>1 and data[1] == "Secondary":
                    outfile.write(line)
                    next(infile)
                    outfile.write("; ")
                    outfile.write(ss_main_new)
                    outfile.write(" \n")
                    line_changed = True

                for i in range(0,len(site_position_main)):

                    if len(data)>6 and data[2] == str(site_position_main[i]) and data[4]=="BB":
                        line_bb = d_bb_atom_lines.popleft()
                        data_bb = line_bb.split()
                        outfile.write(str(int(data_bb[0])+2) + data_bb[1].rjust(11,' ') + str(int(data_bb[2])+2).rjust(5,' ') + data_bb[3].rjust(11,' ') + data_bb[4].rjust(8,' ') + str(int(data_bb[5])+2).rjust(5,' ') + data_bb[6][0].rjust(6,' ') + "\n")

                        atom_indexes_main.append(data[0])
                        atom_indexes_main_elastic.append(data[0])
                        line_changed = True
                        break

                    if len(data)>6 and data[2] == str(int(ss_hvr_end_residue_num)+2) and data[4]=="BB":
                        ss_hvr_end_atom_index_main = data[0]
                    if len(data)>6 and data[2] == str(int(ss_hvr_end_residue_num)+1) and data[4]=="BB":
                        ss_hvr_end_atom_index_onebefore_main = data[0]
                    if len(data)>6 and data[2] == str(int(ss_hvr_end_residue_num_old)+2) and data[4]=="BB":
                        ss_hvr_end_atom_index_main_old = data[0]
                    if len(data)>6 and data[2] == str(int(ss_hvr_end_residue_num_old)+1) and data[4]=="BB":
                        ss_hvr_end_atom_index_main_onebefore_old = data[0]
                                    
                

                if (len(data)>7 and len(data[-1])==27):
                    line_changed_dihedral = False
                    line_changed = True

                if (len(data)==7 and data[4]=='1250' and len(data[-1])==13):
                    line_changed_bond = False
                    line_changed = True

                if (len(data)==8 and len(data[7])==20):
                    line_changed_angle = False
                    line_changed = True

                if (len(data)==6 and len(data[3])==7 and data[4]==";"):
                    line_changed_constraint = False
                    line_changed = True

                
                #print("ATOM INDEXES MAIN")
                #print(atom_indexes_main)

                for i in range(0,len(atom_indexes_main)):

                    if (len(data)==7 and data[4]=='1250' and len(data[-1])==13 and (data[0]==str(atom_indexes_main[i]) or data[1]==str(atom_indexes_main[i]))):

                        line_changed_bond = True

                    elif (len(data)==8 and len(data[7])==20 and (data[0]==str(atom_indexes_main[i]) or data[1]==str(atom_indexes_main[i]) or data[2]==str(atom_indexes_main[i]))):
                        line_changed_angle = True

                    elif (len(data)>7 and len(data[-1])==27 and (data[0]==str(atom_indexes_main[i]) or data[1]==str(atom_indexes_main[i]) or data[2]==str(atom_indexes_main[i]) or data[3]==str(atom_indexes_main[i]))):

                        line_changed_dihedral = True

                    elif (len(data)==6 and len(data[3])==7 and data[4]==";" and (data[0]==str(atom_indexes_main[i]) or data[1]==str(atom_indexes_main[i]))):
                        line_changed_constraint = True



                if (len(data)>7 and len(data[-1])==27) and line_changed_dihedral == False:
                    outfile.write(line)
      
                if (len(data)==7 and data[4]=='1250' and len(data[-1])==13) and line_changed_bond == False:
                    outfile.write(line)

                if (len(data)==8 and len(data[7])==20) and line_changed_angle == False:
                    outfile.write(line)
      
                if (len(data)==6 and len(data[3])==7 and data[4]==";" and line_changed_constraint == False):
                    outfile.write(line)
      


                if len(data)==3 and data[1]=="Backbone" and data[-1] == "bonds":
                    outfile.write(line)
                    #line = next(infile)
                    for i in range(0, len(d_bb_bond_lines)):

                        line_bond = d_bb_bond_lines.popleft()
                        data_bond = line_bond.split()

                        outfile.write(str(int(data_bond[0])+2).rjust(5,' ') + str(int(data_bond[1])+2).rjust(6,' ') + data_bond[2].rjust(7,' ') + data_bond[3].rjust(10,' ') + data_bond[4].rjust(6,' ') + data_bond[5].rjust(2,' ') + data_bond[6].rjust(14,' ') + "\n")


                if len(data)==3 and data[1]=="Backbone" and data[-1] == "angles":
                    outfile.write(line)
                    #line = next(infile)
                    for i in range(0, len(d_bb_angle_lines)):

                        line_ang = d_bb_angle_lines.popleft()
                        #print(line_ang)
                        data_ang = line_ang.split()

                        if (data_ang[7][4] == 'C' or data_ang[7][11] == 'C' or data_ang[7][18] == 'C'):
                            outfile.write(str(int(data_ang[0])+2).rjust(5,' ') + str(int(data_ang[1])+2).rjust(6,' ') + str(int(data_ang[2])+2).rjust(6,' ') + data_ang[3].rjust(7,' ') + str(160).rjust(7, ' ') + str(60).rjust(6, ' ') + data_ang[6].rjust(2,' ') + data_ang[7].rjust(21,' ') + "\n")
                        else:
                            outfile.write(str(int(data_ang[0])+2).rjust(5,' ') + str(int(data_ang[1])+2).rjust(6,' ') + str(int(data_ang[2])+2).rjust(6,' ') + data_ang[3].rjust(7,' ') + str(96).rjust(7, ' ') + str(700).rjust(6, ' ') + data_ang[6].rjust(2,' ') + data_ang[7].rjust(21,' ') + "\n")



                if len(data)==3 and data[1]=="Backbone" and data[-1] == "dihedrals":
                    outfile.write(line)
                    #line = next(infile)
                    for i in range(0,len(d_bb_dihedral_lines)):

                        line_dihed = d_bb_dihedral_lines.popleft()
                        data_dihed = line_dihed.split()

                        outfile.write(str(int(data_dihed[0])+2).rjust(5,' ') + str(int(data_dihed[1])+2).rjust(6,' ') + str(int(data_dihed[2])+2).rjust(6,' ') + str(int(data_dihed[3])+2).rjust(6,' ') + data_dihed[4].rjust(7,' ') + data_dihed[5].rjust(7,' ') + data_dihed[6].rjust(6,' ') + data_dihed[7].rjust(6,' ') + data_dihed[8].rjust(2,' ') + " " + data_dihed[9].rjust(21,' ') + "\n" )

                if len(data)==3 and data[1]=="constraints":
                    outfile.write(line)
                    #line = next(infile)
                    for i in range(0,len(d_constraint_lines)):
                        line_const = d_constraint_lines.popleft()
                        data_const = line_const.split()

                        outfile.write(str(int(data_const[0])+2).rjust(5,' ') + str(int(data_const[1])+2).rjust(6,' ') + data_const[2].rjust(7,' ') + data_const[3].rjust(10,' ') + data_const[4].rjust(2,' ') + data_const[5].rjust(14,' ') + "\n")


                # First we write all the existing bonds that don't have atoms from the protein region of interest.


                #ELASTIC (RUBBER) BANDS

                if (len(data)==5 and len(data[-1])==3 and data[-1]=='500'):
                    line_changed_elastic_band = False
                    line_changed = True

                if (len(data)==5 and len(data[-1])==3 and data[-1]=='500'):
                    
                    atom_indexes_main_elastic.append(ss_hvr_end_atom_index_main_old) if ss_hvr_end_atom_index_main_old not in atom_indexes_main_elastic else atom_indexes_main_elastic
                    atom_indexes_main_elastic.append(ss_hvr_end_atom_index_main_onebefore_old) if ss_hvr_end_atom_index_main_onebefore_old not in atom_indexes_main_elastic else atom_indexes_main_elastic
                    atom_indexes_main_elastic.append(ss_hvr_end_atom_index_onebefore_main) if ss_hvr_end_atom_index_onebefore_main not in atom_indexes_main_elastic else atom_indexes_main_elastic

                    
                #print(atom_indexes_main_elastic)
                for i in range(0,len(atom_indexes_main_elastic)):
                    if (len(data)==5 and len(data[-1])==3 and data[-1]=='500' and (data[0]==str(atom_indexes_main_elastic[i]) or data[1]==str(atom_indexes_main_elastic[i]) )):
                        line_changed_elastic_band = True


                if (len(data)==5 and len(data[-1])==3 and data[-1]=='500') and line_changed_elastic_band == False:
                    if (int(data[0])>=int(ss_hvr_end_atom_index_main) or int(data[1])>=int(ss_hvr_end_atom_index_main)) and (int(data[0]) < 418 and int(data[1]) < 418):
                        continue
                    
                    outfile.write(line)


                #SHORT ELASTIC BONDS

                if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.64000"):
                    line_changed_short_elastic = False
                    line_changed = True

                for i in range(0,len(atom_indexes_main)):
                    if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.64000" and (data[0]==str(atom_indexes_main[i]) or data[1]==str(atom_indexes_main[i]))):

                        line_changed_short_elastic = True


                if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.64000") and line_changed_short_elastic == False:
                    outfile.write(line)


                #LONG ELASTIC BONDS

                if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.97000"):
                    line_changed_short_elastic = False
                    line_changed = True

                for i in range(0,len(atom_indexes_main)):
                    if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.97000" and (data[0]==str(atom_indexes_main[i]) or data[1]==str(atom_indexes_main[i]))):

                        line_changed_short_elastic = True


                if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.97000") and line_changed_short_elastic == False:
                    outfile.write(line)
                    


                # After writing the bonds outside the protein region of interest, we write all the bonds associated with the region (both atoms need to be in the region of interest) at the end of the respective sections.


                if len(data)==1 and data[0] == ";CRD":
                    for i in range(0,len(d_rubber_band_lines)):
                        #outfile.write(d_rubber_band_lines.popleft())
                        line_rubber_band = d_rubber_band_lines.popleft()
                        data_rubber_band = line_rubber_band.split()

                        outfile.write(str(int(data_rubber_band[0])+2).rjust(0,' ') + str(int(data_rubber_band[1])+2).rjust(4,' ') + str(1).rjust(2,' ') + str(data_rubber_band[3]).rjust(8,' ') + str(500).rjust(4,' ') + '\n')


                if len(data)>1 and data[1] == "Long":
                    for i in range(0,len(d_short_elastic_bond_lines)):
                        #outfile.write(d_short_elastic_bond_lines.popleft())
                        line_short_elastic = d_short_elastic_bond_lines.popleft()
                        data_short_elastic = line_short_elastic.split()

                        outfile.write(str(int(data_short_elastic[0])+2).rjust(5,' ') + str(int(data_short_elastic[1])+2).rjust(6,' ') + str(data_short_elastic[2]).rjust(7,' ') + str(data_short_elastic[3]).rjust(10,' ') + str(data_short_elastic[4]).rjust(6,' ') + str(data_short_elastic[5]).rjust(2,' ') + str(data_short_elastic[6]).rjust(14,' ') + str(data_short_elastic[7]).rjust(4,' ') + '\n')




                if len(data)>1 and data[1] == "GTP":
                    for i in range(0,len(d_long_elastic_bond_lines)):
                        #outfile.write(d_long_elastic_bond_lines.popleft())
                        line_long_elastic = d_long_elastic_bond_lines.popleft()
                        data_long_elastic = line_long_elastic.split()

                        outfile.write(str(int(data_long_elastic[0])+2).rjust(5,' ') + str(int(data_long_elastic[1])+2).rjust(6,' ') + str(data_long_elastic[2]).rjust(7,' ') + str(data_long_elastic[3]).rjust(10,' ') + str(data_long_elastic[4]).rjust(6,' ') + str(data_long_elastic[5]).rjust(2,' ') + str(data_long_elastic[6]).rjust(14,' ') + str(data_long_elastic[7]).rjust(4,' ') + '\n')

                '''
                #SCFix
                if (len(data)==10  and data[6]=="25"):
                    line_changed = True
                    line_dihed_cgfix = True

                    for i in range(0,len(atom_indexes_main)):

                        if (data[0] == atom_indexes_main[i] or data[1] == atom_indexes_main[i] or data[2] == atom_indexes_main[i] or data[3] == atom_indexes_main[i]):
                            if (dict_atm_ss[atom_indexes_main[i]] == 'H'):
                                line_dihed_cgfix = False

                    if line_dihed_cgfix == True:
                        outfile.write(line)
                '''

                
                if line_changed == False:
                    data_check = line.split()
                    if (len(data_check)==8 and len(data_check[7])==20):
                        continue
                    if (len(data_check)>7 and len(data_check[-1])==27):
                        continue
                    if (len(data_check)==6 and len(data_check[3])==7 and data_check[4]==";"):
                        continue
                    if (len(data_check)==7 and data_check[4]=='1250' and len(data_check[-1])==13):
                        continue
                    if len(data_check)==3 and data_check[1]=="Backbone" and data_check[-1] == "dihedrals":
                        continue
                    if len(data_check)==3 and data_check[1]=="Backbone" and data_check[-1] == "bonds":
                        continue
                    if len(data_check)==3 and data_check[1]=="Backbone" and data_check[-1] == "angles":
                        continue
                    if len(data_check)==3 and data_check[1]=="constraints":
                        continue
                    
                    outfile.write(line)
                    

        return main_itp_2

    # --------------------------------------------------------------------------
    def update_itp_crd(self):
        LOGGER.info('Starting to update CRD region')

        if self.ras4braf_feedback_done == True:
            itp_type = self.main_itp_ras4braf
        elif self.ras4araf_feedback_done == True:
            itp_type = self.main_itp_ras4araf


        with open(self.revised_itp) as fp:

            for line in fp:
                data = line.split()

                # Find the line corresponding to secondary structure in the temporary itp file
                if len(data)>1 and data[1] == "Secondary":

                    data_ss = next(fp).split()
                    ss_new = data_ss[1]
                    LOGGER.debug('New secondary structure from temporary itp: {}'.format(ss_new)) 
                    #print("New secondary structure: ")
                    #print(ss_new)
                    #print("\n")
                    break


        with open(itp_type) as fp:

            for line in fp:
                data = line.split()

                # Find the line corresponding to secondary structure in the main itp file
                if len(data)>1 and data[1] == "Secondary":

                    data_ss = next(fp).split()
                    ss_main = data_ss[1]
                    LOGGER.debug('Old secondary structure from main itp: {}'.format(ss_main))
                    #print("Old secondary structure from main itp: ")
                    #print(ss_main)
                    #print("\n")
                    break


        # create lists for storing name and position of changed secondary structure
        site_changed = []
        site_position_main = []
        site_position_temp = []

        ss_crd_revised = []

        #CRD = 279-295
        
        if self.ras4braf_feedback_done == True:
            ss_main_crd = ss_main[self.main_crd_index_start_ras4b:self.main_crd_index_end_ras4b]
            ss_new_crd = ss_new[self.temp_crd_index_start_ras4b:self.temp_crd_index_end_ras4b]
            count_main = self.main_crd_index_start_ras4b + 3
            count_revised = self.temp_crd_index_start_ras4b + 1
        elif self.ras4araf_feedback_done == True:
            ss_main_crd = ss_main[self.main_crd_index_start_ras4a:self.main_crd_index_end_ras4a]
            ss_new_crd = ss_new[self.temp_crd_index_start_ras4a:self.temp_crd_index_end_ras4a]
            count_main = self.main_crd_index_start_ras4a + 3
            count_revised = self.temp_crd_index_start_ras4a + 1

            
        #print(ss_main_crd)
        #print(ss_new_crd)

        #count_main = self.main_crd_index_start + 3
        #count_revised = self.temp_crd_index_start + 1

        diff_main_revised = count_main - count_revised

        u = zip(ss_new_crd,ss_main_crd)

        # compare strings for secondary structures from temporary and main itps, and add them to the lists. Only change the residues we want to focus (HVR, CRD)
        for i,j in u:
            if i!=j:
                site_changed.append(j)
                site_position_main.append(count_main)
                site_position_temp.append(count_revised)
                ss_crd_revised.append(i)
            else:
                ss_crd_revised.append(j)
            count_main+=1
            count_revised+=1


        if self.ras4braf_feedback_done == True:
            site_position_elastic_main = list(range(self.main_crd_index_start_ras4b+3, self.main_crd_index_end_ras4b+3))    # 279 to 295
            site_position_elastic_temp = list(range(self.temp_crd_index_start_ras4b+1, self.temp_crd_index_end_ras4b+1))    # 92 to 108
        elif self.ras4araf_feedback_done == True:
            site_position_elastic_main = list(range(self.main_crd_index_start_ras4a+3, self.main_crd_index_end_ras4a+3))    # 280 to 296         
            site_position_elastic_temp = list(range(self.temp_crd_index_start_ras4a+1, self.temp_crd_index_end_ras4a+1))    # 93 to 109      

       
        LOGGER.debug('Revised secondary structure of the residues from crd region: {}'.format(site_changed))
        #print("Revised secondary structure of the residues from crd region: ")
        #print(site_changed)
        #print("\n")

        LOGGER.debug('Positions of the revised residues of the crd in the temporary itp: {}'.format(site_position_temp))
        #print("Positions of the revised residues of the crd in the temporary itp: ")
        #print(site_position_temp)
        #print("\n")

        LOGGER.debug('Positions of the revised residues of the crd in the main itp: {}'.format(site_position_main))
        #print("Positions of the revised residues of the crd in the main itp: ")
        #print(site_position_main)
        #print("\n")


        ss_crd_revised = ''.join(ss_crd_revised)

        if self.ras4braf_feedback_done == True:
            ss_main_new = ss_main[0:self.main_crd_index_start_ras4b] + ss_crd_revised + ss_main[self.main_crd_index_end_ras4b:-1] + ss_main[-1]
        elif self.ras4araf_feedback_done == True:
            ss_main_new = ss_main[0:self.main_crd_index_start_ras4a] + ss_crd_revised + ss_main[self.main_crd_index_end_ras4a:-1] + ss_main[-1]


        LOGGER.debug('New secondary structure with changes in the protein region of interest: {}'.format(ss_main_new))
        #print("New secondary structure with changes in the protein region of interest: ")
        #print(ss_main_new)
        #print("\n")


        #create lists to keep track of changed atom, angle and dihedral lines from temporary, revised itp
        bb_atom_lines = []
        bb_bond_lines = []
        bb_angle_lines = []
        bb_dihedral_lines = []

        constraint_lines = []
        rubber_band_lines = []
        short_elastic_bond_lines = []
        long_elastic_bond_lines = []

        atom_indexes = []
        atom_indexes_main = []

        #atom indexes for elastic bands
        atom_indexes_elastic_temp = []
        atom_indexes_elastic_main = []

        #dict_atm_ss = {}

        #find the changed lines in the temporary, revised itp file by using the SS positions stored in site_position list
        with open(self.revised_itp) as fp:

            for line in fp:
                data = line.split()

                for i in range(0,len(site_position_temp)):

                    if len(data)>7 and data[2] == str(site_position_temp[i]) and data[4]=="BB":
                        bb_atom_lines.append(line)
                        atom_indexes.append(data[0])
                        #dict_atm_ss[str(int(data[0])+2)] = ss_main[(site_position_temp[i]+2)-1]
                        break

                #find atom indexes for elastic band update
                for i in range(0,len(site_position_elastic_temp)):
                    if len(data)>7 and data[2] == str(site_position_elastic_temp[i]) and data[4]=="BB":
                        atom_indexes_elastic_temp.append(data[0])
                        break


                # STORE THE PARAMETERS FROM DIFFERENT SECTIONS

                # BONDS
                for i in range(0,len(atom_indexes)):
                    if (len(data)==7 and data[4]=='1250' and len(data[-1])==13 and (data[0]==str(atom_indexes[i]) or data[1]==str(atom_indexes[i]))):
                        bb_bond_lines.append(line)
                        break

                # ANGLES
                    elif (len(data)==8 and len(data[7])==20 and (data[0]==str(atom_indexes[i]) or data[1]==str(atom_indexes[i]) or data[2]==str(atom_indexes[i]))):
                        bb_angle_lines.append(line)
                        break

                # DIHEDRALS
                    elif (len(data)>7 and len(data[-1])==27 and (data[0]==str(atom_indexes[i]) or data[1]==str(atom_indexes[i]) or data[2]==str(atom_indexes[i]) or data[3]==str(atom_indexes[i]))):
                        bb_dihedral_lines.append(line)
                        break

                # SHORT ELASTIC BONDS
                for i in range(0,len(atom_indexes)):
                    if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.64000" and (data[0]==str(atom_indexes[i]) or data[1]==str(atom_indexes[i]))):
                        short_elastic_bond_lines.append(line)
                        break


                # LONG ELASTIC BONDS
                for i in range(0,len(atom_indexes)):
                    if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.97000" and (data[0]==str(atom_indexes[i]) or data[1]==str(atom_indexes[i]))):
                        long_elastic_bond_lines.append(line)
                        break


                # ELASTIC BANDS
                for i in range(0,len(atom_indexes_elastic_temp)):
                    if (len(data)==5 and len(data[-1])==18 and (data[0]==str(atom_indexes_elastic_temp[i]) or data[1]==str(atom_indexes_elastic_temp[i]))):
                        rubber_band_lines.append(line)
                        break

                                        
                # CONSTRAINTS                                                         
                for i in range(0,len(atom_indexes)):
                    if (len(data)==6 and len(data[3])==7 and data[4]==";" and (data[0]==str(atom_indexes[i]) or data[1]==str(atom_indexes[i]))):
                        constraint_lines.append(line)
                        break



        #print(rubber_band_lines)
        #print(short_elastic_bond_lines)
        #print(long_elastic_bond_lines)
        #print("dihedrals from revised itp")
        #print(bb_dihedral_lines)

        #print("bonds from revised itp")
        #print(bb_bond_lines)

        d_bb_atom_lines = deque(bb_atom_lines)
        d_bb_bond_lines = deque(bb_bond_lines)
        d_bb_angle_lines = deque(bb_angle_lines)
        d_bb_dihedral_lines = deque(bb_dihedral_lines)
        d_short_elastic_bond_lines = deque(short_elastic_bond_lines)
        d_long_elastic_bond_lines = deque(long_elastic_bond_lines)
        d_rubber_band_lines = deque(rubber_band_lines)
        d_constraint_lines = deque(constraint_lines)

        
        if self.ras4braf_feedback_done == True:
            main_itp_2 = Naming.fname_fb_itp(self.itpkey_ras4braf[:-4],'temp_crd')
            crd_index_diff = self.crd_4b_index_diff
        elif self.ras4araf_feedback_done == True:
            main_itp_2 = Naming.fname_fb_itp(self.itpkey_ras4araf[:-4],'temp_crd')
            crd_index_diff = self.crd_4a_index_diff
        
        # write and return new itp file
        with open(itp_type) as infile, open(main_itp_2,"w") as outfile:

            for line in infile:
                line_changed = False

                data = line.split()

                if len(data)>1 and data[1] == "Secondary":
                    outfile.write(line)
                    next(infile)
                    outfile.write("; ")
                    outfile.write(ss_main_new)
                    outfile.write(" \n")
                    line_changed = True

                for i in range(0,len(site_position_main)):

                    if len(data)>6 and data[2] == str(site_position_main[i]) and data[4]=="BB":
                        line_bb = d_bb_atom_lines.popleft()
                        data_bb = line_bb.split()
                        outfile.write(str(int(data_bb[0])+crd_index_diff) + data_bb[1].rjust(11,' ') + str(int(data_bb[2])+diff_main_revised).rjust(5,' ') + data_bb[3].rjust(11,' ') + data_bb[4].rjust(8,' ') + str(int(data_bb[5])+crd_index_diff).rjust(5,' ') + data_bb[6][0].rjust(6,' ') + "\n")

                        atom_indexes_main.append(data[0])
                        line_changed = True
                        break

                # find atom indexes for elastic bands
                for i in range(0,len(site_position_elastic_main)):
                    
                    if len(data)>6 and data[2] == str(site_position_elastic_main[i]) and data[4]=="BB":
                        atom_indexes_elastic_main.append(data[0])
                        break
                        

                if (len(data)>7 and len(data[-1])==27):
                    line_changed_dihedral = False
                    line_changed = True

                if (len(data)==7 and data[4]=='1250' and len(data[-1])==13):
                    line_changed_bond = False
                    line_changed = True

                if (len(data)==8 and len(data[7])==20):
                    line_changed_angle = False
                    line_changed = True

                if (len(data)==6 and len(data[3])==7 and data[4]==";"):
                    line_changed_constraint = False
                    line_changed = True



                for i in range(0,len(atom_indexes_main)):

                    if (len(data)==7 and data[4]=='1250' and len(data[-1])==13 and (data[0]==str(atom_indexes_main[i]) or data[1]==str(atom_indexes_main[i]))):

                        line_changed_bond = True

                    elif (len(data)==8 and len(data[7])==20 and (data[0]==str(atom_indexes_main[i]) or data[1]==str(atom_indexes_main[i]) or data[2]==str(atom_indexes_main[i]))):

                        line_changed_angle = True

                    elif (len(data)>7 and len(data[-1])==27 and (data[0]==str(atom_indexes_main[i]) or data[1]==str(atom_indexes_main[i]) or data[2]==str(atom_indexes_main[i]) or data[3]==str(atom_indexes_main[i]))):

                        line_changed_dihedral = True

                    elif (len(data)==6 and len(data[3])==7 and data[4]==";" and (data[0]==str(atom_indexes_main[i]) or data[1]==str(atom_indexes_main[i]))):
                        line_changed_constraint = True



                if (len(data)>7 and len(data[-1])==27) and line_changed_dihedral == False:
                    outfile.write(line)

                if (len(data)==7 and data[4]=='1250' and len(data[-1])==13) and line_changed_bond == False:
                    outfile.write(line)

                if (len(data)==8 and len(data[7])==20) and line_changed_angle == False:
                    outfile.write(line)

                if (len(data)==6 and len(data[3])==7 and data[4]==";" and line_changed_constraint == False):
                    outfile.write(line)



                if len(data)==3 and data[1]=="Backbone" and data[-1] == "bonds":
                    outfile.write(line)
                    line = next(infile)
                    for i in range(0, len(d_bb_bond_lines)):

                        line_bond = d_bb_bond_lines.popleft()
                        data_bond = line_bond.split()

                        outfile.write(str(int(data_bond[0])+crd_index_diff).rjust(5,' ') + str(int(data_bond[1])+crd_index_diff).rjust(6,' ') + data_bond[2].rjust(7,' ') + data_bond[3].rjust(10,' ') + data_bond[4].rjust(6,' ') + data_bond[5].rjust(2,' ') + data_bond[6].rjust(14,' ') + "\n")


                if len(data)==3 and data[1]=="Backbone" and data[-1] == "angles":
                    outfile.write(line)
                    line = next(infile)
                    for i in range(0, len(d_bb_angle_lines)):

                        line_ang = d_bb_angle_lines.popleft()
                        data_ang = line_ang.split()

                        outfile.write(str(int(data_ang[0])+crd_index_diff).rjust(5,' ') + str(int(data_ang[1])+crd_index_diff).rjust(6,' ') + str(int(data_ang[2])+crd_index_diff).rjust(6,' ') + data_ang[3].rjust(7,' ') + str(data_ang[4]).rjust(7, ' ') + str(data_ang[5]).rjust(6, ' ') + data_ang[6].rjust(2,' ') + data_ang[7].rjust(21,' ') + "\n")



                if len(data)==3 and data[1]=="Backbone" and data[-1] == "dihedrals":
                    outfile.write(line)
                    line = next(infile)
                    for i in range(0,len(d_bb_dihedral_lines)):

                        line_dihed = d_bb_dihedral_lines.popleft()
                        data_dihed = line_dihed.split()

                        outfile.write(str(int(data_dihed[0])+crd_index_diff).rjust(5,' ') + str(int(data_dihed[1])+crd_index_diff).rjust(6,' ') + str(int(data_dihed[2])+crd_index_diff).rjust(6,' ') + str(int(data_dihed[3])+crd_index_diff).rjust(6,' ') + data_dihed[4].rjust(7,' ') + data_dihed[5].rjust(7,' ') + data_dihed[6].rjust(6,' ') + data_dihed[7].rjust(6,' ') + data_dihed[8].rjust(2,' ') + " " + data_dihed[9].rjust(21,' ') + "\n" )


                if len(data)==3 and data[1]=="constraints":
                    outfile.write(line)
                    line = next(infile)
                    for i in range(0,len(d_constraint_lines)):
                        line_const = d_constraint_lines.popleft()
                        data_const = line_const.split()

                        outfile.write(str(int(data_const[0])+crd_index_diff).rjust(5,' ') + str(int(data_const[1])+crd_index_diff).rjust(6,' ') + data_const[2].rjust(7,' ') + data_const[3].rjust(10,' ') + data_const[4].rjust(2,' ') + data_const[5].rjust(14,' ') + "\n")


                # First we write all the existing bonds that don't have atoms from the protein region of interest.


                #ELASTIC (RUBBER) BANDS

                if (len(data)==5 and len(data[-1])==3 and data[-1]=='500'):
                    line_changed_elastic_band = False
                    line_changed = True

                for i in range(0,len(atom_indexes_elastic_main)):
                    if (len(data)==5 and len(data[-1])==3 and data[-1]=='500' and (data[0]==str(atom_indexes_elastic_main[i]) or data[1]==str(atom_indexes_elastic_main[i]))):

                        line_changed_elastic_band = True


                if (len(data)==5 and len(data[-1])==3 and data[-1]=='500') and line_changed_elastic_band == False:
                    outfile.write(line)


                #SHORT ELASTIC BONDS

                if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.64000"):
                    line_changed_short_elastic = False
                    line_changed = True

                for i in range(0,len(atom_indexes_main)):
                    if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.64000" and (data[0]==str(atom_indexes_main[i]) or data[1]==str(atom_indexes_main[i]))):

                        line_changed_short_elastic = True


                if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.64000") and line_changed_short_elastic == False:
                    outfile.write(line)


                #LONG ELASTIC BONDS

                if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.97000"):
                    line_changed_short_elastic = False
                    line_changed = True

                for i in range(0,len(atom_indexes_main)):
                    if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.97000" and (data[0]==str(atom_indexes_main[i]) or data[1]==str(atom_indexes_main[i]))):

                        line_changed_short_elastic = True


                if (len(data)==8 and len(data[-1])==3 and len(data[-3])==1 and data[3]=="0.97000") and line_changed_short_elastic == False:
                    outfile.write(line)



                # After writing the bonds outside the protein region of interest, we write all the bonds associated with the region (both atoms need to be in the region of interest) at the end of the respective sections.


                if len(data)==1 and data[0] == ";CRD":
                    outfile.write(line)
                    line = next(infile)
                    for i in range(0,len(d_rubber_band_lines)):
                        #outfile.write(d_rubber_band_lines.popleft())
                        line_rubber_band = d_rubber_band_lines.popleft()
                        data_rubber_band = line_rubber_band.split()

                        outfile.write(str(int(data_rubber_band[0])+crd_index_diff).rjust(0,' ') + str(int(data_rubber_band[1])+crd_index_diff).rjust(4,' ') + str(1).rjust(2,' ') + str(data_rubber_band[3]).rjust(8,' ') + str(500).rjust(4,' ') + '\n')


                if len(data)>1 and data[1] == "Long":
                    for i in range(0,len(d_short_elastic_bond_lines)):
                        #outfile.write(d_short_elastic_bond_lines.popleft())
                        line_short_elastic = d_short_elastic_bond_lines.popleft()
                        data_short_elastic = line_short_elastic.split()

                        outfile.write(str(int(data_short_elastic[0])+crd_index_diff).rjust(5,' ') + str(int(data_short_elastic[1])+crd_index_diff).rjust(6,' ') + str(data_short_elastic[2]).rjust(7,' ') + str(data_short_elastic[3]).rjust(10,' ') + str(data_short_elastic[4]).rjust(6,' ') + str(data_short_elastic[5]).rjust(2,' ') + str(data_short_elastic[6]).rjust(14,' ') + str(data_short_elastic[7]).rjust(4,' ') + '\n')




                if len(data)>1 and data[1] == "GTP":
                    for i in range(0,len(d_long_elastic_bond_lines)):
                        #outfile.write(d_long_elastic_bond_lines.popleft())
                        line_long_elastic = d_long_elastic_bond_lines.popleft()
                        data_long_elastic = line_long_elastic.split()

                        outfile.write(str(int(data_long_elastic[0])+crd_index_diff).rjust(5,' ') + str(int(data_long_elastic[1])+crd_index_diff).rjust(6,' ') + str(data_long_elastic[2]).rjust(7,' ') + str(data_long_elastic[3]).rjust(10,' ') + str(data_long_elastic[4]).rjust(6,' ') + str(data_long_elastic[5]).rjust(2,' ') + str(data_long_elastic[6]).rjust(14,' ') + str(data_long_elastic[7]).rjust(4,' ') + '\n')

                '''
                #SCFix
                if (len(data)==10  and data[6]=="25"):
                    line_changed = True
                    line_dihed_cgfix = True

                    for i in range(0,len(atom_indexes_main)):

                        if (data[0] == atom_indexes_main[i] or data[1] == atom_indexes_main[i] or data[2] == atom_indexes_main[i] or data[3] == atom_indexes_main[i]):
                            if (dict_atm_ss[atom_indexes_main[i]] == 'H'):
                                line_dihed_cgfix = False

                    if line_dihed_cgfix == True:
                        outfile.write(line)
                '''



                if line_changed == False:
                    outfile.write(line)



        return main_itp_2

    # --------------------------------------------------------------------------
    # load the data (do we need this as a separate function?)
    def load(self):
        print('AAtoCGFeedback_Worker.load() is not implemented')
        return

    # --------------------------------------------------------------------------
    # aggregate the data
    def aggregate(self):

        # ----------------------------------------------------------------------
        LOGGER.debug(f'Aggregating data ({self.name}): naggregates = {self.naggregates}')

        # read all the keys from the feedback-aa directory
        # keys for ras4b-raf
        fb_key_ras4braf = Naming.fb_aa_key_ras4braf('worker')
        keys_ras4braf = self.iointerface.list_keys(self.fbpath, fb_key_ras4braf, interfaces=['simple'])
        # keys for ras4b
        fb_key_ras4b = Naming.fb_aa_key_ras4b('worker')
        keys_ras4b = self.iointerface.list_keys(self.fbpath, fb_key_ras4b, interfaces=['simple'])
        # keys for ras4a-raf
        fb_key_ras4araf = Naming.fb_aa_key_ras4araf('worker')
        keys_ras4araf = self.iointerface.list_keys(self.fbpath, fb_key_ras4araf, interfaces=['simple'])
        # keys for ras4a
        fb_key_ras4a = Naming.fb_aa_key_ras4a('worker')
        keys_ras4a = self.iointerface.list_keys(self.fbpath, fb_key_ras4a, interfaces=['simple'])


        # TODO: fix this rht right way
        # We should revisit the use of the 'done' namespace, because it requires
        # explicit filtering over the whole dataset
        
        # ras4b-raf
        keys_ras4braf = [_ for _ in keys_ras4braf if not _.startswith('done')]
        nkeys_ras4braf = len(keys_ras4braf)
        # ras4b
        keys_ras4b = [_ for _ in keys_ras4b if not _.startswith('done')]
        nkeys_ras4b = len(keys_ras4b)
        # ras4a-raf
        keys_ras4araf = [_ for _ in keys_ras4araf if not _.startswith('done')]
        nkeys_ras4araf = len(keys_ras4araf)
        # ras4a
        keys_ras4a = [_ for _ in keys_ras4a if not _.startswith('done')]
        nkeys_ras4a = len(keys_ras4a)


        LOGGER.debug('found {} keys for ras4b-raf in ({})'.format(nkeys_ras4braf, self.fbpath))
        LOGGER.debug('found {} keys for ras4b in ({})'.format(nkeys_ras4b, self.fbpath))
        LOGGER.debug('found {} keys for ras4a-raf in ({})'.format(nkeys_ras4araf, self.fbpath))
        LOGGER.debug('found {} keys for ras4a in ({})'.format(nkeys_ras4a, self.fbpath))
        
        if nkeys_ras4braf == 0 and nkeys_ras4b == 0 and nkeys_ras4araf == 0 and nkeys_ras4a == 0:
            return

        # RAS4B-RAF
        if nkeys_ras4braf != 0:
            # for all keys, read the rmsf and save to tarfile
            tarfile_ras4braf = os.path.join(self.fbpath, 'done_ras4braf.tar')

            # ------------------------------------------------------------------
            # 2021/01/04: specialized code to handle the processing faster
            # assuming the interfaces to be simple and tar
            LOGGER.debug('loading {} feedback files for ras4b-raf'.format(len(keys_ras4braf)))
            _data_ras4braf = []
            valid_keys_ras4braf = []
            files_ras4braf = []
            for k in keys_ras4braf:
                npz = f"{os.path.join(self.fbpath, k)}"
                try:
                    _ = FeedbackManager_AA2CG.load_npz(npz)
                    _data_ras4braf.append(_)
                    valid_keys_ras4braf.append(k)
                    files_ras4braf.append(npz)
                except Exception as _excp:
                    LOGGER.warning('Failed to load RMSF + coordinates from ({}): {}'.format(npz, _excp))
                    continue

            rmsfs_ras4braf = [_[0] for _ in _data_ras4braf]
            coords_ras4braf = [_[1] for _ in _data_ras4braf]
            for i in range(len(_data_ras4braf)):
                np.append(self.aggregate_rmsf_ras4braf, rmsfs_ras4braf[i])
                np.append(self.aggregate_coord_ras4braf, coords_ras4braf[i])


            # ------------------------------------------------------------------
            # parallel execution of check_ss
            t0 = time.time()

            pool = multiprocessing.Pool(
                    processes=8,
                    initializer=sig_ign_and_rename_proc,
                    initargs=('pool_check_ss_ras4braf',)
                )
            results_ras4braf = pool.map(self.check_ss_ras4braf, coords_ras4braf)
            pool.close()
            
            # collect the files to write
            pdbs_hvr_to_tar_ras4braf = []
            pdbs_crd_to_tar_ras4braf = []
            pdbs_full_to_tar_ras4braf = []

            _append_idx_to_file_ras4braf = lambda _,i: _[:-4] + '-{}'.format(i) + _[-4:]
            idx = 0
            for r in results_ras4braf:
                idx += 1
                _pdb_hvr_ras4braf = _append_idx_to_file_ras4braf(r['pdb_hvr'], idx)
                _pdb_crd_ras4braf = _append_idx_to_file_ras4braf(r['pdb_crd'], idx)
                _pdb_full_ras4braf = _append_idx_to_file_ras4braf(r['pdb_full'], idx)

                pdbs_hvr_to_tar_ras4braf.append( (_pdb_hvr_ras4braf, r['pdb_hvr_data']) )
                pdbs_crd_to_tar_ras4braf.append( (_pdb_crd_ras4braf, r['pdb_crd_data']) )
                pdbs_full_to_tar_ras4braf.append( (_pdb_full_ras4braf, r['pdb_full_data']) )

                ss_hvr_fixed = r['ss_hvr_fixed']
                ss_crd = r['ss_crd']
                if r['do_hvr_append']:
                    self.ss_hvr_all_ras4braf.append(ss_hvr_fixed)
                    self.dict_hvr_ss_pdb_ras4braf[ss_hvr_fixed] = _pdb_hvr_ras4braf
                    self.dict_full_ss_pdb_ras4braf[ss_hvr_fixed] = _pdb_full_ras4braf
                if r['do_crd_append']:
                    self.ss_crd_all_ras4braf.append(ss_crd)
                    self.dict_crd_ss_pdb_ras4braf[ss_crd] = _pdb_crd_ras4braf
                    self.dict_full_ss_pdb_ras4braf[ss_crd] = _pdb_full_ras4braf

            LOGGER.info('Completed parallel execution of check_ss_ras4braf in {} sec'.format(time.time()-t0))

            # ------------------------------------------------------------------
            # need to tar up all the pdbs
            tario = mummi_core.get_io('taridx')
            FeedbackManager_AA2CG.save_tuple(tario, self.pdb_hvr_tar_ras4braf, pdbs_hvr_to_tar_ras4braf)
            FeedbackManager_AA2CG.save_tuple(tario, self.pdb_crd_tar_ras4braf, pdbs_crd_to_tar_ras4braf)
            FeedbackManager_AA2CG.save_tuple(tario, self.pdb_full_tar_ras4braf, pdbs_full_to_tar_ras4braf)

            # now delete all the data that has been processed
            tario.save_rmsfs_all(tarfile_ras4braf, valid_keys_ras4braf, rmsfs_ras4braf, coords_ras4braf)
            LOGGER.debug('END AGGREGATION: removing %d files of %d keys', len(files_ras4braf), len(keys_ras4braf))
            remove_parallel(files_ras4braf)
            # ------------------------------------------------------------------
            self.write_ss_ras4braf()
            self.naggregates += 1

        # RAS4B
        if nkeys_ras4b != 0:

            # ras4b
            tarfile_ras4b = os.path.join(self.fbpath, 'done_ras4b.tar')

            # same process for ras4b
            LOGGER.debug('loading {} feedback files for ras4b'.format(len(keys_ras4b)))
            _data_ras4b = []
            valid_keys_ras4b = []
            files_ras4b = []
            for k in keys_ras4b:
                npz = f"{os.path.join(self.fbpath, k)}"
                try:
                    _ = FeedbackManager_AA2CG.load_npz(npz)
                    _data_ras4b.append(_)
                    valid_keys_ras4b.append(k)
                    files_ras4b.append(npz)
                except Exception as _excp:
                    LOGGER.warning('Failed to load RMSF + coordinates from ({}): {}'.format(npz, _excp))
                    continue

            rmsfs_ras4b = [_[0] for _ in _data_ras4b]
            coords_ras4b = [_[1] for _ in _data_ras4b]
            for i in range(len(_data_ras4b)):
                np.append(self.aggregate_rmsf_ras4b, rmsfs_ras4b[i])
                np.append(self.aggregate_coord_ras4b, coords_ras4b[i])


            # ------------------------------------------------------------------
            # parallel execution of check_ss_ras4b
            t0 = time.time()
            if False:
               pool = multiprocessing.Pool(
                processes=8,
                initializer=sig_ign_and_rename_proc,
                initargs=('pool_check_ss_ras4b',)
                )
               results_ras4b = pool.map(self.check_ss_ras4b, coords_ras4b)
               pool.close()
            else:
               results_ras4b = [self.check_ss_ras4b(s) for s in coords_ras4b]


            # collect the files to write for ras4b (no crd region)
            pdbs_hvr_to_tar_ras4b = []
            pdbs_full_to_tar_ras4b = []

            _append_idx_to_file_ras4b = lambda _,i: _[:-4] + '-{}'.format(i) + _[-4:]
            idx = 0
            for r in results_ras4b:
                idx += 1
                _pdb_hvr_ras4b = _append_idx_to_file_ras4b(r['pdb_hvr'], idx)
                _pdb_full_ras4b = _append_idx_to_file_ras4b(r['pdb_full'], idx)

                pdbs_hvr_to_tar_ras4b.append( (_pdb_hvr_ras4b, r['pdb_hvr_data']) )
                pdbs_full_to_tar_ras4b.append( (_pdb_full_ras4b, r['pdb_full_data']) )

                ss_hvr_fixed = r['ss_hvr_fixed']
                if r['do_hvr_append']:
                    self.ss_hvr_all_ras4b.append(ss_hvr_fixed)
                    self.dict_hvr_ss_pdb_ras4b[ss_hvr_fixed] = _pdb_hvr_ras4b
                    self.dict_full_ss_pdb_ras4b[ss_hvr_fixed] = _pdb_full_ras4b


            LOGGER.info('Completed parallel execution of check_ss_ras4b in {} sec'.format(time.time()-t0))

            # ------------------------------------------------------------------
            # need to tar up all the pdbs for ras4b
            tario = mummi_core.get_io('taridx')
            FeedbackManager_AA2CG.save_tuple(tario, self.pdb_hvr_tar_ras4b, pdbs_hvr_to_tar_ras4b)
            FeedbackManager_AA2CG.save_tuple(tario, self.pdb_full_tar_ras4b, pdbs_full_to_tar_ras4b)

            # now delete all the data that has been processed for ras4b
            tario.save_rmsfs_all(tarfile_ras4b, valid_keys_ras4b, rmsfs_ras4b, coords_ras4b)
            LOGGER.debug('END AGGREGATION: removing %d files of %d keys for ras4b', len(files_ras4b), len(keys_ras4b))
            remove_parallel(files_ras4b)
            # ------------------------------------------------------------------
 
            self.write_ss_ras4b()
            self.naggregates += 1
            
        # RAS4A-RAF    
        if nkeys_ras4araf != 0:
            # for all keys, read the rmsf and save to tarfile
            tarfile_ras4araf = os.path.join(self.fbpath, 'done_ras4araf.tar')

            # ------------------------------------------------------------------
            # 2021/01/04: specialized code to handle the processing faster
            # assuming the interfaces to be simple and tar
            LOGGER.debug('loading {} feedback files for ras4a-raf'.format(len(keys_ras4araf)))
            _data_ras4araf = []
            valid_keys_ras4araf = []
            files_ras4araf = []
            for k in keys_ras4araf:
                npz = f"{os.path.join(self.fbpath, k)}"
                try:
                    _ = FeedbackManager_AA2CG.load_npz(npz)
                    _data_ras4araf.append(_)
                    valid_keys_ras4araf.append(k)
                    files_ras4araf.append(npz)
                except Exception as _excp:
                    LOGGER.warning('Failed to load RMSF + coordinates from ({}): {}'.format(npz, _excp))
                    continue

            rmsfs_ras4araf = [_[0] for _ in _data_ras4araf]
            coords_ras4araf = [_[1] for _ in _data_ras4araf]
            for i in range(len(_data_ras4araf)):
                np.append(self.aggregate_rmsf_ras4araf, rmsfs_ras4araf[i])
                np.append(self.aggregate_coord_ras4araf, coords_ras4araf[i])


            # ------------------------------------------------------------------
            # parallel execution of check_ss
            t0 = time.time()

            pool = multiprocessing.Pool(
                    processes=8,
                    initializer=sig_ign_and_rename_proc,
                    initargs=('pool_check_ss_ras4araf',)
            )

            results_ras4araf = pool.map(self.check_ss_ras4araf, coords_ras4araf)
            pool.close()
        
            
            # collect the files to write
            pdbs_hvr_to_tar_ras4araf = []
            pdbs_crd_to_tar_ras4araf = []
            pdbs_full_to_tar_ras4araf = []

            _append_idx_to_file_ras4araf = lambda _,i: _[:-4] + '-{}'.format(i) + _[-4:]
            idx = 0
            for r in results_ras4araf:
                idx += 1
                _pdb_hvr_ras4araf = _append_idx_to_file_ras4araf(r['pdb_hvr'], idx)
                _pdb_crd_ras4araf = _append_idx_to_file_ras4araf(r['pdb_crd'], idx)
                _pdb_full_ras4araf = _append_idx_to_file_ras4araf(r['pdb_full'], idx)

                pdbs_hvr_to_tar_ras4araf.append( (_pdb_hvr_ras4araf, r['pdb_hvr_data']) )
                pdbs_crd_to_tar_ras4araf.append( (_pdb_crd_ras4araf, r['pdb_crd_data']) )
                pdbs_full_to_tar_ras4araf.append( (_pdb_full_ras4araf, r['pdb_full_data']) )

                ss_hvr_fixed = r['ss_hvr_fixed']
                ss_crd = r['ss_crd']
                if r['do_hvr_append']:
                    self.ss_hvr_all_ras4araf.append(ss_hvr_fixed)
                    self.dict_hvr_ss_pdb_ras4araf[ss_hvr_fixed] = _pdb_hvr_ras4araf
                    self.dict_full_ss_pdb_ras4araf[ss_hvr_fixed] = _pdb_full_ras4araf
                if r['do_crd_append']:
                    self.ss_crd_all_ras4araf.append(ss_crd)
                    self.dict_crd_ss_pdb_ras4araf[ss_crd] = _pdb_crd_ras4araf
                    self.dict_full_ss_pdb_ras4araf[ss_crd] = _pdb_full_ras4araf

            LOGGER.info('Completed parallel execution of check_ss_ras4araf in {} sec'.format(time.time()-t0))

            # ------------------------------------------------------------------
            # need to tar up all the pdbs
            tario = mummi_core.get_io('taridx')
            FeedbackManager_AA2CG.save_tuple(tario, self.pdb_hvr_tar_ras4araf, pdbs_hvr_to_tar_ras4araf)
            FeedbackManager_AA2CG.save_tuple(tario, self.pdb_crd_tar_ras4araf, pdbs_crd_to_tar_ras4araf)
            FeedbackManager_AA2CG.save_tuple(tario, self.pdb_full_tar_ras4araf, pdbs_full_to_tar_ras4araf)

            # now delete all the data that has been processed
            tario.save_rmsfs_all(tarfile_ras4araf, valid_keys_ras4araf, rmsfs_ras4araf, coords_ras4araf)
            LOGGER.debug('END AGGREGATION: removing %d files of %d keys', len(files_ras4araf), len(keys_ras4araf))
            remove_parallel(files_ras4araf)
            # ------------------------------------------------------------------
            self.write_ss_ras4araf()
            self.naggregates += 1

        # RAS4A
        if nkeys_ras4a != 0:

            # ras4a
            tarfile_ras4a = os.path.join(self.fbpath, 'done_ras4a.tar')

            # same process for ras4a
            LOGGER.debug('loading {} feedback files for ras4a'.format(len(keys_ras4a)))
            _data_ras4a = []
            valid_keys_ras4a = []
            files_ras4a = []
            for k in keys_ras4a:
                npz = f"{os.path.join(self.fbpath, k)}"
                try:
                    _ = FeedbackManager_AA2CG.load_npz(npz)
                    _data_ras4a.append(_)
                    valid_keys_ras4a.append(k)
                    files_ras4a.append(npz)
                except Exception as _excp:
                    LOGGER.error('Failed to load RMSF + coordinates from ({}): {}'.format(npz, _excp))
                    continue

            rmsfs_ras4a = [_[0] for _ in _data_ras4a]
            coords_ras4a = [_[1] for _ in _data_ras4a]
            for i in range(len(_data_ras4a)):
                np.append(self.aggregate_rmsf_ras4a, rmsfs_ras4a[i])
                np.append(self.aggregate_coord_ras4a, coords_ras4a[i])


            # ------------------------------------------------------------------
            # parallel execution of check_ss_ras4b
            t0 = time.time()
            pool = multiprocessing.Pool(
                    processes=8,
                    initializer=sig_ign_and_rename_proc,
                    initargs=('pool_check_ss_ras4a',)
            )

            results_ras4a = pool.map(self.check_ss_ras4a, coords_ras4a)
            pool.close()


            # collect the files to write for ras4a (no crd region)
            pdbs_hvr_to_tar_ras4a = []
            pdbs_full_to_tar_ras4a = []

            _append_idx_to_file_ras4a = lambda _,i: _[:-4] + '-{}'.format(i) + _[-4:]
            idx = 0
            for r in results_ras4a:
                idx += 1
                _pdb_hvr_ras4a = _append_idx_to_file_ras4a(r['pdb_hvr'], idx)
                _pdb_full_ras4a = _append_idx_to_file_ras4a(r['pdb_full'], idx)

                pdbs_hvr_to_tar_ras4a.append( (_pdb_hvr_ras4a, r['pdb_hvr_data']) )
                pdbs_full_to_tar_ras4a.append( (_pdb_full_ras4a, r['pdb_full_data']) )

                ss_hvr_fixed = r['ss_hvr_fixed']
                if r['do_hvr_append']:
                    self.ss_hvr_all_ras4a.append(ss_hvr_fixed)
                    self.dict_hvr_ss_pdb_ras4a[ss_hvr_fixed] = _pdb_hvr_ras4a
                    self.dict_full_ss_pdb_ras4a[ss_hvr_fixed] = _pdb_full_ras4a


            LOGGER.info('Completed parallel execution of check_ss_ras4a in {} sec'.format(time.time()-t0))

            # ------------------------------------------------------------------
            # need to tar up all the pdbs for ras4a
            tario = mummi_core.get_io('taridx')
            FeedbackManager_AA2CG.save_tuple(tario, self.pdb_hvr_tar_ras4a, pdbs_hvr_to_tar_ras4a)
            FeedbackManager_AA2CG.save_tuple(tario, self.pdb_full_tar_ras4a, pdbs_full_to_tar_ras4a)

            # now delete all the data that has been processed for ras4a
            tario.save_rmsfs_all(tarfile_ras4a, valid_keys_ras4a, rmsfs_ras4a, coords_ras4a)
            LOGGER.debug('END AGGREGATION: removing %d files of %d keys for ras4a', len(files_ras4a), len(keys_ras4a))
            remove_parallel(files_ras4a)
            # ------------------------------------------------------------------
 
            self.write_ss_ras4a()
            self.naggregates += 1    

        LOGGER.info('Aggregated data ({}): naggregates = {}, naggs_at_last_report = {}'.format(self.name, self.naggregates, self.naggs_at_last_report))

    # --------------------------------------------------------------------------
    @staticmethod
    def _extract_ss(_):
        with open(_) as fp:
            for line in fp:
                data = line.split()
                if len(data) > 1 and data[1] == "Secondary":
                    data_ss = next(fp).split()
                    ss = data_ss[1]
                    # print(ss)
                    return ss
        return ''

    def check_ss_ras4braf(self, coord_ras4braf):

        # ----------------------------------------------------------------------
        cmd = 'python martinize_2.6_3_modified.py -f {} -o temp-system.top' \
              ' -x cg.pdb -logfile martinize_log.log -p backbone  -ff martini22 -nt -elastic -dssp mkdssp' \
              ' > /dev/null'


        region1 = 'hvr'
        region2 = 'crd'
        region3 = 'full'

        ss_hvr_fixed = ''
        ss_crd = ''
        do_hvr_append = False
        do_crd_append = False

        def _write(_fname, _data):
            with open(_fname, 'w') as fp:
                fp.write(_data)

        # ----------------------------------------------------------------------
        p = multiprocessing.current_process()

        # create pdb strings using the coords
        pdb_hvr_data = FeedbackManager_AA2CG.create_pdb(self.main_pdb_hvr_ras4braf, coord_ras4braf, region1, 'ras4b')
        pdb_crd_data = FeedbackManager_AA2CG.create_pdb(self.main_pdb_crd_ras4braf, coord_ras4braf, region2, 'ras4b')
        pdb_full_data = FeedbackManager_AA2CG.create_pdb(self.main_pdb_full_ras4braf, coord_ras4braf, region3, 'ras4b')

        # ---------------------------------------------------------------------
        # these are the counts before the current process pool started
        _ts = add_timestamp('')[1:]
        pdbid_hvr = '{}-p{}'.format(_ts, p.pid)
        pdbid_crd = '{}-p{}'.format(_ts, p.pid)
        pdbid_full = '{}-p{}'.format(_ts, p.pid)

        pdb_hvr = os.path.basename(Naming.fname_fb_pdb(self.pdbkey_hvr_ras4braf[:-4], region1, pdbid_hvr))
        pdb_crd = os.path.basename(Naming.fname_fb_pdb(self.pdbkey_crd_ras4braf[:-4], region2, pdbid_crd))
        pdb_full = os.path.basename(Naming.fname_fb_pdb(self.pdbkey_full_ras4braf[:-9], region3, pdbid_full))

        # ---------------------------------------------------------------------
        # create a temp directory to allow concurrent runs
        tmp_dir = os.path.join(self.workspace, 'check_ss_ras4braf_{}'.format(p.pid))
        #LOGGER.debug('Creating temp directory for check_ss_ras4braf ({})'.format(tmp_dir))
        mkdir(tmp_dir)
        os.chdir(tmp_dir)

        _write(os.path.join(tmp_dir, pdb_hvr), pdb_hvr_data)
        _write(os.path.join(tmp_dir, pdb_crd), pdb_crd_data)
        _write(os.path.join(tmp_dir, pdb_full), pdb_full_data)

        # copy relevant files to the directory
        files_to_copy = ['martinize_2.6_3_modified.py']
        for _ in files_to_copy:
           shutil.copy(os.path.join(self.workspace, _), tmp_dir)

        _write(os.path.join(tmp_dir, pdb_hvr), pdb_hvr_data)
        _write(os.path.join(tmp_dir, pdb_crd), pdb_crd_data)
        _write(os.path.join(tmp_dir, pdb_full), pdb_full_data)

        # -----------------------------------------------------------------
        # execute first martinize
        #LOGGER.debug('executing ('+cmd.format(pdb_hvr)+')')
        os.system(cmd.format(pdb_hvr))

        # fetch the secondary structure
        ss = FeedbackManager_AA2CG._extract_ss('Protein.itp')
        if ss != ' ':
            #ss_hvr = ss[161:182]
            ss_hvr = ss[self.temp_hvr_index_start:self.temp_hvr_index_end_ras4b]
            bad_ss = False

            for i in range(0,len(ss_hvr)):
                if ss_hvr[i]=='T' or ss_hvr[i]=='S' or ss_hvr[i]=='3':
                    ss_hvr_fixed = ss_hvr_fixed + 'C'
                else:
                    ss_hvr_fixed = ss_hvr_fixed + ss_hvr[i]

            for i in range(0,len(ss_hvr_fixed)-1):
                if (ss_hvr_fixed[i]=='C' and ss_hvr_fixed[i+1]=='1') or (ss_hvr_fixed[i]=='2' and ss_hvr_fixed[i+1]=='1'):
                    bad_ss = True

            #print(ss_hvr_fixed)
            # Prevent shrinking of helix below 164
            if ss_hvr_fixed[0] == 'C':
                bad_ss = True

            if bad_ss != True:
                do_hvr_append = True

        # -------------------------------------------------------------
        # execute second martinize
        #LOGGER.debug('executing ('+cmd.format(pdb_crd)+')')
        os.system(cmd.format(pdb_crd))

        # fetch the seconday structure
        ss = FeedbackManager_AA2CG._extract_ss('Protein.itp')
        if ss != ' ':
            #ss_crd = ss[91:108]
            ss_crd = ss[self.temp_crd_index_start_ras4b:self.temp_crd_index_end_ras4b]
            do_crd_append = True

            #if (ss_crd == 'EEEEEEECTTTCSEEEE'):
            #    print("CRD PDB ------------------>>>>>>>")
            #    print(pdb_crd)


        # cleanup the temp
        os.chdir(self.workspace)
        shutil.rmtree(tmp_dir)

        # return value!
        return {'pdb_hvr': pdb_hvr, 'pdb_crd': pdb_crd, 'pdb_full': pdb_full,
                'pdb_hvr_data': pdb_hvr_data, 'pdb_crd_data': pdb_crd_data, 'pdb_full_data': pdb_full_data,
                'ss_hvr_fixed': ss_hvr_fixed, 'ss_crd': ss_crd,
                'do_hvr_append': do_hvr_append, 'do_crd_append': do_crd_append}

    # --------------------------------------------------------------------------

    def check_ss_ras4b(self, coord_ras4b):

        # ----------------------------------------------------------------------
        cmd = 'python martinize_2.6_3_modified.py -f {} -o temp-system.top' \
            ' -x cg.pdb -logfile martinize_log.log -p backbone  -ff martini22 -nt -elastic -dssp mkdssp' \
            ' > /dev/null'


        region1 = 'hvr'
        region3 = 'full'

        ss_hvr_fixed = ''
        do_hvr_append = False

        def _write(_fname, _data):
            with open(_fname, 'w') as fp:
                fp.write(_data)


        p = multiprocessing.current_process()

        # create pdb strings using the coords
        pdb_hvr_data = FeedbackManager_AA2CG.create_pdb(self.main_pdb_hvr_ras4b, coord_ras4b, region1, 'ras4b')
        pdb_full_data = FeedbackManager_AA2CG.create_pdb(self.main_pdb_full_ras4b, coord_ras4b, region3, 'ras4b')

        # ---------------------------------------------------------------------                                         
        # these are the counts before the current process pool started
        _ts = add_timestamp('')[1:]
        pdbid_hvr = '{}-p{}'.format(_ts, p.pid)
        pdbid_full = '{}-p{}'.format(_ts, p.pid)

        pdb_hvr = os.path.basename(Naming.fname_fb_pdb(self.pdbkey_hvr_ras4b[:-4], region1, pdbid_hvr))
        pdb_full = os.path.basename(Naming.fname_fb_pdb(self.pdbkey_full_ras4b[:-9], region3, pdbid_full))

        # ---------------------------------------------------------------------
        # create a temp directory to allow concurrent runs
        tmp_dir = os.path.join(self.workspace, 'check_ss_ras4b_{}'.format(p.pid))
        #LOGGER.debug('Creating temp directory for check_ss_ras4b ({})'.format(tmp_dir))
        mkdir(tmp_dir)
        os.chdir(tmp_dir)

        _write(os.path.join(tmp_dir, pdb_hvr), pdb_hvr_data)
        _write(os.path.join(tmp_dir, pdb_full), pdb_full_data)

        # copy relevant files to the directory
        files_to_copy = ['martinize_2.6_3_modified.py']
        for _ in files_to_copy:
           shutil.copy(os.path.join(self.workspace, _), tmp_dir)

        _write(os.path.join(tmp_dir, pdb_hvr), pdb_hvr_data)
        _write(os.path.join(tmp_dir, pdb_full), pdb_full_data)
        
        # -----------------------------------------------------------------
        # execute first martinize
        #LOGGER.debug('executing ('+cmd.format(pdb_hvr)+')')
        os.system(cmd.format(pdb_hvr))

        # fetch the secondary structure
        ss = FeedbackManager_AA2CG._extract_ss('Protein.itp')
        if ss != ' ':
            ss_hvr = ss[self.temp_hvr_index_start:self.temp_hvr_index_end_ras4b]
            bad_ss = False

            for i in range(0,len(ss_hvr)):
                if ss_hvr[i]=='T' or ss_hvr[i]=='S' or ss_hvr[i]=='3':
                    ss_hvr_fixed = ss_hvr_fixed + 'C'
                else:
                    ss_hvr_fixed = ss_hvr_fixed + ss_hvr[i]

            for i in range(0,len(ss_hvr_fixed)-1):
                if (ss_hvr_fixed[i]=='C' and ss_hvr_fixed[i+1]=='1') or (ss_hvr_fixed[i]=='2' and ss_hvr_fixed[i+1]=='1'):
                    bad_ss = True

            #print(ss_hvr_fixed)                                                                                                           
            # Prevent shrinking of helix below 164                                                                                        
            if ss_hvr_fixed[0] == 'C':
                bad_ss = True

            if bad_ss != True:
                do_hvr_append = True
                

        # cleanup the temp
        os.chdir(self.workspace)
        shutil.rmtree(tmp_dir)

        # return value!
        return {'pdb_hvr': pdb_hvr, 'pdb_full': pdb_full,
                'pdb_hvr_data': pdb_hvr_data, 'pdb_full_data': pdb_full_data,
                'ss_hvr_fixed': ss_hvr_fixed, 'do_hvr_append': do_hvr_append}
                

    def check_ss_ras4araf(self, coord_ras4araf):

        # ----------------------------------------------------------------------
        cmd = 'python martinize_2.6_3_modified.py -f {} -o temp-system.top' \
              ' -x cg.pdb -logfile martinize_log.log -p backbone  -ff martini22 -nt -elastic -dssp mkdssp' \
              ' > /dev/null'


        region1 = 'hvr'
        region2 = 'crd'
        region3 = 'full'

        ss_hvr_fixed = ''
        ss_crd = ''
        do_hvr_append = False
        do_crd_append = False

        def _write(_fname, _data):
            with open(_fname, 'w') as fp:
                fp.write(_data)

        # ----------------------------------------------------------------------
        p = multiprocessing.current_process()

        # create pdb strings using the coords
        pdb_hvr_data = FeedbackManager_AA2CG.create_pdb(self.main_pdb_hvr_ras4araf, coord_ras4araf, region1, 'ras4a')
        pdb_crd_data = FeedbackManager_AA2CG.create_pdb(self.main_pdb_crd_ras4araf, coord_ras4araf, region2, 'ras4a')
        pdb_full_data = FeedbackManager_AA2CG.create_pdb(self.main_pdb_full_ras4araf, coord_ras4araf, region3, 'ras4a')

        # ---------------------------------------------------------------------
        # these are the counts before the current process pool started
        _ts = add_timestamp('')[1:]
        pdbid_hvr = '{}-p{}'.format(_ts, p.pid)
        pdbid_crd = '{}-p{}'.format(_ts, p.pid)
        pdbid_full = '{}-p{}'.format(_ts, p.pid)

        pdb_hvr = os.path.basename(Naming.fname_fb_pdb(self.pdbkey_hvr_ras4araf[:-4], region1, pdbid_hvr))
        pdb_crd = os.path.basename(Naming.fname_fb_pdb(self.pdbkey_crd_ras4araf[:-4], region2, pdbid_crd))
        pdb_full = os.path.basename(Naming.fname_fb_pdb(self.pdbkey_full_ras4araf[:-9], region3, pdbid_full))

        # ---------------------------------------------------------------------
        # create a temp directory to allow concurrent runs
        tmp_dir = os.path.join(self.workspace, 'check_ss_ras4araf_{}'.format(p.pid))
        #LOGGER.debug('Creating temp directory for check_ss_ras4araf ({})'.format(tmp_dir))
        mkdir(tmp_dir)
        os.chdir(tmp_dir)

        _write(os.path.join(tmp_dir, pdb_hvr), pdb_hvr_data)
        _write(os.path.join(tmp_dir, pdb_crd), pdb_crd_data)
        _write(os.path.join(tmp_dir, pdb_full), pdb_full_data)

        # copy relevant files to the directory
        files_to_copy = ['martinize_2.6_3_modified.py']
        for _ in files_to_copy:
           shutil.copy(os.path.join(self.workspace, _), tmp_dir)

        _write(os.path.join(tmp_dir, pdb_hvr), pdb_hvr_data)
        _write(os.path.join(tmp_dir, pdb_crd), pdb_crd_data)
        _write(os.path.join(tmp_dir, pdb_full), pdb_full_data)

        # -----------------------------------------------------------------
        # execute first martinize
        #LOGGER.debug('executing ('+cmd.format(pdb_hvr)+')')
        os.system(cmd.format(pdb_hvr))

        # fetch the secondary structure
        ss = FeedbackManager_AA2CG._extract_ss('Protein.itp')
        if ss != ' ':
            #ss_hvr = ss[161:182]
            ss_hvr = ss[self.temp_hvr_index_start:self.temp_hvr_index_end_ras4a]
            bad_ss = False

            for i in range(0,len(ss_hvr)):
                if ss_hvr[i]=='T' or ss_hvr[i]=='S' or ss_hvr[i]=='3':
                    ss_hvr_fixed = ss_hvr_fixed + 'C'
                else:
                    ss_hvr_fixed = ss_hvr_fixed + ss_hvr[i]

            for i in range(0,len(ss_hvr_fixed)-1):
                if (ss_hvr_fixed[i]=='C' and ss_hvr_fixed[i+1]=='1') or (ss_hvr_fixed[i]=='2' and ss_hvr_fixed[i+1]=='1'):
                    bad_ss = True

            #print(ss_hvr_fixed)
            # Prevent shrinking of helix below 164
            if ss_hvr_fixed[0] == 'C':
                bad_ss = True

            if bad_ss != True:
                do_hvr_append = True

        # -------------------------------------------------------------
        # execute second martinize
        #LOGGER.debug('executing ('+cmd.format(pdb_crd)+')')
        os.system(cmd.format(pdb_crd))

        # fetch the seconday structure
        ss = FeedbackManager_AA2CG._extract_ss('Protein.itp')
        if ss != ' ':
            #ss_crd = ss[91:108]
            ss_crd = ss[self.temp_crd_index_start_ras4a:self.temp_crd_index_end_ras4a]
            do_crd_append = True

            #if (ss_crd == 'EEEEEEECTTTCSEEEE'):
            #    print("CRD PDB ------------------>>>>>>>")
            #    print(pdb_crd)


        # cleanup the temp
        os.chdir(self.workspace)
        shutil.rmtree(tmp_dir)

        # return value!
        return {'pdb_hvr': pdb_hvr, 'pdb_crd': pdb_crd, 'pdb_full': pdb_full,
                'pdb_hvr_data': pdb_hvr_data, 'pdb_crd_data': pdb_crd_data, 'pdb_full_data': pdb_full_data,
                'ss_hvr_fixed': ss_hvr_fixed, 'ss_crd': ss_crd,
                'do_hvr_append': do_hvr_append, 'do_crd_append': do_crd_append}

    # --------------------------------------------------------------------------

    def check_ss_ras4a(self, coord_ras4a):

        # ----------------------------------------------------------------------
        cmd = 'python martinize_2.6_3_modified.py -f {} -o temp-system.top' \
            ' -x cg.pdb -logfile martinize_log.log -p backbone  -ff martini22 -nt -elastic -dssp mkdssp' \
            ' > /dev/null'


        region1 = 'hvr'
        region3 = 'full'

        ss_hvr_fixed = ''
        do_hvr_append = False

        def _write(_fname, _data):
            with open(_fname, 'w') as fp:
                fp.write(_data)


        p = multiprocessing.current_process()

        # create pdb strings using the coords
        pdb_hvr_data = FeedbackManager_AA2CG.create_pdb(self.main_pdb_hvr_ras4a, coord_ras4a, region1, 'ras4a')
        pdb_full_data = FeedbackManager_AA2CG.create_pdb(self.main_pdb_full_ras4a, coord_ras4a, region3, 'ras4a')

        # ---------------------------------------------------------------------                                         
        # these are the counts before the current process pool started
        _ts = add_timestamp('')[1:]
        pdbid_hvr = '{}-p{}'.format(_ts, p.pid)
        pdbid_full = '{}-p{}'.format(_ts, p.pid)

        pdb_hvr = os.path.basename(Naming.fname_fb_pdb(self.pdbkey_hvr_ras4a[:-4], region1, pdbid_hvr))
        pdb_full = os.path.basename(Naming.fname_fb_pdb(self.pdbkey_full_ras4a[:-9], region3, pdbid_full))

        # ---------------------------------------------------------------------
        # create a temp directory to allow concurrent runs
        tmp_dir = os.path.join(self.workspace, 'check_ss_ras4a_{}'.format(p.pid))
        #LOGGER.debug('Creating temp directory for check_ss_ras4a ({})'.format(tmp_dir))
        mkdir(tmp_dir)
        os.chdir(tmp_dir)

        _write(os.path.join(tmp_dir, pdb_hvr), pdb_hvr_data)
        _write(os.path.join(tmp_dir, pdb_full), pdb_full_data)

        # copy relevant files to the directory
        files_to_copy = ['martinize_2.6_3_modified.py']
        for _ in files_to_copy:
           shutil.copy(os.path.join(self.workspace, _), tmp_dir)

        _write(os.path.join(tmp_dir, pdb_hvr), pdb_hvr_data)
        _write(os.path.join(tmp_dir, pdb_full), pdb_full_data)
        
        # -----------------------------------------------------------------
        # execute first martinize
        #LOGGER.debug('executing ('+cmd.format(pdb_hvr)+')')
        os.system(cmd.format(pdb_hvr))

        # fetch the secondary structure
        ss = FeedbackManager_AA2CG._extract_ss('Protein.itp')
        if ss != ' ':
            ss_hvr = ss[self.temp_hvr_index_start:self.temp_hvr_index_end_ras4a]
            bad_ss = False

            for i in range(0,len(ss_hvr)):
                if ss_hvr[i]=='T' or ss_hvr[i]=='S' or ss_hvr[i]=='3':
                    ss_hvr_fixed = ss_hvr_fixed + 'C'
                else:
                    ss_hvr_fixed = ss_hvr_fixed + ss_hvr[i]

            for i in range(0,len(ss_hvr_fixed)-1):
                if (ss_hvr_fixed[i]=='C' and ss_hvr_fixed[i+1]=='1') or (ss_hvr_fixed[i]=='2' and ss_hvr_fixed[i+1]=='1'):
                    bad_ss = True

            #print(ss_hvr_fixed)                                                                                                           
            # Prevent shrinking of helix below 164                                                                                        
            if ss_hvr_fixed[0] == 'C':
                bad_ss = True

            if bad_ss != True:
                do_hvr_append = True
                

        # cleanup the temp
        os.chdir(self.workspace)
        shutil.rmtree(tmp_dir)

        # return value!
        return {'pdb_hvr': pdb_hvr, 'pdb_full': pdb_full,
                'pdb_hvr_data': pdb_hvr_data, 'pdb_full_data': pdb_full_data,
                'ss_hvr_fixed': ss_hvr_fixed, 'do_hvr_append': do_hvr_append}


    # report the analysis
    def report(self):

        # ----------------------------------------------------------------------
        s = '({}): naggregates = {}, naggs_at_last_report = {}'
        s = s.format(self.name, self.naggregates, self.naggs_at_last_report)
        LOGGER.info('> entering report function: {}'.format(s))
        LOGGER.info('type of agg: {} {}'.format(type(self.naggregates), type(self.naggs_at_last_report)))

        if self.naggregates == self.naggs_at_last_report:
            LOGGER.info('Skip reporting stale analysis'.format(s))
            return

        # ----------------------------------------------------------------------
        key = Naming.fb_aa_key_ras4a('manager', nreports=self.nreports)
        LOGGER.info('Reporting analysis {} fbpath = ({}), key ({})'.format(s, self.outpath, key))

        # should be iointerface here for final output
        self.iointerface.save_npz(self.outpath, key,
                                  {'aggregated_rmsf_ras4braf': self.aggregate_rmsf_ras4braf,
                                   'aggregated_coord_ras4braf': self.aggregate_coord_ras4braf,
                                   'aggregated_rmsf_ras4b': self.aggregate_rmsf_ras4b,
                                   'aggregated_coord_ras4b': self.aggregate_coord_ras4b,
                                   'aggregated_rmsf_ras4araf': self.aggregate_rmsf_ras4araf,
                                   'aggregated_coord_ras4araf': self.aggregate_coord_ras4araf,
                                   'aggregated_rmsf_ras4a': self.aggregate_rmsf_ras4a,
                                   'aggregated_coord_ras4a': self.aggregate_coord_ras4a,
                                   'num_aggs': self.naggregates,
                                   'num_reports': self.nreports},
                                   interfaces=['simple'])

        # ----------------------------------------------------------------------
        self.naggs_at_last_report = self.naggregates
        self.nreports += 1

        # ----------------------------------------------------------------------
        s = 'Reported analysis ({}): #reports = {}; #aggregates = {}'
        s = s.format(self.name, self.nreports, self.naggregates)
        LOGGER.info(s)

        if self.tot_num_hvr_ras4braf >= self.nframes_fb_threshold_ras4braf and self.tot_num_hvr_ras4araf >= self.nframes_fb_threshold_ras4araf and self.tot_num_hvr_ras4b >= self.nframes_fb_threshold_ras4b and self.tot_num_hvr_ras4a >= self.nframes_fb_threshold_ras4a:            
            
            self.nframes_fb_threshold_ras4braf = self.tot_num_hvr_ras4braf + self.nframes_fb_increment
            self.nframes_fb_threshold_ras4araf = self.tot_num_hvr_ras4araf + self.nframes_fb_increment
            self.nframes_fb_threshold_ras4b = self.tot_num_hvr_ras4b + self.nframes_fb_increment
            self.nframes_fb_threshold_ras4a = self.tot_num_hvr_ras4a + self.nframes_fb_increment

            self.ras4braf_feedback()
            self.ras4b_feedback()
            self.ras4araf_feedback()
            self.ras4a_feedback()



    # ---------------------------------------------------------------------------
    def ras4braf_feedback(self):

        LOGGER.info('Doing Feedback for RAS4B-RAF')
        crd_updated = False
        hvr_updated = False

        # To update the itp always based on the initial structure (no history of previous updates)
        self.main_itp_ras4braf = os.path.join(self.respath, 'RAS_RAF_TERNARY_scfix-CYFpos.itp')
 
        
        # ----------------------------------------------------------------------
        #if self.tot_num_hvr > self.nframes_fb and self.perc_hvr > 0.0 and self.most_hvr != self.prev_hvr:
        if self.perc_hvr_ras4braf > self.hvr_th_ras4braf and self.most_hvr_ras4braf != self.prev_hvr_ras4braf:

            self.ras4braf_feedback_done = True
            
            self.pdb_most_hvr_ss_ras4braf = self.dict_hvr_ss_pdb_ras4braf.get(self.most_hvr_ras4braf)
            
            out_itp_1_ras4braf = 'Protein_hvr1_ras4braf'
            self._martinize_with_dssp(self.pdb_hvr_tar_ras4braf, self.pdb_most_hvr_ss_ras4braf, out_itp_1_ras4braf)

            out_itp_1_ras4braf = out_itp_1_ras4braf + '.itp'
            with open(out_itp_1_ras4braf) as fp:
                for line in fp:
                    data = line.split()
                    if len(data) > 1 and data[1] == "Secondary":
                        data_ss = next(fp).split()
                        ss = data_ss[1]
                        break



            ss_revised = ss[0:self.temp_hvr_index_start] + self.most_hvr_ras4braf + ss[self.temp_hvr_index_end_ras4b:-1] + ss[-1]
            


            for i in range(0,len(ss_revised)-1):
                resid = i+3
                if (ss_revised[i] == '2' and ss_revised[i+1] == 'C' and resid > 160 and resid < 185):
                    self.hvr_end_main_ras4braf = resid



            ss_hvr_revised_fixed = ''

            for i in range(0,len(ss_revised)):
                if ss_revised[i] in ['1', '2']:
                    ss_hvr_revised_fixed = ss_hvr_revised_fixed + 'H'
                else:
                    ss_hvr_revised_fixed = ss_hvr_revised_fixed + ss_revised[i]

            # use new AA structure if helix grows beyond residue 172, otherwise use the initial structure
            print(self.hvr_end_main_ras4braf)
            if self.hvr_end_main_ras4braf > 172:
                out_itp_2_ras4braf = 'Protein_hvr2_ras4braf'
                self._martinize_with_ss(self.pdb_hvr_tar_ras4braf, self.pdb_most_hvr_ss_ras4braf, ss_hvr_revised_fixed, out_itp_2_ras4braf)
            else:
                out_itp_2_ras4braf = 'Protein_hvr2_ras4braf'
                self._martinize_with_ss(self.pdb_hvr_tar_ras4braf, self.pdbkey_hvr_ras4braf, ss_hvr_revised_fixed, out_itp_2_ras4braf)

            out_itp_2_ras4braf = out_itp_2_ras4braf + '.itp'
            self.revised_itp = out_itp_2_ras4braf
            self.main_itp_ras4braf = self.update_itp_hvr()
            self.prev_hvr_ras4braf = self.most_hvr_ras4braf
            self.pdb_full_hvr_ras4braf = self.dict_full_ss_pdb_ras4braf.get(self.most_hvr_ras4braf)
            hvr_updated = True

            self.ras4braf_feedback_done = False

        else:
            LOGGER.info('Threshold for HVR is not exceeded for RAS4B-RAF, continue with current itp file')
            hvr_updated = False

        if self.perc_crd_ras4braf > self.crd_th_ras4braf and self.most_crd_ras4braf != self.prev_crd_ras4braf:

            self.ras4braf_feedback_done = True

            self.pdb_most_crd_ss_ras4braf = self.dict_crd_ss_pdb_ras4braf.get(self.most_crd_ras4braf)
            
            LOGGER.debug("most_crd=%s\npdb_crd_tar=%s\npdb_most_crd_ss=%s", self.most_crd_ras4braf, self.pdb_crd_tar_ras4braf, self.pdb_most_crd_ss_ras4braf)
            out_itp_3_ras4braf = 'Protein_crd1_ras4braf'
            self._martinize_with_dssp(self.pdb_crd_tar_ras4braf, self.pdb_most_crd_ss_ras4braf, out_itp_3_ras4braf)
            out_itp_3_ras4braf = out_itp_3_ras4braf + '.itp'
            self.revised_itp = out_itp_3_ras4braf
            self.main_itp_ras4braf = self.update_itp_crd()
            self.prev_crd_ras4braf = self.most_crd_ras4braf
            self.pdb_full_crd_ras4braf = self.dict_full_ss_pdb_ras4braf.get(self.most_crd_ras4braf)
            crd_updated = True

            self.ras4braf_feedback_done = False

        else:
            LOGGER.info('Threshold for CRD is not exceeded for RAS4B-RAF, continue with current itp file')
            crd_updated = False


        if not hvr_updated and not crd_updated:
            return

        if hvr_updated == True and crd_updated == False and self.pdb_most_crd_ss_ras4braf != ' ':
            self.ras4braf_feedback_done = True
            out_itp_4_ras4braf = 'Protein_crd2_ras4braf'
            self._martinize_with_dssp(self.pdb_crd_tar_ras4braf, self.pdb_most_crd_ss_ras4braf, out_itp_4_ras4braf)
            out_itp_4_ras4braf = out_itp_4_ras4braf + '.itp'
            self.revised_itp = out_itp_4_ras4braf
            self.main_itp_ras4braf = self.update_itp_crd()
            self.ras4braf_feedback_done = False

        elif hvr_updated == False and crd_updated == True and self.pdb_most_hvr_ss_ras4braf != ' ':
            self.ras4braf_feedback_done = True
            out_itp_6_ras4braf = 'Protein_hvr4_ras4braf'
            self._martinize_with_dssp(self.pdb_hvr_tar_ras4braf, self.pdb_most_hvr_ss_ras4braf, out_itp_6_ras4braf)

            out_itp_6_ras4braf = out_itp_6_ras4braf + '.itp'
            with open(out_itp_6_ras4braf) as fp:
                for line in fp:
                    data = line.split()
                    if len(data) > 1 and data[1] == "Secondary":
                        data_ss = next(fp).split()
                        ss = data_ss[1]
                        break


            ss_revised = ss[0:self.temp_hvr_index_start] + self.most_hvr_ras4braf + ss[self.temp_hvr_index_end_ras4b:-1] + ss[-1]


            for i in range(0,len(ss_revised)-1):
                resid = i+3
                if (ss_revised[i] == '2' and ss_revised[i+1] == 'C' and resid > 160 and resid < 185):
                    self.hvr_end_main_ras4braf = resid



            ss_hvr_revised_fixed = ''

            for i in range(0,len(ss_revised)):                                                               
                if ss_revised[i] in ['1', '2']:
                    ss_hvr_revised_fixed = ss_hvr_revised_fixed + 'H'
                else:
                    ss_hvr_revised_fixed = ss_hvr_revised_fixed + ss_revised[i]

            
            if self.hvr_end_main_ras4braf > 172:
                out_itp_5_ras4braf = 'Protein_hvr3_ras4braf'
                self._martinize_with_ss(self.pdb_hvr_tar_ras4braf, self.pdb_most_hvr_ss_ras4braf, ss_hvr_revised_fixed, out_itp_5_ras4braf)
            else:
                out_itp_5 = 'Protein_hvr3'
                self._martinize_with_ss(self.pdb_hvr_tar_ras4braf, self.pdbkey_hvr_ras4braf, ss_hvr_revised_fixed, out_itp_5_ras4braf)

            out_itp_5_ras4braf = out_itp_5_ras4braf + '.itp'
            self.revised_itp = out_itp_5_ras4braf
            self.main_itp_ras4braf = self.update_itp_hvr()
            self.ras4braf_feedback_done = False

        if hvr_updated == True or crd_updated == True:

            ts = add_timestamp('')
            main_itp_final_ras4braf = self.main_itp_ras4braf[:-13] + ts + '.itp'

            fullFile = self.main_itp_ras4braf[:-13] + ts + '-PULL' + '.itp'
            if not os.path.isfile(fullFile):
                # PULL file is missing creat
                with open(self.main_itp_ras4braf, "r") as sources:
                    lines = sources.readlines()
                    with open(fullFile, "w") as sources:
                        for line in lines:
                            line = re.sub(r'RAS_RAF   1', 'RASR_PULL   1', line)
                            line = re.sub(r'  414   415      1', ';  414   415      1', line)
                            line = re.sub(r'  413   414   415      2', ';  413   414   415      2', line)
                            line = re.sub(r'  414   415   416      2', ';  414   415   416      2', line)
                            sources.write(line)
                LOGGER.warning("Feedback AA to CG - PULL file for .itp {} not found - create file {}".format(main_itp_final_ras4braf, fullFile))


        # ----------------------------------------------------------------------
        # write the file
        with open(self.main_itp_ras4braf) as fp:
            for line in fp:
                data = line.split()
                if len(data)>1 and data[1] == "Secondary":
                    data_ss = next(fp).split()
                    ss_full = data_ss[1]
                    break


        if hvr_updated:
            filename_hvr = self._create_tmp(self.pdb_full_tar_ras4braf, self.pdb_full_hvr_ras4braf)
            with open(filename_hvr,'r+') as fp:
                lines = fp.readlines()
                lines.insert(0, "REMARK " + "300 " + ss_full + "\n")
                fp.seek(0)
                fp.writelines(lines)

        if crd_updated:
            filename_crd = self._create_tmp(self.pdb_full_tar_ras4braf, self.pdb_full_crd_ras4braf)
            if hvr_updated and filename_hvr != filename_crd:
                with open(filename_crd,'r+') as fp:
                    lines = fp.readlines()
                    lines.insert(0, "REMARK " + "300 " + ss_full + "\n")
                    fp.seek(0)
                    fp.writelines(lines)
            elif not hvr_updated:
                with open(filename_crd,'r+') as fp:
                    lines = fp.readlines()
                    lines.insert(0, "REMARK " + "300 " + ss_full + "\n")
                    fp.seek(0)
                    fp.writelines(lines)


        if hvr_updated and crd_updated:
            pdb_full_hvr_ts = 'RAS_RAF_TERNARY_scfix-CYFpos' + '-HVR' + ts + '.pdb'
            pdb_full_crd_ts = 'RAS_RAF_TERNARY_scfix-CYFpos' + '-CRD' + ts + '.pdb'
            shutil.copy(self.main_itp_ras4braf, main_itp_final_ras4braf)
            shutil.copy(filename_crd, pdb_full_crd_ts)
            shutil.copy(filename_hvr, pdb_full_hvr_ts)
        elif not hvr_updated and crd_updated:
            pdb_full_crd_ts = 'RAS_RAF_TERNARY_scfix-CYFpos' + '-CRD' + ts + '.pdb'
            shutil.copy(self.main_itp_ras4braf, main_itp_final_ras4braf)
            shutil.copy(filename_crd, pdb_full_crd_ts)
        elif hvr_updated and not crd_updated:
            pdb_full_hvr_ts = 'RAS_RAF_TERNARY_scfix-CYFpos' + '-HVR' + ts + '.pdb'
            shutil.copy(self.main_itp_ras4braf, main_itp_final_ras4braf)
            shutil.copy(filename_hvr, pdb_full_hvr_ts)
        else:
            LOGGER.info('There is no update in both HVR and CRD regions')


        if hvr_updated:
            os.remove(filename_hvr)

        if crd_updated:
            if crd_updated and hvr_updated and filename_hvr != filename_crd:
                os.remove(filename_crd)
            elif not hvr_updated:
                os.remove(filename_crd)

        # ----------------------------------------------------------------------
        ######### Helgi's protocol ##########
        # ----------------------------------------------------------------------
        if hvr_updated or crd_updated:

            os.chdir(self.workspace)

            os.mkdir('ddcMDParams' + ts)
            os.chdir('ddcMDParams'+ ts)

            files_copy = Naming.dir_res('ddcMD')+'/martini2objfiles/*'
            for f in glob.glob(files_copy):
                shutil.copy(f,os.getcwd())

            #TODO: Fikret, please check if this needs a so/q;p
            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + main_itp_final_ras4braf + " > temp_RAS4b_RAF.itp")

            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + self.main_itp_latest_ras4b + " > temp_RAS4b.itp")
            
            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + self.main_itp_latest_ras4araf + " > temp_RAS4a_RAF.itp")
            
            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + self.main_itp_latest_ras4a + " > temp_RAS4a.itp")
            

            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4b_RAF.itp > RAS_RAF.itp')

            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4b.itp > RAS.itp')
            
            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4a_RAF.itp > RAS4a_RAF.itp')
            
            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4a.itp > RAS4a.itp')

            os.system('rm temp_RAS4b_RAF.itp temp_RAS4b.itp temp_RAS4a.itp temp_RAS4a_RAF.itp')

            os.system('martini2obj -i itpList -t proItpList -p  martini_v2.1-dna.itp -o martini.data -l molecule.data > itp_to_ddcMD.out')

            os.system("sed -i -n '/moleculeClass MOLECULECLASS/,$p' molecule.data")

        # ----------------------------------------------------------------------
        self.main_itp_latest_ras4braf = main_itp_final_ras4braf
    
        # ras-raf
        crd_updated = False
        hvr_updated = False

    # --------------------------------------------------------------------------
    
    def ras4b_feedback(self):

        LOGGER.info('Doing Feedback for RAS4B')
        hvr_updated = False

        # To update the itp always based on the initial structure (no history of previous updates)
        self.main_itp_ras4b = os.path.join(self.respath,'RAS_scfix-CYFpos.itp')

        
        # ----------------------------------------------------------------------
        if self.perc_hvr_ras4b > self.hvr_th_ras4b and self.most_hvr_ras4b != self.prev_hvr_ras4b:

            self.ras4b_feedback_done = True

            self.pdb_most_hvr_ss_ras4b = self.dict_hvr_ss_pdb_ras4b.get(self.most_hvr_ras4b)
            
            out_itp_1_ras4b = 'Protein_hvr1_ras4b'
            self._martinize_with_dssp(self.pdb_hvr_tar_ras4b, self.pdb_most_hvr_ss_ras4b, out_itp_1_ras4b)

            out_itp_1_ras4b = out_itp_1_ras4b + '.itp'
            with open(out_itp_1_ras4b) as fp:
                for line in fp:
                    data = line.split()
                    if len(data) > 1 and data[1] == "Secondary":
                        data_ss = next(fp).split()
                        ss = data_ss[1]
                        break


            ss_revised = ss[0:self.temp_hvr_index_start] + self.most_hvr_ras4b + ss[self.temp_hvr_index_end_ras4b:-1] + ss[-1]

            #ss_revised = ss[0:self.temp_hvr_index_start] + 'HHHHHH2222CCCCCCCCCCC' + ss[self.temp_hvr_index_end_ras4b:-1] + ss[-1]


            for i in range(0,len(ss_revised)-1):
                resid = i+3
                if (ss_revised[i] == '2' and ss_revised[i+1] == 'C' and resid > 160 and resid < 185):
                    self.hvr_end_main_ras4b = resid


            ss_hvr_revised_fixed = ''

            for i in range(0,len(ss_revised)):
                if ss_revised[i] in ['1', '2']:
                    ss_hvr_revised_fixed = ss_hvr_revised_fixed + 'H'
                else:
                    ss_hvr_revised_fixed = ss_hvr_revised_fixed + ss_revised[i]

            # use new AA structure if helix grows beyond residue 172, otherwise use the initial structure
            print(self.hvr_end_main_ras4b)
            if self.hvr_end_main_ras4b > 172:
                out_itp_2_ras4b = 'Protein_hvr2_ras4b'
                self._martinize_with_ss(self.pdb_hvr_tar_ras4b, self.pdb_most_hvr_ss_ras4b, ss_hvr_revised_fixed, out_itp_2_ras4b)
            else:
                out_itp_2_ras4b = 'Protein_hvr2_ras4b'
                self._martinize_with_ss(self.pdb_hvr_tar_ras4b, self.pdbkey_hvr_ras4b, ss_hvr_revised_fixed, out_itp_2_ras4b)

            out_itp_2_ras4b = out_itp_2_ras4b + '.itp'
            self.revised_itp = out_itp_2_ras4b
            self.main_itp_ras4b = self.update_itp_hvr()
            self.prev_hvr_ras4b = self.most_hvr_ras4b
            self.pdb_full_hvr_ras4b = self.dict_full_ss_pdb_ras4b.get(self.most_hvr_ras4b)
            hvr_updated = True

            self.ras4b_feedback_done = False
            
        else:
            LOGGER.info('Threshold for HVR is not exceeded for RAS4B, continue with current itp file')
            hvr_updated = False



        if not hvr_updated:
            return

 

        if hvr_updated:

            ts = add_timestamp('')

            main_itp_final_ras4b = self.main_itp_ras4b[:-13] + ts + '.itp'

            fullFile = self.main_itp_ras4b[:-13] + ts + '-PULL' + '.itp'
            if not os.path.isfile(fullFile):
                # PULL file is missing creat                                                                                                  
                with open(self.main_itp_ras4b, "r") as sources:
                    lines = sources.readlines()
                    with open(fullFile, "w") as sources:
                        for line in lines:
                            line = re.sub(r'RAS   1', 'RAS_PULL   1', line)
                            line = re.sub(r'  414   415      1', ';  414   415      1', line)
                            line = re.sub(r'  413   414   415      2', ';  413   414   415      2', line)
                            line = re.sub(r'  414   415   416      2', ';  414   415   416      2', line)
                            sources.write(line)
                LOGGER.warning("Feedback AA to CG - PULL file for .itp {} not found - create file {}".format(main_itp_final_ras4b, fullFile))



        # ----------------------------------------------------------------------
        # write the file
        with open(self.main_itp_ras4b) as fp:
            for line in fp:
                data = line.split()
                if len(data)>1 and data[1] == "Secondary":
                    data_ss = next(fp).split()
                    ss_full = data_ss[1]
                    break


        if hvr_updated:
            filename_hvr = self._create_tmp(self.pdb_full_tar_ras4b, self.pdb_full_hvr_ras4b)
            with open(filename_hvr,'r+') as fp:
                lines = fp.readlines()
                lines.insert(0, "REMARK " + "300 " + ss_full + "\n")
                fp.seek(0)
                fp.writelines(lines)



        if hvr_updated:
            pdb_full_hvr_ts = 'RAS_scfix-CYFpos' + '-HVR' + ts + '.pdb'
            shutil.copy(self.main_itp_ras4b, main_itp_final_ras4b)
            shutil.copy(filename_hvr, pdb_full_hvr_ts)
        else:
            LOGGER.info('There is no update in HVR region for RAS4B')


        if hvr_updated:
            os.remove(filename_hvr)


        # ----------------------------------------------------------------------
        ######### Helgi's protocol ##########
        # ----------------------------------------------------------------------
        if hvr_updated:

            os.chdir(self.workspace)

            os.mkdir('ddcMDParams' + ts)
            os.chdir('ddcMDParams'+ ts)

            files_copy = Naming.dir_res('ddcMD')+'/martini2objfiles/*'
            for f in glob.glob(files_copy):
                shutil.copy(f,os.getcwd())


            #TODO: Fikret, please check if this needs a so/q;p
            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + self.main_itp_latest_ras4braf + " > temp_RAS4b_RAF.itp")
            
            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + main_itp_final_ras4b + " > temp_RAS4b.itp")
            
            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + self.main_itp_latest_ras4araf + " > temp_RAS4a_RAF.itp")
            
            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + self.main_itp_latest_ras4a + " > temp_RAS4a.itp")

            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4b_RAF.itp > RAS_RAF.itp')
            
            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4b.itp > RAS.itp')
            
            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4a_RAF.itp > RAS4a_RAF.itp')
            
            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4a.itp > RAS4a.itp')


            os.system('rm temp_RAS4b_RAF.itp temp_RAS4b.itp temp_RAS4a.itp temp_RAS4a_RAF.itp')

            os.system('martini2obj -i itpList -t proItpList -p  martini_v2.1-dna.itp -o martini.data -l molecule.data > itp_to_ddcMD.out')

            os.system("sed -i -n '/moleculeClass MOLECULECLASS/,$p' molecule.data")

        # ----------------------------------------------------------------------
        
        self.main_itp_latest_ras4b = main_itp_final_ras4b
        
        # ras4b
        hvr_updated = False

    # --------------------------------------------------------------------------
        
    def ras4araf_feedback(self):

        LOGGER.info('Doing Feedback for RAS4A-RAF')
        crd_updated = False
        hvr_updated = False

        # To update the itp always based on the initial structure (no history of previous updates)
        self.main_itp_ras4araf = os.path.join(self.respath, 'RAS4A_RAF_scfix-CYFpos.itp')
 
        
        # ----------------------------------------------------------------------
        if self.perc_hvr_ras4araf > self.hvr_th_ras4araf and self.most_hvr_ras4araf != self.prev_hvr_ras4araf:

            self.ras4araf_feedback_done = True
            
            self.pdb_most_hvr_ss_ras4araf = self.dict_hvr_ss_pdb_ras4araf.get(self.most_hvr_ras4araf)
            
            out_itp_1_ras4araf = 'Protein_hvr1_ras4araf'
            self._martinize_with_dssp(self.pdb_hvr_tar_ras4araf, self.pdb_most_hvr_ss_ras4araf, out_itp_1_ras4araf)

            out_itp_1_ras4araf = out_itp_1_ras4araf + '.itp'
            with open(out_itp_1_ras4araf) as fp:
                for line in fp:
                    data = line.split()
                    if len(data) > 1 and data[1] == "Secondary":
                        data_ss = next(fp).split()
                        ss = data_ss[1]
                        break


            ss_revised = ss[0:self.temp_hvr_index_start] + self.most_hvr_ras4araf + ss[self.temp_hvr_index_end_ras4a:-1] + ss[-1]
            

            
            for i in range(0,len(ss_revised)-1):
                resid = i+3
                if (ss_revised[i] == '2' and ss_revised[i+1] == 'C' and resid > 160 and resid < 185):
                    self.hvr_end_main_ras4araf = resid



            ss_hvr_revised_fixed = ''

            for i in range(0,len(ss_revised)):
                #if ss_revised[i] == '1' or ss_revised[i] == '2':
                if ss_revised[i] in ['1', '2']:
                    ss_hvr_revised_fixed = ss_hvr_revised_fixed + 'H'
                else:
                    ss_hvr_revised_fixed = ss_hvr_revised_fixed + ss_revised[i]

            
            # use new AA structure if helix grows beyond residue 172, otherwise use the initial structure
            print(self.hvr_end_main_ras4araf)
            if self.hvr_end_main_ras4araf > 172:
                out_itp_2_ras4araf = 'Protein_hvr2_ras4araf'
                self._martinize_with_ss(self.pdb_hvr_tar_ras4a_elastic, self.pdbkey_hvr_ras4a_elastic, ss_hvr_revised_fixed, out_itp_2_ras4araf)
            else:
                out_itp_2_ras4araf = 'Protein_hvr2_ras4araf'
                self._martinize_with_ss(self.pdb_hvr_tar_ras4a_elastic, self.pdbkey_hvr_ras4a_elastic, ss_hvr_revised_fixed, out_itp_2_ras4araf)

            out_itp_2_ras4araf = out_itp_2_ras4araf + '.itp'
            self.revised_itp = out_itp_2_ras4araf
            self.main_itp_ras4araf = self.update_itp_hvr()
            self.prev_hvr_ras4araf = self.most_hvr_ras4araf
            self.pdb_full_hvr_ras4araf = self.dict_full_ss_pdb_ras4araf.get(self.most_hvr_ras4araf)
            hvr_updated = True

            self.ras4araf_feedback_done = False

        else:
            LOGGER.info('Threshold for HVR is not exceeded for RAS4A-RAF, continue with current itp file')
            hvr_updated = False


        if self.perc_crd_ras4araf > self.crd_th_ras4araf and self.most_crd_ras4araf != self.prev_crd_ras4araf:

            self.ras4araf_feedback_done = True

            self.pdb_most_crd_ss_ras4araf = self.dict_crd_ss_pdb_ras4araf.get(self.most_crd_ras4araf)
            
            LOGGER.debug("most_crd=%s\npdb_crd_tar=%s\npdb_most_crd_ss=%s", self.most_crd_ras4araf, self.pdb_crd_tar_ras4araf, self.pdb_most_crd_ss_ras4araf)
            out_itp_3_ras4araf = 'Protein_crd1_ras4araf'
            self._martinize_with_dssp(self.pdb_crd_tar_ras4araf, self.pdb_most_crd_ss_ras4araf, out_itp_3_ras4araf)
            out_itp_3_ras4araf = out_itp_3_ras4araf + '.itp'
            self.revised_itp = out_itp_3_ras4araf
            self.main_itp_ras4araf = self.update_itp_crd()
            self.prev_crd_ras4araf = self.most_crd_ras4araf
            self.pdb_full_crd_ras4araf = self.dict_full_ss_pdb_ras4araf.get(self.most_crd_ras4araf)
            crd_updated = True

            self.ras4araf_feedback_done = False

        else:
            LOGGER.info('Threshold for CRD is not exceeded for RAS4A-RAF, continue with current itp file')
            crd_updated = False


        if not hvr_updated and not crd_updated:
            return

        if hvr_updated == True and crd_updated == False and self.pdb_most_crd_ss_ras4araf != ' ':
            self.ras4araf_feedback_done = True
            out_itp_4_ras4araf = 'Protein_crd2_ras4araf'
            self._martinize_with_dssp(self.pdb_crd_tar_ras4araf, self.pdb_most_crd_ss_ras4araf, out_itp_4_ras4araf)
            out_itp_4_ras4araf = out_itp_4_ras4araf + '.itp'
            self.revised_itp = out_itp_4_ras4araf
            self.main_itp_ras4araf = self.update_itp_crd()
            self.ras4araf_feedback_done = False

        elif hvr_updated == False and crd_updated == True and self.pdb_most_hvr_ss_ras4araf != ' ':
            self.ras4araf_feedback_done = True
            out_itp_6_ras4araf = 'Protein_hvr4_ras4araf'
            self._martinize_with_dssp(self.pdb_hvr_tar_ras4araf, self.pdb_most_hvr_ss_ras4araf, out_itp_6_ras4araf)

            out_itp_6_ras4araf = out_itp_6_ras4araf + '.itp'
            with open(out_itp_6_ras4araf) as fp:
                for line in fp:
                    data = line.split()
                    if len(data) > 1 and data[1] == "Secondary":
                        data_ss = next(fp).split()
                        ss = data_ss[1]
                        break
                                                                       

            ss_revised = ss[0:self.temp_hvr_index_start] + self.most_hvr_ras4araf + ss[self.temp_hvr_index_end_ras4a:-1] + ss[-1]


            for i in range(0,len(ss_revised)-1):
                resid = i+3
                if (ss_revised[i] == '2' and ss_revised[i+1] == 'C' and resid > 160 and resid < 185):
                    self.hvr_end_main_ras4araf = resid



            ss_hvr_revised_fixed = ''

            for i in range(0,len(ss_revised)):                                                             
                if ss_revised[i] in ['1', '2']:
                    ss_hvr_revised_fixed = ss_hvr_revised_fixed + 'H'
                else:
                    ss_hvr_revised_fixed = ss_hvr_revised_fixed + ss_revised[i]

            
            if self.hvr_end_main_ras4araf > 172:
                out_itp_5_ras4araf = 'Protein_hvr3_ras4araf'
                self._martinize_with_ss(self.pdb_hvr_tar_ras4a_elastic, self.pdbkey_hvr_ras4a_elastic, ss_hvr_revised_fixed, out_itp_5_ras4araf)
            else:
                out_itp_5_ras4araf = 'Protein_hvr3_ras4araf'
                self._martinize_with_ss(self.pdb_hvr_tar_ras4a_elastic, self.pdbkey_hvr_ras4a_elastic, ss_hvr_revised_fixed, out_itp_5_ras4araf)

            out_itp_5_ras4araf = out_itp_5_ras4araf + '.itp'
            self.revised_itp = out_itp_5_ras4araf
            self.main_itp_ras4araf = self.update_itp_hvr()
            self.ras4araf_feedback_done = False

        if hvr_updated == True or crd_updated == True:

            ts = add_timestamp('')
            main_itp_final_ras4araf = self.main_itp_ras4araf[:-13] + ts + '.itp'

            fullFile = self.main_itp_ras4araf[:-13] + ts + '-PULL' + '.itp'
            if not os.path.isfile(fullFile):
                # PULL file is missing creat
                with open(self.main_itp_ras4araf, "r") as sources:
                    lines = sources.readlines()
                    with open(fullFile, "w") as sources:
                        for line in lines:
                            line = re.sub(r'RAS4A_RAF         1', 'RAS4A_RAF_PULL    1', line)
                            #line = re.sub(r'  414   415      1', ';  414   415      1', line)
                            #line = re.sub(r'  413   414   415      2', ';  413   414   415      2', line)
                            #line = re.sub(r'  414   415   416      2', ';  414   415   416      2', line)
                            line = re.sub(r'  399   400      1', ';  399   400      1', line)
                            line = re.sub(r'  418   419      1', ';  418   419      1', line)
                            line = re.sub(r'  398   399   400      2', ';  398   399   400      2', line)
                            line = re.sub(r'  399   400   401      2', ';  399   400   401      2', line)
                            line = re.sub(r'  417   418   419      2', ';  417   418   419      2', line)
                            line = re.sub(r'  418   419   420      2', ';  418   419   420      2', line)

                            sources.write(line)
                LOGGER.warning("Feedback AA to CG - PULL file for .itp {} not found - create file {}".format(main_itp_final_ras4araf, fullFile))


        # ----------------------------------------------------------------------
        # write the file
        with open(self.main_itp_ras4araf) as fp:
            for line in fp:
                data = line.split()
                if len(data)>1 and data[1] == "Secondary":
                    data_ss = next(fp).split()
                    ss_full = data_ss[1]
                    break


        if hvr_updated:
            filename_hvr = self._create_tmp(self.pdb_full_tar_ras4araf, self.pdb_full_hvr_ras4araf)
            with open(filename_hvr,'r+') as fp:
                lines = fp.readlines()
                lines.insert(0, "REMARK " + "300 " + ss_full + "\n")
                fp.seek(0)
                fp.writelines(lines)

        if crd_updated:
            filename_crd = self._create_tmp(self.pdb_full_tar_ras4araf, self.pdb_full_crd_ras4araf)
            if hvr_updated and filename_hvr != filename_crd:
                with open(filename_crd,'r+') as fp:
                    lines = fp.readlines()
                    lines.insert(0, "REMARK " + "300 " + ss_full + "\n")
                    fp.seek(0)
                    fp.writelines(lines)
            elif not hvr_updated:
                with open(filename_crd,'r+') as fp:
                    lines = fp.readlines()
                    lines.insert(0, "REMARK " + "300 " + ss_full + "\n")
                    fp.seek(0)
                    fp.writelines(lines)


        if hvr_updated and crd_updated:
            pdb_full_hvr_ts = 'RAS4A_RAF_scfix-CYFpos' + '-HVR' + ts + '.pdb'
            pdb_full_crd_ts = 'RAS4A_RAF_scfix-CYFpos' + '-CRD' + ts + '.pdb'
            shutil.copy(self.main_itp_ras4araf, main_itp_final_ras4araf)
            shutil.copy(filename_crd, pdb_full_crd_ts)
            shutil.copy(filename_hvr, pdb_full_hvr_ts)
        elif not hvr_updated and crd_updated:
            pdb_full_crd_ts = 'RAS4A_RAF_scfix-CYFpos' + '-CRD' + ts + '.pdb'
            shutil.copy(self.main_itp_ras4araf, main_itp_final_ras4araf)
            shutil.copy(filename_crd, pdb_full_crd_ts)
        elif hvr_updated and not crd_updated:
            pdb_full_hvr_ts = 'RAS4A_RAF_scfix-CYFpos' + '-HVR' + ts + '.pdb'
            shutil.copy(self.main_itp_ras4araf, main_itp_final_ras4araf)
            shutil.copy(filename_hvr, pdb_full_hvr_ts)
        else:
            LOGGER.info('There is no update in both HVR and CRD regions for RAS4A-RAF')


        if hvr_updated:
            os.remove(filename_hvr)

        if crd_updated:
            if crd_updated and hvr_updated and filename_hvr != filename_crd:
                os.remove(filename_crd)
            elif not hvr_updated:
                os.remove(filename_crd)

        # ----------------------------------------------------------------------
        ######### Helgi's protocol ##########
        # ----------------------------------------------------------------------
        if hvr_updated or crd_updated:

            os.chdir(self.workspace)

            os.mkdir('ddcMDParams' + ts)
            os.chdir('ddcMDParams'+ ts)

            files_copy = Naming.dir_res('ddcMD')+'/martini2objfiles/*'
            for f in glob.glob(files_copy):
                shutil.copy(f,os.getcwd())

            #TODO: Fikret, please check if this needs a so/q;p
            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + main_itp_final_ras4araf + " > temp_RAS4a_RAF.itp")

            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + self.main_itp_latest_ras4b + " > temp_RAS4b.itp")
            
            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + self.main_itp_latest_ras4braf + " > temp_RAS4b_RAF.itp")
            
            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + self.main_itp_latest_ras4a + " > temp_RAS4a.itp")
            

            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4b_RAF.itp > RAS_RAF.itp')

            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4b.itp > RAS.itp')
            
            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4a_RAF.itp > RAS4a_RAF.itp')
            
            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4a.itp > RAS4a.itp')

            os.system('rm temp_RAS4b_RAF.itp temp_RAS4b.itp temp_RAS4a.itp temp_RAS4a_RAF.itp')

            os.system('martini2obj -i itpList -t proItpList -p  martini_v2.1-dna.itp -o martini.data -l molecule.data > itp_to_ddcMD.out')

            os.system("sed -i -n '/moleculeClass MOLECULECLASS/,$p' molecule.data")

        # ----------------------------------------------------------------------
        self.main_itp_latest_ras4araf = main_itp_final_ras4araf
    
        # ras-raf
        crd_updated = False
        hvr_updated = False

    # --------------------------------------------------------------------------
    
    def ras4a_feedback(self):

        LOGGER.info('Doing Feedback for RAS4A')
        hvr_updated = False

        # To update the itp always based on the initial structure (no history of previous updates)
        self.main_itp_ras4a = os.path.join(self.respath,'RAS4A_scfix-CYFpos.itp')

        
        # ----------------------------------------------------------------------
        if self.perc_hvr_ras4a > self.hvr_th_ras4a and self.most_hvr_ras4a != self.prev_hvr_ras4a:

            self.ras4a_feedback_done = True

            self.pdb_most_hvr_ss_ras4a = self.dict_hvr_ss_pdb_ras4a.get(self.most_hvr_ras4a)
            
            out_itp_1_ras4a = 'Protein_hvr1_ras4a'
            self._martinize_with_dssp(self.pdb_hvr_tar_ras4a, self.pdb_most_hvr_ss_ras4a, out_itp_1_ras4a)

            out_itp_1_ras4a = out_itp_1_ras4a + '.itp'
            with open(out_itp_1_ras4a) as fp:
                for line in fp:
                    data = line.split()
                    if len(data) > 1 and data[1] == "Secondary":
                        data_ss = next(fp).split()
                        ss = data_ss[1]
                        break


            ss_revised = ss[0:self.temp_hvr_index_start] + self.most_hvr_ras4a + ss[self.temp_hvr_index_end_ras4a:-1] + ss[-1]
            


            for i in range(0,len(ss_revised)-1):
                resid = i+3
                if (ss_revised[i] == '2' and ss_revised[i+1] == 'C' and resid > 160 and resid < 185):
                    self.hvr_end_main_ras4a = resid


            ss_hvr_revised_fixed = ''

            for i in range(0,len(ss_revised)):
                if ss_revised[i] in ['1', '2']:
                    ss_hvr_revised_fixed = ss_hvr_revised_fixed + 'H'
                else:
                    ss_hvr_revised_fixed = ss_hvr_revised_fixed + ss_revised[i]

            # use new AA structure if helix grows beyond residue 172, otherwise use the initial structure
            print(self.hvr_end_main_ras4a)
            if self.hvr_end_main_ras4a > 172:
                out_itp_2_ras4a = 'Protein_hvr2_ras4a'
                self._martinize_with_ss(self.pdb_hvr_tar_ras4a_elastic, self.pdbkey_hvr_ras4a_elastic, ss_hvr_revised_fixed, out_itp_2_ras4a)
            else:
                out_itp_2_ras4a = 'Protein_hvr2_ras4a'
                self._martinize_with_ss(self.pdb_hvr_tar_ras4a_elastic, self.pdbkey_hvr_ras4a_elastic, ss_hvr_revised_fixed, out_itp_2_ras4a)
                

            out_itp_2_ras4a = out_itp_2_ras4a + '.itp'
            self.revised_itp = out_itp_2_ras4a
            self.main_itp_ras4a = self.update_itp_hvr()
            self.prev_hvr_ras4a = self.most_hvr_ras4a
            self.pdb_full_hvr_ras4a = self.dict_full_ss_pdb_ras4a.get(self.most_hvr_ras4a)
            hvr_updated = True

            self.ras4a_feedback_done = False
            
        else:
            LOGGER.info('Threshold for HVR is not exceeded for RAS4A, continue with current itp file')
            hvr_updated = False



        if not hvr_updated:
            return
 

        if hvr_updated:

            ts = add_timestamp('')

            main_itp_final_ras4a = self.main_itp_ras4a[:-13] + ts + '.itp'

            fullFile = self.main_itp_ras4a[:-13] + ts + '-PULL' + '.itp'
            if not os.path.isfile(fullFile):
                # PULL file is missing creat                                                                                                  
                with open(self.main_itp_ras4a, "r") as sources:
                    lines = sources.readlines()
                    with open(fullFile, "w") as sources:
                        for line in lines:
                            line = re.sub(r'RAS4A         1', 'RAS4A_PULL   1', line)
                            #line = re.sub(r'  414   415      1', ';  414   415      1', line)
                            #line = re.sub(r'  413   414   415      2', ';  413   414   415      2', line)
                            #line = re.sub(r'  414   415   416      2', ';  414   415   416      2', line)
                            line = re.sub(r'  399   400      1', ';  399   400      1', line)
                            line = re.sub(r'  418   419      1', ';  418   419      1', line)
                            line = re.sub(r'  398   399   400      2', ';  398   399   400      2', line)
                            line = re.sub(r'  399   400   401      2', ';  399   400   401      2', line)
                            line = re.sub(r'  417   418   419      2', ';  417   418   419      2', line)
                            line = re.sub(r'  418   419   420      2', ';  418   419   420      2', line)
                            

                            sources.write(line)
                LOGGER.warning("Feedback AA to CG - PULL file for .itp {} not found - create file {}".format(main_itp_final_ras4a, fullFile))



        # ----------------------------------------------------------------------
        # write the file
        with open(self.main_itp_ras4a) as fp:
            for line in fp:
                data = line.split()
                if len(data)>1 and data[1] == "Secondary":
                    data_ss = next(fp).split()
                    ss_full = data_ss[1]
                    break


        if hvr_updated:
            filename_hvr = self._create_tmp(self.pdb_full_tar_ras4a, self.pdb_full_hvr_ras4a)
            with open(filename_hvr,'r+') as fp:
                lines = fp.readlines()
                lines.insert(0, "REMARK " + "300 " + ss_full + "\n")
                fp.seek(0)
                fp.writelines(lines)



        if hvr_updated:
            pdb_full_hvr_ts = 'RAS4A_scfix-CYFpos' + '-HVR' + ts + '.pdb'
            shutil.copy(self.main_itp_ras4a, main_itp_final_ras4a)
            shutil.copy(filename_hvr, pdb_full_hvr_ts)
        else:
            LOGGER.info('There is no update in HVR region for RAS4A')


        if hvr_updated:
            os.remove(filename_hvr)


        # ----------------------------------------------------------------------
        ######### Helgi's protocol ##########
        # ----------------------------------------------------------------------
        if hvr_updated:

            os.chdir(self.workspace)

            os.mkdir('ddcMDParams' + ts)
            os.chdir('ddcMDParams'+ ts)

            files_copy = Naming.dir_res('ddcMD')+'/martini2objfiles/*'
            for f in glob.glob(files_copy):
                shutil.copy(f,os.getcwd())

            #TODO: Fikret, please check if this needs a so/q;p
            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + self.main_itp_latest_ras4braf + " > temp_RAS4b_RAF.itp")
            
            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + self.main_itp_latest_ras4b + " > temp_RAS4b.itp")
            
            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + self.main_itp_latest_ras4araf + " > temp_RAS4a_RAF.itp")
            
            os.system("sed -n '/; Add Z pos res to farnizil/q;p' " + main_itp_final_ras4a + " > temp_RAS4a.itp")

            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4b_RAF.itp > RAS_RAF.itp')
            
            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4b.itp > RAS.itp')
            
            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4a_RAF.itp > RAS4a_RAF.itp')
            
            os.system(Naming.dir_res('ddcMD') + '/dihedrals_itp2ddcmd.awk temp_RAS4a.itp > RAS4a.itp')


            os.system('rm temp_RAS4b_RAF.itp temp_RAS4b.itp temp_RAS4a.itp temp_RAS4a_RAF.itp')

            os.system('martini2obj -i itpList -t proItpList -p  martini_v2.1-dna.itp -o martini.data -l molecule.data > itp_to_ddcMD.out')

            os.system("sed -i -n '/moleculeClass MOLECULECLASS/,$p' molecule.data")

        # ----------------------------------------------------------------------
        
        self.main_itp_latest_ras4a = main_itp_final_ras4a
        
        # ras4a
        hvr_updated = False

    # --------------------------------------------------------------------------
    
    # --------------------------------------------------------------------------
    def checkpoint(self):

        LOGGER.debug('tying to checkpoint: nreports = {} {}'.format(self.nreports, type(self.nreports)))
        LOGGER.debug('tying to checkpoint: naggregates = {} {}'.format(self.naggregates, type(self.naggregates)))
        LOGGER.debug('tying to checkpoint: naggs_at_last_report = {} {}'.format(self.naggs_at_last_report, type(self.naggs_at_last_report)))

        filename = os.path.join(self.workspace, 'checkpoint.npz')
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        data = {'ts': st,
                'data_rmsf_ras4braf': self.aggregate_rmsf_ras4braf,
                'data_rmsf_ras4b': self.aggregate_rmsf_ras4b,
                'data_coord_ras4braf': self.aggregate_coord_ras4braf,
                'data_coord_ras4b': self.aggregate_coord_ras4b,
                'data_hvr_ss_ras4braf': self.ss_hvr_all_ras4braf,
                'data_hvr_ss_ras4b': self.ss_hvr_all_ras4b,
                'data_crd_ss_ras4braf': self.ss_crd_all_ras4braf,
                'previous_hvr_ras4braf': self.prev_hvr_ras4braf,
                'previous_hvr_ras4b': self.prev_hvr_ras4b,
                'previous_crd_ras4braf':self.prev_crd_ras4braf,
                'pdb_full_hvr_ras4braf':self.pdb_full_hvr_ras4braf,
                'pdb_full_hvr_ras4b':self.pdb_full_hvr_ras4b,
                'pdb_full_crd_ras4braf':self.pdb_full_crd_ras4braf,
                'pdb_most_hvr_ss_ras4braf':self.pdb_most_hvr_ss_ras4braf,
                'pdb_most_hvr_ss_ras4b':self.pdb_most_hvr_ss_ras4b,
                'pdb_most_crd_ss_ras4braf':self.pdb_most_crd_ss_ras4braf,
                'hvr_end_main_ras4braf':self.hvr_end_main_ras4braf,
                'hvr_end_main_ras4b':self.hvr_end_main_ras4b,
                'dict_hvr_ss_pdb_ras4braf':self.dict_hvr_ss_pdb_ras4braf,
                'dict_hvr_ss_pdb_ras4b':self.dict_hvr_ss_pdb_ras4b,
                'dict_crd_ss_pdb_ras4braf':self.dict_crd_ss_pdb_ras4braf,
                'dict_full_ss_pdb_ras4braf':self.dict_full_ss_pdb_ras4braf,
                'dict_full_ss_pdb_ras4b':self.dict_full_ss_pdb_ras4b,
                'data_rmsf_ras4araf': self.aggregate_rmsf_ras4araf,
                'data_rmsf_ras4a': self.aggregate_rmsf_ras4a,
                'data_coord_ras4araf': self.aggregate_coord_ras4araf,
                'data_coord_ras4a': self.aggregate_coord_ras4a,
                'data_hvr_ss_ras4araf': self.ss_hvr_all_ras4araf,
                'data_hvr_ss_ras4a': self.ss_hvr_all_ras4a,
                'data_crd_ss_ras4araf': self.ss_crd_all_ras4araf,
                'previous_hvr_ras4araf': self.prev_hvr_ras4araf,
                'previous_hvr_ras4a': self.prev_hvr_ras4a,
                'previous_crd_ras4araf':self.prev_crd_ras4araf,
                'pdb_full_hvr_ras4araf':self.pdb_full_hvr_ras4araf,
                'pdb_full_hvr_ras4a':self.pdb_full_hvr_ras4a,
                'pdb_full_crd_ras4araf':self.pdb_full_crd_ras4araf,
                'pdb_most_hvr_ss_ras4araf':self.pdb_most_hvr_ss_ras4araf,
                'pdb_most_hvr_ss_ras4a':self.pdb_most_hvr_ss_ras4a,
                'pdb_most_crd_ss_ras4araf':self.pdb_most_crd_ss_ras4araf,
                'hvr_end_main_ras4araf':self.hvr_end_main_ras4araf,
                'hvr_end_main_ras4a':self.hvr_end_main_ras4a,
                'dict_hvr_ss_pdb_ras4araf':self.dict_hvr_ss_pdb_ras4araf,
                'dict_hvr_ss_pdb_ras4a':self.dict_hvr_ss_pdb_ras4a,
                'dict_crd_ss_pdb_ras4araf':self.dict_crd_ss_pdb_ras4araf,
                'dict_full_ss_pdb_ras4araf':self.dict_full_ss_pdb_ras4araf,
                'dict_full_ss_pdb_ras4a':self.dict_full_ss_pdb_ras4a,
                'nframes_fb_threshold_ras4a':self.nframes_fb_threshold_ras4a,
                'nframes_fb_threshold_ras4b':self.nframes_fb_threshold_ras4b,
                'nframes_fb_threshold_ras4araf':self.nframes_fb_threshold_ras4araf,
                'nframes_fb_threshold_ras4braf':self.nframes_fb_threshold_ras4braf,
                'main_itp_latest_ras4a':self.main_itp_latest_ras4a,
                'main_itp_latest_ras4b':self.main_itp_latest_ras4b,
                'main_itp_latest_ras4araf':self.main_itp_latest_ras4araf,
                'main_itp_latest_ras4braf':self.main_itp_latest_ras4braf,
                'naggregates': self.naggregates,
                'nreports': self.nreports,
                'naggs_at_last_report': self.naggs_at_last_report}

        self.iointerface.take_backup(filename)
        self.iointerface.save_npz(self.workspace, 'checkpoint', data, interfaces=['simple'])

        LOGGER.debug('after checkpoint: nreports = {} {}'.format(self.nreports, type(self.nreports)))
        LOGGER.debug('after checkpoint: naggregates = {} {}'.format(self.naggregates, type(self.naggregates)))
        LOGGER.debug('after checkpoint: naggs_at_last_report = {} {}'.format(self.naggs_at_last_report, type(self.naggs_at_last_report)))

    def restore(self, filename=''):

        if filename == '':
            filename = os.path.join(self.workspace, 'checkpoint.npz')

        if not os.path.isfile(filename):
            LOGGER.info('Could not find file "%s". Returning.', filename)
            return

        data = self.iointerface.load_npz(self.workspace, 'checkpoint', interfaces = ['simple'])
        if not data or len(data) <= 0:
            LOGGER.info('Could not restore from "%s". Returning.', filename)
            return

        try:
            ts = data['ts']
            self.aggregate_rmsf_ras4braf = data['data_rmsf_ras4braf']
            self.aggregate_rmsf_ras4b = data['data_rmsf_ras4b']
            self.aggregate_coord_ras4braf = data['data_coord_ras4braf']
            self.aggregate_coord_ras4b = data['data_coord_ras4b']
            self.ss_hvr_all_ras4braf = data['data_hvr_ss_ras4braf']
            self.ss_hvr_all_ras4b = data['data_hvr_ss_ras4b']
            self.ss_crd_all_ras4braf = data['data_crd_ss_ras4braf']
            self.prev_hvr_ras4braf = data['previous_hvr_ras4braf']
            self.prev_hvr_ras4b = data['previous_hvr_ras4b']
            self.prev_crd_ras4braf = data['previous_crd_ras4braf']
            self.pdb_full_hvr_ras4braf = data['pdb_full_hvr_ras4braf']
            self.pdb_full_hvr_ras4b = data['pdb_full_hvr_ras4b']
            self.pdb_full_crd_ras4braf = data['pdb_full_crd_ras4braf']
            self.pdb_most_hvr_ss_ras4braf = data['pdb_most_hvr_ss_ras4braf']
            self.pdb_most_hvr_ss_ras4b = data['pdb_most_hvr_ss_ras4b']
            self.pdb_most_crd_ss_ras4braf = data['pdb_most_crd_ss_ras4braf']
            self.hvr_end_main_ras4braf = data['hvr_end_main_ras4braf']
            self.hvr_end_main_ras4b = data['hvr_end_main_ras4b']
            self.dict_hvr_ss_pdb_ras4braf = data['dict_hvr_ss_pdb_ras4braf']
            self.dict_hvr_ss_pdb_ras4b = data['dict_hvr_ss_pdb_ras4b']
            self.dict_crd_ss_pdb_ras4braf = data['dict_crd_ss_pdb_ras4braf']
            self.dict_full_ss_pdb_ras4braf = data['dict_full_ss_pdb_ras4braf']
            self.dict_full_ss_pdb_ras4b = data['dict_full_ss_pdb_ras4b']
            self.aggregate_rmsf_ras4araf = data['data_rmsf_ras4araf']
            self.aggregate_rmsf_ras4a = data['data_rmsf_ras4a']
            self.aggregate_coord_ras4araf = data['data_coord_ras4araf']
            self.aggregate_coord_ras4a = data['data_coord_ras4a']
            self.ss_hvr_all_ras4araf = data['data_hvr_ss_ras4araf']
            self.ss_hvr_all_ras4a = data['data_hvr_ss_ras4a']
            self.ss_crd_all_ras4araf = data['data_crd_ss_ras4araf']
            self.prev_hvr_ras4araf = data['previous_hvr_ras4araf']
            self.prev_hvr_ras4a = data['previous_hvr_ras4a']
            self.prev_crd_ras4araf = data['previous_crd_ras4araf']
            self.pdb_full_hvr_ras4araf = data['pdb_full_hvr_ras4araf']
            self.pdb_full_hvr_ras4a = data['pdb_full_hvr_ras4a']
            self.pdb_full_crd_ras4araf = data['pdb_full_crd_ras4araf']
            self.pdb_most_hvr_ss_ras4araf = data['pdb_most_hvr_ss_ras4araf']
            self.pdb_most_hvr_ss_ras4a = data['pdb_most_hvr_ss_ras4a']
            self.pdb_most_crd_ss_ras4araf = data['pdb_most_crd_ss_ras4araf']
            self.hvr_end_main_ras4araf = data['hvr_end_main_ras4araf']
            self.hvr_end_main_ras4a = data['hvr_end_main_ras4a']
            self.dict_hvr_ss_pdb_ras4araf = data['dict_hvr_ss_pdb_ras4araf']
            self.dict_hvr_ss_pdb_ras4a = data['dict_hvr_ss_pdb_ras4a']
            self.dict_crd_ss_pdb_ras4araf = data['dict_crd_ss_pdb_ras4araf']
            self.dict_full_ss_pdb_ras4araf = data['dict_full_ss_pdb_ras4araf']
            self.dict_full_ss_pdb_ras4a = data['dict_full_ss_pdb_ras4a']
            self.nframes_fb_threshold_ras4a = data['nframes_fb_threshold_ras4a']
            self.nframes_fb_threshold_ras4b = data['nframes_fb_threshold_ras4b']
            self.nframes_fb_threshold_ras4araf = data['nframes_fb_threshold_ras4araf']
            self.nframes_fb_threshold_ras4braf = data['nframes_fb_threshold_ras4braf']
            self.naggregates = data['naggregates']
            self.nreports = data['nreports']
            self.naggs_at_last_report = data['naggs_at_last_report']
            self.main_itp_latest_ras4a = data['main_itp_latest_ras4a']
            self.main_itp_latest_ras4b = data['main_itp_latest_ras4b']
            self.main_itp_latest_ras4araf = data['main_itp_latest_ras4araf']
            self.main_itp_latest_ras4braf = data['main_itp_latest_ras4braf']
            
            LOGGER.info('Restored checkpoint file {} from {}'.format(filename, ts))

            # ras4b-raf
            if isinstance(self.ss_hvr_all_ras4braf, np.ndarray):
                self.ss_hvr_all_ras4braf = self.ss_hvr_all_ras4braf.tolist()
            if isinstance(self.ss_crd_all_ras4braf, np.ndarray):
                self.ss_crd_all_ras4braf = self.ss_crd_all_ras4braf.tolist()
            if isinstance(self.aggregate_rmsf_ras4braf, np.ndarray):
                self.aggregate_rmsf_ras4braf = self.aggregate_rmsf_ras4braf.tolist()
            if isinstance(self.aggregate_coord_ras4braf, np.ndarray):
                self.aggregate_coord_ras4braf = self.aggregate_coord_ras4braf.tolist()

            # ras4b
            if isinstance(self.ss_hvr_all_ras4b, np.ndarray):
                self.ss_hvr_all_ras4b = self.ss_hvr_all_ras4b.tolist()
            if isinstance(self.aggregate_rmsf_ras4b, np.ndarray):
                self.aggregate_rmsf_ras4b = self.aggregate_rmsf_ras4b.tolist()
            if isinstance(self.aggregate_coord_ras4b, np.ndarray):
                self.aggregate_coord_ras4b = self.aggregate_coord_ras4b.tolist()
                
            # ras4a-raf
            if isinstance(self.ss_hvr_all_ras4araf, np.ndarray):
                self.ss_hvr_all_ras4araf = self.ss_hvr_all_ras4araf.tolist()
            if isinstance(self.ss_crd_all_ras4araf, np.ndarray):
                self.ss_crd_all_ras4araf = self.ss_crd_all_ras4araf.tolist()
            if isinstance(self.aggregate_rmsf_ras4araf, np.ndarray):
                self.aggregate_rmsf_ras4araf = self.aggregate_rmsf_ras4araf.tolist()
            if isinstance(self.aggregate_coord_ras4araf, np.ndarray):
                self.aggregate_coord_ras4araf = self.aggregate_coord_ras4araf.tolist()

            # ras4a
            if isinstance(self.ss_hvr_all_ras4a, np.ndarray):
                self.ss_hvr_all_ras4a = self.ss_hvr_all_ras4a.tolist()
            if isinstance(self.aggregate_rmsf_ras4a, np.ndarray):
                self.aggregate_rmsf_ras4a = self.aggregate_rmsf_ras4a.tolist()
            if isinstance(self.aggregate_coord_ras4a, np.ndarray):
                self.aggregate_coord_ras4a = self.aggregate_coord_ras4a.tolist()    
            
            #--------------------------------------------------------------------------            

            # ras4b-raf
            if isinstance(self.dict_hvr_ss_pdb_ras4braf, np.ndarray):
                self.dict_hvr_ss_pdb_ras4braf = self.dict_hvr_ss_pdb_ras4braf.item()
            if isinstance(self.dict_crd_ss_pdb_ras4braf, np.ndarray):
                self.dict_crd_ss_pdb_ras4braf = self.dict_crd_ss_pdb_ras4braf.item()
            if isinstance(self.dict_full_ss_pdb_ras4braf, np.ndarray):
                self.dict_full_ss_pdb_ras4braf = self.dict_full_ss_pdb_ras4braf.item()

            # ras4b
            if isinstance(self.dict_hvr_ss_pdb_ras4b, np.ndarray):
                self.dict_hvr_ss_pdb_ras4b = self.dict_hvr_ss_pdb_ras4b.item()
            if isinstance(self.dict_full_ss_pdb_ras4b, np.ndarray):
                self.dict_full_ss_pdb_ras4b = self.dict_full_ss_pdb_ras4b.item()

            # ras4a-raf
            if isinstance(self.dict_hvr_ss_pdb_ras4araf, np.ndarray):
                self.dict_hvr_ss_pdb_ras4araf = self.dict_hvr_ss_pdb_ras4araf.item()
            if isinstance(self.dict_crd_ss_pdb_ras4araf, np.ndarray):
                self.dict_crd_ss_pdb_ras4araf = self.dict_crd_ss_pdb_ras4araf.item()
            if isinstance(self.dict_full_ss_pdb_ras4araf, np.ndarray):
                self.dict_full_ss_pdb_ras4araf = self.dict_full_ss_pdb_ras4araf.item()

            # ras4a
            if isinstance(self.dict_hvr_ss_pdb_ras4a, np.ndarray):
                self.dict_hvr_ss_pdb_ras4a = self.dict_hvr_ss_pdb_ras4a.item()
            if isinstance(self.dict_full_ss_pdb_ras4a, np.ndarray):
                self.dict_full_ss_pdb_ras4a = self.dict_full_ss_pdb_ras4a.item()

            #--------------------------------------------------------------------------
            
            if isinstance(self.nreports, np.ndarray):
                self.nreports = int(self.nreports)
            if isinstance(self.naggregates, np.ndarray):
                self.naggregates = int(self.naggregates)
            if isinstance(self.naggs_at_last_report, np.ndarray):
                self.naggs_at_last_report = int(self.naggs_at_last_report)

            #--------------------------------------------------------------------------
            if isinstance(self.nframes_fb_threshold_ras4a, np.ndarray):
                self.nframes_fb_threshold_ras4a = int(self.nframes_fb_threshold_ras4a)
            if isinstance(self.nframes_fb_threshold_ras4b, np.ndarray):
                self.nframes_fb_threshold_ras4b = int(self.nframes_fb_threshold_ras4b)
            if isinstance(self.nframes_fb_threshold_ras4araf, np.ndarray):
                self.nframes_fb_threshold_ras4araf = int(self.nframes_fb_threshold_ras4araf)
            if isinstance(self.nframes_fb_threshold_ras4braf, np.ndarray):
                self.nframes_fb_threshold_ras4braf = int(self.nframes_fb_threshold_ras4braf)

            #--------------------------------------------------------------------------    

            # ras4b-raf
            if isinstance(self.hvr_end_main_ras4braf, np.ndarray):
                self.hvr_end_main_ras4braf = int(self.hvr_end_main_ras4braf)

            # ras4b
            if isinstance(self.hvr_end_main_ras4b, np.ndarray):
                self.hvr_end_main_ras4b = int(self.hvr_end_main_ras4b)
                
            # ras4a-raf
            if isinstance(self.hvr_end_main_ras4araf, np.ndarray):
                self.hvr_end_main_ras4araf = int(self.hvr_end_main_ras4araf)

            # ras4a
            if isinstance(self.hvr_end_main_ras4a, np.ndarray):
                self.hvr_end_main_ras4a = int(self.hvr_end_main_ras4a)    
                               
            #--------------------------------------------------------------------------    

            # ras4b-raf
            if isinstance(self.pdb_most_hvr_ss_ras4braf, np.ndarray):
                self.pdb_most_hvr_ss_ras4braf = str(self.pdb_most_hvr_ss_ras4braf)
            if isinstance(self.pdb_most_crd_ss_ras4braf, np.ndarray):
                self.pdb_most_crd_ss_ras4braf = str(self.pdb_most_crd_ss_ras4braf)
            if isinstance(self.prev_hvr_ras4braf, np.ndarray):
                self.prev_hvr_ras4braf = str(self.prev_hvr_ras4braf)
            if isinstance(self.prev_crd_ras4braf, np.ndarray):
                self.prev_crd_ras4braf = str(self.prev_crd_ras4braf)
            if isinstance(self.pdb_full_hvr_ras4braf, np.ndarray):
                self.pdb_full_hvr_ras4braf = str(self.pdb_full_hvr_ras4braf)
            if isinstance(self.pdb_full_crd_ras4braf, np.ndarray):
                self.pdb_full_crd_ras4braf = str(self.pdb_full_crd_ras4braf)
            if isinstance(self.main_itp_latest_ras4braf, np.ndarray):
                self.main_itp_latest_ras4braf = str(self.main_itp_latest_ras4braf)

            # ras4b
            if isinstance(self.pdb_most_hvr_ss_ras4b, np.ndarray):
                self.pdb_most_hvr_ss_ras4b = str(self.pdb_most_hvr_ss_ras4b)
            if isinstance(self.prev_hvr_ras4b, np.ndarray):
                self.prev_hvr_ras4b = str(self.prev_hvr_ras4b)
            if isinstance(self.pdb_full_hvr_ras4b, np.ndarray):
                self.pdb_full_hvr_ras4b = str(self.pdb_full_hvr_ras4b)
            if isinstance(self.main_itp_latest_ras4b, np.ndarray):
                self.main_itp_latest_ras4b =str(self.main_itp_latest_ras4b)
                
            # ras4a-raf
            if isinstance(self.pdb_most_hvr_ss_ras4araf, np.ndarray):
                self.pdb_most_hvr_ss_ras4araf = str(self.pdb_most_hvr_ss_ras4araf)
            if isinstance(self.pdb_most_crd_ss_ras4araf, np.ndarray):
                self.pdb_most_crd_ss_ras4araf = str(self.pdb_most_crd_ss_ras4araf)
            if isinstance(self.prev_hvr_ras4araf, np.ndarray):
                self.prev_hvr_ras4araf = str(self.prev_hvr_ras4araf)
            if isinstance(self.prev_crd_ras4araf, np.ndarray):
                self.prev_crd_ras4araf = str(self.prev_crd_ras4araf)
            if isinstance(self.pdb_full_hvr_ras4araf, np.ndarray):
                self.pdb_full_hvr_ras4araf = str(self.pdb_full_hvr_ras4araf)
            if isinstance(self.pdb_full_crd_ras4araf, np.ndarray):
                self.pdb_full_crd_ras4araf = str(self.pdb_full_crd_ras4araf)
            if isinstance(self.main_itp_latest_ras4araf, np.ndarray):
                self.main_itp_latest_ras4araf =str(self.main_itp_latest_ras4araf)

            # ras4a
            if isinstance(self.pdb_most_hvr_ss_ras4a, np.ndarray):
                self.pdb_most_hvr_ss_ras4a = str(self.pdb_most_hvr_ss_ras4a)
            if isinstance(self.prev_hvr_ras4a, np.ndarray):
                self.prev_hvr_ras4a = str(self.prev_hvr_ras4a)
            if isinstance(self.pdb_full_hvr_ras4a, np.ndarray):
                self.pdb_full_hvr_ras4a = str(self.pdb_full_hvr_ras4a)
            if isinstance(self.main_itp_latest_ras4a, np.ndarray):
                self.main_itp_latest_ras4a =str(self.main_itp_latest_ras4a)
                
            #--------------------------------------------------------------------------    


        except Exception:
            LOGGER.error('Found corrupt checkpoint file!')
            return

    @staticmethod
    def save_tuple(io, namespace, tuple):
        filenames = [t[0] for t in tuple]
        data = [str.encode(t[1]) for t in tuple]
        io.save_files(namespace, filenames, data)

    @staticmethod
    def test(path, key):
        return True
