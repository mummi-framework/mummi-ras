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

import numpy as np
import os
from logging import getLogger
import torch
import yaml
import time
from typing import List, Union, Any, Dict
import glob
import MDAnalysis as mda
from scipy.spatial import distance_matrix
from operator import itemgetter

from .generator import Generator
from ..models.model_regression_torch import CustomRegressionModelBackmapping as UCGEncoder

# ------------------------------------------------------------------------------
LOGGER = getLogger(__name__)
FRAMEWORK = "torch"
# ------------------------------------------------------------------------------

class UCGGenerator(Generator):
    """
    Generate protein structures from ultra-caorse grained (UCG) simulations (LAMMPS).

    :param trajectories: List of path to trajectory files in the .xyz format
    :type trajectories: List[str]
    :param data_path: Path to backmapping data
    :type data_path: str
    :param model_path: Path to the ML model (a file `config.yaml` must be present)
    :type model_path: str
    :param outpath: Directory where patches will be written out
    :type outpath: str
    """
    def __init__(self, iteration_id: int, trajectory_path: str, data_path: str, model_path: str, outpath: str) -> None:
        super().__init__()
        self.framework = FRAMEWORK
        is_framework_supported = self.setup()
        assert is_framework_supported, f"Invalid framework: {self.framework}"
        self.model_path = model_path
        self.config_fname = os.path.join(self.model_path, "config.yaml")
        with open(self.config_fname, 'r') as data:
            self.config = yaml.load(data, Loader = yaml.Loader)
        self.config["model"]["name"] = os.path.basename(self.model_path)

        self.outpath = outpath
        self.iteration_id_str = "iter{:02d}".format(iteration_id)
        self.data_path = data_path
        self.trajectory_path = trajectory_path
        self.inshape = tuple(self.config["data"]["inshape"])
        self.outshape = tuple(self.config["model"]["z_dim"])
        self.config["model"]["type"] = "xyz_reconstruction"
        self.config["model"]["inshape"] = self.inshape
        self.config["outpath"] = self.outpath

        LOGGER.info(f"> Initialization of generator with")
        LOGGER.info(f"   > trajectories                 = {self.trajectory_path}")
        LOGGER.info(f"   > UCG datapath                 = {self.data_path}")
        LOGGER.info(f"   > encoder                      = {self.model_path}")
        LOGGER.info(f"   > dimensions(input, output)    = ({self.inshape}, {self.outshape})")
        LOGGER.info(f"   > outpath                      = {self.outpath}")
        LOGGER.info(f"   > iter                         = {self.iteration_id_str}")

        self.load_model()
        # Store the last position read for each trajectory ID
        self.last_pos_read: Dict[str,str] = {}

        # Define the number of beads per UCG structure (per frame)
        self.n_beads = 40
        # Number of elements of the strictly upper triangular distance matrix
        self.num_elem_upper = 780

        self.trajectories: List[Union[str,str,str]] = self.read_trajectories_file(self.trajectory_path)
        LOGGER.info(f"   > Found {len(self.trajectories)} running UCG simulations in {self.trajectory_path}")
        retry = 5
        while len(self.trajectories) == 0 and retry >= 0:
            LOGGER.warning(f"No UCG simulations are running. UCGServer will fail. Waiting and retry ...")
            retry -= 1
            time.sleep(30)
            self.trajectories: List[Union[str,str,str]] = self.read_trajectories_file(self.trajectory_path)
            LOGGER.info(f"   > Found {len(self.trajectories)} running UCG simulations in {self.trajectory_path}")

        if len(self.trajectories) == 0:
            LOGGER.critical(f"No UCG simulations are running. UCGServer will fail. Aborting ...")
        
        # Contains the latest frame ID generated in ucg/iterXX/structure_iterXX_trajYY_f00000000K for each trajectory file
        self.frame_id: Dict[str,int] = {}

    # ------------------------------------------------------------------------------
    def read_trajectories_file(self, path: str) -> List[Union[str,str,str]]:
        """
        Read a file listing information (name and path to trajectory file) running `k` UCG simulations.
        This file called ucg_trajectories.txt located in mummi_ras.Naming.dir_root('ucg') and
        follows this pattern:
            name_ucgsim_1 seed_1 last_pos_read_1 traj_path_1
            ...
            name_ucgsim_K seed_K last_pos_read_1 traj_path_K

        where last_pos_read_K represents the number of bytes read from file traj_path_K.

        :param path: Path to file listing different paths to trajectories generated by UCG simulations
        :type path: str
        :return: 
        :rtype: List[Union[str,str]]
        """
        res = []
        with open(path, 'r') as f:
            for line in f:
                id_str = line.split(' ')[0].rstrip()
                # seed = line.split(' ')[1].rstrip()
                read_pos = line.split(' ')[2].rstrip()
                self.last_pos_read[id_str] = read_pos
                path = line.split(' ')[-1].rstrip()
                # id = ''.join([x for x in id_str if x.isdigit()])
                res.append((id_str, read_pos, path))
        return res

    # ------------------------------------------------------------------------------
    def check_existing_ucg_structures(self, trajectory_ids: List[Union[str,str]]) -> Dict[str,int]:
        """
        Check if we have existing structures for each trajectory that are running.
        Example: If frame ID for trajectory X is 15 then we want to write the first structure 
        with frame ID 16 to not erase previously generated structure.

        Return a hashmap that give the ID of the first frame that we should generate for 
        each trajectory ID (e.g., ID of the UCG simulation):
            - by default 0 if no UCG structure has been generated for a given trajectory ID
            - otherwise, we look into the outpath to find the latest generated frame ID 
            for a given trajectory ID and we add one

        :param trajectory_ids: A list of tuple (trajectory ID, trajectory file path)
        :type trajectory_ids: List[Union[str,str]]

        :return: Hashmap of the largest existing frame ID for each UCG simulations
        :rtype: Dict[str,int]
        """
        res = {}
        for traj_id, _, _ in trajectory_ids:
            res[traj_id] = self._latest_frame_generated(os.path.join(self.outpath, self.iteration_id_str), traj_id)
        return res

    # ------------------------------------------------------------------------------
    def _latest_frame_generated(self, outpath: str, traj_id: str) -> int:
        """
        Find in the directory outpath containing already some structures
        like structure_iterXX_trajYY_f00000000K.
        - XX: iteration ID
        - YY: trajectory ID
        - K: Structure ID (frame number)
        This function finds the highest K for a given YY (traj_id).
        """
        traj_id = ''.join([x for x in str(traj_id) if x.isdigit()])
        try:
            traj_id = int(traj_id)
        except ValueError as e:
            return 0

        samples = glob.glob(outpath + '/structure*.npz')
        nsamples = len(samples)
        LOGGER.info(f"#samples = {len(samples)}")

        if nsamples == 0:
            return 0
        samples = [x.split('/')[-1].split('.')[0] for x in samples]
        # List of list, each sublist is of the form ['iterX', 'trajY', 'f000000000K']
        samples = [x.split('_')[1:] for x in samples]
        clean_samples = []
        try:
            # List of list, each sublist is of the form [Y, K]
            for elem in samples:
                clean_samples.append([int(self._remove_but_digits(x)) for x in elem])
        except ValueError as e:
            LOGGER.error(f"{e}")
            return 0

        traj_samples = [x for x in clean_samples if x[1] == traj_id]
        if len(traj_samples) == 0:
            return 0

        LOGGER.debug(f"for {traj_id} = traj_samples = {samples}")

        # find the max [Y, K] (maximum K)
        index_frame = -1 # The frame ID is the last element of the list
        max_frame_id = max(traj_samples, key=itemgetter(index_frame))
        return max_frame_id[index_frame]

    def setup(self) -> bool:
        """Setup the ML framework used.

        :return: True if torch is supported
        :rtype: bool
        """
        if self.framework == "torch":
            LOGGER.info(f" > Using torch v{torch.__version__} ({torch.__file__})")
            return True
        return False

    def load_model(self) -> bool:
        """Load the ML model and setup the devices, use GPUs if possible.

        :raises ValueError: if ML framework is not supported

        :return: True if the model is loaded
        :rtype: bool
        """
        self.model = UCGEncoder(
            inshape = self.inshape,
            outshape = self.outshape, 
            config = self.config,
            tstamp = '0',
            create_subdirs = False
        )

        if self.framework == "torch":
            # TODO: this should be computed while setting GPUs
            # TODO: ideally we should remove hard-coded paths
            device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
            # device = torch.device('cpu')
            LOGGER.info(f" > Device used : {device}")
            if torch.cuda.is_available():
                LOGGER.info(f"   > CUDA total devices         : {torch.cuda.device_count()}")
                LOGGER.info(f"   > CUDA active device         : device:{torch.cuda.current_device()}")
                if torch.backends.cudnn.is_available() and torch.backends.cudnn.enabled:
                    LOGGER.info(f"   > cudnn version              : {torch.backends.cudnn.version()}")

            self.model = self.model.to(device)
            self.model.device = device   # so the model can move the data to device later
            model_name = os.path.join(self.model_path, "model.pth")
            std_fname = os.path.join(self.model_path, "UCG-data_summary.npz")
            self.model.load_model(model_name, std_fname) #, map_location = device)
        else:
            ValueError(f"Only PyTorch is supported, {self.framework} is not supported")
        return True

    def get_atoms(self, mdu, residxs):
        resids = " ".join([str(mdu.atoms.types[0 + i]) for i in residxs])
        return mdu.select_atoms(f"type {resids}")

    # get periodic displacement vector
    def get_displacement_2d(self, mdu, a, b):
        d = b - a
        for i in range(2):
            if d[i] > 0.5 * mdu.dimensions[i]:
                d[i] -= mdu.dimensions[i]
            if d[i] < -0.5 * mdu.dimensions[i]:
                d[i] += mdu.dimensions[i]
        return d

    def get_anchors(self, anchor_origin, anchor_xaxis, anchor_xyaxis):
        xy = None
        if anchor_xyaxis is not None:
            xy = anchor_xyaxis.center_of_mass(unwrap=True)

        return (
            anchor_origin.center_of_mass(unwrap=True),
            anchor_xaxis.center_of_mass(unwrap=True),
            xy,
        )

    # --------------------------------------------------------------------------
    def align(self, positions):
        residxs_origin = np.arange(0, 13, dtype=int)  # RESIDXS_GDOMAIN
        residxs_xaxis = np.arange(32, 34, dtype=int)  # RESIDXS_CRD_4B
        residxs_xyaxis = np.array([0])

        ras_raf_state_data = os.path.join(self.data_path, "ras_raf_state_ma_0.data")
        mdu = mda.Universe(ras_raf_state_data)
        mdu.atoms.positions = positions

        anchor_origin = self.get_atoms(mdu, residxs_origin)
        anchor_xaxis = self.get_atoms(mdu, residxs_xaxis)
        anchor_xyaxis = self.get_atoms(mdu, residxs_xyaxis)
        zaxis = np.array([0.0, 0.0, 1.0])

        # ----------------------------------------------------------------------
        # variable to hold return angle
        angle = None
        # membrane alignment
        # ------------------------------------------------------------------
        # step 1. bring anchor 1 to the center (GDomain)
        x0, x1, _ = self.get_anchors(anchor_origin, anchor_xaxis, anchor_xyaxis)
        ref_crd = mdu.dimensions[0:3] / 2
        dx = ref_crd - x0
        mdu.atoms.translate(dx)

        # ------------------------------------------------------------------
        # step 2. bring anchor 2 to +x axis (CRD_4B)
        x0, x1, _ = self.get_anchors(anchor_origin, anchor_xaxis, anchor_xyaxis)

        dvec = self.get_displacement_2d(mdu, x0, x1)
        angle = np.rad2deg(np.arctan2(dvec[1], dvec[0]))

        mdu.atoms.rotateby(-angle, zaxis, point=x0)

        # now move
        x0, x1, _ = self.get_anchors(anchor_origin, anchor_xaxis, anchor_xyaxis)
        dx = -x0
        dx[2] = 0  # we dont want to move z
        mdu.atoms.translate(dx)

        # print(mdu.atoms.positions)
        new_positions = mdu.atoms.positions
        # CYF SC is located in the 20th UCG bead  (19 because the index starts from 0)
        # 413 bead id in CG
        new_positions[:, 2] -= new_positions[np.array([19]), 2]
        return new_positions

    # --------------------------------------------------------------------------
    def standardize(self, data, label, dmean = None, dstd = None):
        LOGGER.debug(f"Standardizing {data.shape} {label}")
        data -= dmean
        data /= dstd
        return data, dmean, dstd

    def destandardize(self, data, label, dmean = None, dstd = None):
        """ Destandardize an array """
        LOGGER.debug(f"De-standardizing {data.shape} {label}")
        data *= dstd
        data += dmean
        return data

    def filter_traj_xyz(self, fname: str, last_position: int = 0, retry: int = 10, sleep_duration: int = 5, max_traj_load: int = None) -> Union[List[str], int]:
        """
        Filter out content from trajectory files using .xyz format. 
        :param fname: File name
        :type fname: str
        :param last_position: Last position read in that trajectory file (where to start reading)
        :type last_position: int, optional
        :param retry: Number of retries
        :type retry: int, optional
        :param sleep_duration: Sleep duration in secs between each retry
        :type sleep_duration: int, optional
        :param max_traj_load: maximum trajectories to load from that file (None means all trajectories)
        :type max_traj_load: int, optional

        :return: List that contains each line of the file which is a valid position and 
        :rtype: Union[List[str], int]
        """
        if retry <= 0:
            return
        i = 0
        stream = []
        last_pos_file = 0
        while i < retry:
            stream, last_pos_file = self._filter_traj_xyz(fname, last_position, max_traj_load)
            if len(stream) > 0:
                break
            i += 1
            time.sleep(sleep_duration)
        return stream, last_pos_file

    def _filter_traj_xyz(self, fname: str, last_position: int = 0, max_traj_load: int = None) -> Union[List[str], int]:
        """
        Internal function to filter out content from trajectory files using .xyz format.
        :param fname: File name
        :type fname: str
        :param last_position: Last position read in that trajectory file (where to start reading)
        :type last_position: int, optional
        :param max_traj_load: maximum trajectories to load from that file (None means all trajectories)
        :type max_traj_load: int, optional

        :return: List that contains each line of the file which is a valid position and 
        :rtype: Union[List[str], int]
        """

        if not max_traj_load:
            max_traj_load = float('inf')

        data = []
        new_last_position = 0
        current_traj = 0
        with open(fname) as f:
            # Jump to the part of the file not read yet
            f.seek(last_position)
            while line := f.readline():
                if not line.strip("\n").startswith("Atoms") and not line.startswith("40\n") and not line.startswith(" generated by VMD\n"):
                    data.append(line)
                if line.startswith("40 "):
                    # we have read one trajectory 
                    current_traj += 1
                    new_last_position = f.tell()
                if current_traj >= max_traj_load:
                    break
            # Store new position in file for the next call of that funcion
            print(f"Reading file {fname} starting from position {last_position}. New last position is {new_last_position}. max_traj_load = {max_traj_load}")
        return data, new_last_position

    def load_trajectories(self, max_traj_load: int = None) -> tuple[np.array, np.array]:
        """
        Read LAMMPS trajectories defined in :py:attr:`.trajectories`.

        :param max_traj_load: The maximum number of position to read. The number of trajectories
        can be <= max_traj_load.
        :type max_traj_load: int, optional
        
        :return: The trajectory IDs and the raw trajectories data for each ID
        :rtype: (np.array, np.array)
        """
        n = len(self.trajectories)
        raw_data = []
        trajectories_ids = []

        if isinstance(max_traj_load, int):
            nb_traj_to_load = max_traj_load // n + (max_traj_load % n)
        else:
            nb_traj_to_load = None # We read all trajectories

        # We read the existing structures in ucg/iterXX to find the latest one generated so 
        # we do not erase existing one
        self.frame_id: Dict[str,int] = self.check_existing_ucg_structures(trajectory_ids = self.trajectories)
        LOGGER.info(f"   > Latest existing frame IDs for UCG structures generated for each trajectory = {self.frame_id}")

        for tridx in range(n):
            trajectoryid = self.trajectories[tridx][0]
            read_pos = int(self.trajectories[tridx][1])
            traj_file = self.trajectories[tridx][-1]

            LOGGER.info(f"Reading {nb_traj_to_load} trajectories from file {traj_file}")
            # Stream contains the data and last_pos_file the position where we stop reading from traj_file
            stream, last_pos_file = self.filter_traj_xyz(fname = traj_file, last_position = read_pos, retry = 10, sleep_duration = 5, max_traj_load = nb_traj_to_load)
            # If a trajectory file is empty (could happen if a UCG simulation crashed)
            if len(stream) == 0:
                continue
            npz = np.loadtxt(stream)
            # Each trajectory has 40 beads (40 lines of 3D coordinates) and a trajectory file
            # contains an unknown number of trajectories (called frames)
            # npz contains n_timesteps * n_beads lines (in this case 40) with 3D coordinates
            # 1 x_1 y_1 z_1
            # ...
            # 40 x_40 y_40 z_40
            # ...
            # 1 x_1 y_1 z_1
            # ...
            # 40 x_40 y_40 z_40
            # ...
            # We remove the first column (line number within each trajectory)
            data = npz[:, 1:]

            LOGGER.info(f"     Data trajectory                                   : {data.shape}")
            n_timesteps = int(data.shape[0] / self.n_beads)
            ucg_positions = np.empty((n_timesteps,  self.n_beads, 3), dtype=np.float32)
            trajectoryid_t = np.empty(n_timesteps, dtype='object')

            LOGGER.debug(f"     Trajectory id                                    : {trajectoryid}")
            for i in range(n_timesteps):
                # We break down data per frame ucg_positions[i] contains coordinates for frame i
                ucg_positions[i, :, :] = data[i*self.n_beads : (i*self.n_beads)+self.n_beads, :]
                # We store a unique ID (we remove the letter from the trajectoryid and keep only digits)
                tid = "{:03d}".format(int(''.join([x for x in trajectoryid if x.isdigit()])))
                structure_id = "{:012d}".format(self.frame_id[trajectoryid]+i+1)
                trajectoryid_t[i] = tid + "_f" + str(structure_id)

            LOGGER.info(f"     UCG positions (#frames, #beads, #dimensions)      : {ucg_positions.shape}")
            raw_data.append(ucg_positions)
            trajectories_ids.append(trajectoryid_t)
            # Store last_pos_file to update ucg_trajectories.txt
            self.last_pos_read[trajectoryid] = last_pos_file
            self.update_pos_traj_file(self.trajectory_path, trajectoryid, last_pos_file)

        raw_data = np.concatenate(raw_data, axis=0)
        trajectories_ids = np.concatenate(trajectories_ids, axis=0)
        return trajectories_ids, raw_data

    def update_pos_traj_file(self, fname: str, traj_id: str, new_pos: str):
        """
        Update the last position read for a given trajectory 
        in the file summarizing the trajectories generated by 
        UCG simulations. This file is (usually) named `ucg_trajectories.txt`.

        :param fname: The name of the file
        :type fname: str
        :param traj_id: The ID of the trajectory, for example "ucg.1" 
        :type traj_id: str
        :param new_pos: The new position in the file (byte position, like the one given by tell())
        :type new_pos: str
        """
        new_content =[]
        with open(fname) as f:
            for line in f:
                data = line.split(' ')
                if data[0] == traj_id:
                    new_line = f"{data[0]} {data[1]} {new_pos} {data[3]}"
                    new_content.append(new_line)
                else:
                    new_content.append(line)

        with open(fname, 'w') as f:
            f.writelines(new_content)

        LOGGER.info(f"Updated file {fname}, traj_id = {traj_id} with new_pos = {new_pos}")

    def process_trajectories(self, raw_data) -> np.array:
        """
        Process the raw trajectories, align the data, standardize etc.
        It returns the processed data that can be sent to the model.

        :param raw_data: The raw trajectories returned by :py:meth:`.load_trajectories`
        :type raw_data: np.array
        :return: The processed trajectories
        :rtype: np.array 
        """
        LOGGER.debug(f"Total ucg_positions : {raw_data.shape}")
        # Apply alignment
        for i in range(raw_data.shape[0]):
            raw_data[i] = self.align(raw_data[i])
        # Apply standardization
        data, _, _ = self.standardize(
            raw_data, "ucg_positions", 
            self.model.datamean, self.model.datastd
        )

        n_patches = data.shape[0]
        n_beads = data.shape[1]
        assert self.n_beads == n_beads, "Beads should be the same, data might be corrupted"
        # # Number of elements of the strictly upper triangular distance matrix
        # self.num_elem_upper = 780
        dist_ms = np.empty((n_patches, self.num_elem_upper), dtype=np.float32)
        flatten_data = np.empty((n_patches, n_beads * 3), dtype=np.float32)
        # for each patch
        for pidx in range(n_patches):
            patch_pos_ucg = data[pidx]
            distance_m = distance_matrix(patch_pos_ucg, patch_pos_ucg)
            dist_ms[pidx] = distance_m[np.triu_indices(n_beads, 1)]
            flatten_data[pidx] = patch_pos_ucg.flatten()
        
        data = np.concatenate((flatten_data, dist_ms), axis=1)
        return data

    def new_positions(self, coords: np.array) -> np.array:
        """
        Generate new positions from UCG trajectories.

        :param coords: The coordinates (positions) of the initial data (either Latent space or UCG simulations)
        :type coords: np.array

        :return: The new positions
        :rtype: np.array 
        """
        # trajectories_ids, raw_data = self.load_trajectories()
        # ucg_coords = self.process_trajectories(raw_data)
        new_positions = self.model.predict(coords)
        LOGGER.info(f"new_positions= {new_positions.shape} from ucg_coords = {coords.shape}")
        return new_positions

    def write_patches(self, iteration_id: int, new_positions: np.array, trajectories_ids: np.array) -> Union[List[str], List[Any]]:
        """
        From positions, write a patch as a npz file and  return its path.

        :param iteration_id: The current iteration ID (which will define the name of the structures generated)
        :type iteration_id: int
        :param new_positions: The new positions that will be written out
        :type new_positions: np.array
        :param trajectories_ids: The IDs of each trajectory, usually from 1 to K (e.g., from which UCG simulations these positions come from)
        :type trajectories_ids: np.array

        :return: The new positions
        :rtype: Union[List[str], List[Any]]
        """
        start = time.time()
        iteration_id = "{:02d}".format(iteration_id)
        RPATH = os.path.join(self.outpath, f"iter{iteration_id}")
        if not os.path.isdir(RPATH):
            os.makedirs(RPATH, exist_ok=True)

        LOGGER.info(f"Write UCG patches in {RPATH}")
        npatches = new_positions.shape[0]
        positionsArray = []
        namesArray = []

        # We want to avoid erasing previously generated structures
        # offset = self._latest_patch_id_generated(RPATH)

        # Create an .npz for each new patch
        for sidx in range(npatches):
            positions = new_positions[sidx]
            traj_id = "{}".format(trajectories_ids[sidx])
            # structure_id = "{:012d}".format(offset+sidx)
            # Structure ID is encoded into traj_id for UCG
            structure_name = f"structure_iter{iteration_id}_traj{traj_id}"
            outfile_positions = os.path.join(RPATH, f"{structure_name}.npz")
            LOGGER.debug(f"Writing {outfile_positions}")
            np.savez_compressed(outfile_positions, data = positions)
            positionsArray.append(positions)
            namesArray.append(structure_name)

        end = time.time() - start
        LOGGER.info(f"Wrote {npatches} patches in {end:.3f} secondes")
        return namesArray, positionsArray

    def _remove_but_digits(self, s: str) -> str:
        """ Remove everything but digits from a string"""
        return ''.join([x for x in s if x.isdigit()])

    def number_struct_generated(self, RPATH: str) -> tuple[int,int]:
        """
        return the total number of structures generated, and the number of valid one
        """
        nsamples = len(glob.glob(RPATH + '/structure*.npz'))
        nsamples_valid = len(glob.glob(RPATH + '/structure*.gro'))
        # Take max to avoid returning a negative number if RPATH was empty
        return max(nsamples, 0), max(nsamples_valid, 0)

