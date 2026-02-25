"""
Sampler class for drawing new points from the latent space.
"""

import collections
import random
import os
import re
import time
from typing import List, Union, Type
import numpy.typing as npt

import numpy as np

import mummi_core
from mummi_core.utils import Naming
import mummi_ras

from logging import getLogger
LOGGER = getLogger(__name__)

class LSSampler_C3:
    """
    LSSampler does the following tasks:
        * encapuslates a generative model for drawing new proposals (currently based on optimal transport interpolation with feedback from validation). This includes:
            1) IO routines for initializing an OT schedule based on A and B sets (reads the database for thsi)
            2) routines for building or rebuilding a OT schedule based on new definitions of A and B and/or new validation information.
            3) routines for sampling new points from the OT schedule and storing those points in a FIFO queue for the ML server (or any other reasonable data structure...)
        * queries the database to get feedback on validation status of points
        * updates proposal distribution based on feedback via validation data in database (feedback loop)
        * maintains a queue of new LS points to pass to the ML server or directly to validation code
        * has a user-facing function to return the first k new LS points from this list

    Contents of feedback file:
    1. structure_names[i] contains the name of the structure i
    2. v_coords[i] contains the coords (4 floats) for structure i: shape of one element (4,)
    3. v2ucg_coords[i] contains the coords (750 floats) for structure i: shape of one element (250, 3)
    4. ucg2cg_coords[i] contains the coords (14,148 floats) for structure i: shape of one element (4716, 3)
    5. validation_status[i] is True if structure i passed the validation step
    6. createsims_status[i] is equals
        * to 2 if the status of the createsims running structure i is unknown (not failed nor successful)
        * to 1 if the status of createsims is success
        * to 0 if the createsims failed

    Contents of feedback frame file:
    The keys "structures" contains a dictionnary of structures and their positions in LS for each frames
        >>> import nump as np
        >>> data = np.load("db-feedback-frames.npz", allow_pickle=True)
        >>> structures = data["structures"].tolist()  # The function .tolist() is VERY important here
        >>> structures.keys()
        dict_keys(['structure_iter00_000000000803', 'structure_iter00_000000000032', 'structure_iter00_000000000017'])
        >>> structures["structure_iter00_000000000803"][0].shape #We select the first frame of structure 803
        (1, 64)

    This file is updated in real times with new structures. Please open it in read-only

    """

    def __init__(self,
            structures_db: str,
            feedback_file: str = None,
            feedback_frame: str = None,
            epsilon = 10,
            variance_threshold = 1e-8,
            proposal_limit = 1e6,
            log_diagnostics = False,
            stopping_buffer_size = 10, # set to 1 for very aggressive stopping 
            minimum_sim_length = 200. # 100 ns for a 0.5 ns framerate.
        ) -> None:
        """
        Setup a new LS sampler.
        """
        # This contains all the sructure pre-generated
        # the fields are and shapes are (for 1000 points):
        #   model_v2ucg (1,)
        #   model_ucg2cg (1,)
        #   model_ucg_encoder (1,)
        #   v (1000, 4)
        #   pred_ucg (1000, 250, 3)
        #   pred_cg (1000, 4716, 3)
        #   ls_coords (1000, 1, 20)
        #   selected (1000,)
        # ls_coords are the UCG encodings of pred_cg
        self.structure_db = structures_db
        # this is a NPZ file that contains informations about the status of createsims runs
        # contains a field named createsims_status where 0 = Fail, 1 = Success and 2 = Unknown (createsim is still running)
        self.feedback_file = feedback_file  # expects a string which is a path
        self.feedback_frame_path = feedback_frame  # Path to DB that contains frames from cganalysis encoded into LS coordinates
        self.feedback_db = {}  # place to store the dictionary read from feedback file
        self.epsilon = epsilon  # this is the neighborhood distance parameter
        self.variance_threshold = variance_threshold  # this parameter determines when to stop sampling this local region.
        self.proposal_limit = proposal_limit
        self.log_diagnostics = log_diagnostics

        self.feedback_db["v_coords"] = np.array([]) # Log v
        self.feedback_db["v2ucg_coords"] = np.array([])
        self.feedback_db["ucg2cg_coords"] = np.array([])
        self.feedback_db["ls_coords"] = np.array([])

        self.feedback_db["createsims_status"] = np.array([])
        self.feedback_db["validation_status"] = np.array([])

        if not self.structure_db_exists():
            raise ValueError(f"The DB {self.structure_db} does not exist")

        self.stopping_log = {}
        self.stopping_flags = []  # a list of currently running simulations that need to be stopped
        self.stopping_buffer_size = stopping_buffer_size
        self.minimum_sim_length = minimum_sim_length  # how long each simulation trajector must be before checking stopping 
        self.dt = 0.5 # rate at which a frame is produced in ns (500 ps for MuMMI) @TODO not hardcode that in the code

        if self.log_diagnostics:
            self.variance_criterion = []
            self.bernouli_prob = []
            self.sample_size = []

        mummi_core.init() # Required to use Naming function (Naming.dir_root)
        # IOinterface to stop simulations
        self.iointerface = mummi_ras.get_io("mummi")

    def structure_db_exists(self) -> bool:
        """
        Return true if feedback file is available
        """
        return self.structure_db is not None and os.path.isfile(self.structure_db)

    def feedback_exists(self) -> bool:
        """
        Return true if feedback file is available
        """
        return self.feedback_file is not None and os.path.isfile(self.feedback_file)

    def feedback_frame_exists(self) -> bool:
        """
        Return true if feedback files for frame is available
        """
        return self.feedback_frame_path is not None and os.path.isdir(self.feedback_frame_path)


    def _random_selection(self, k_samples: int, array: npt.NDArray) -> npt.NDArray:
        """
        Select randomly k_samples from a numpy array
        """
        return array[np.random.choice(len(array), size = k_samples, replace = False)]

    def _count_success(self, neighbors) -> int:
        """
        Return the number of successful points based on feedback.

            If validation_status = True
                if createsims_status = 1 => this createsims is successful
                if createsims_status = 0 => this createsims is a failure
                if createsims_status = 2 => this createsims is possibly running or is queued
            If validation_status = False
                Then no matter what createsims_status = 2

        """
        if self.feedback_exists(): 
            nsuccess = np.sum(self.feedback_db["createsims_status"][neighbors] == 1)
            nfailures = np.sum(self.feedback_db["createsims_status"][neighbors] == 0)
            ntrials = nsuccess + nfailures
        else:
            nsuccess = 0
            ntrials = 0

        return nsuccess, ntrials

    def get_ls_neighborhood(self, lsp) -> list:
        """
        Return an array of neighboring LS points.
        """
        if self.feedback_exists() and len(self.feedback_db["ls_coords"]) > 0:
            dvec = np.linalg.norm(self.feedback_db["ls_coords"] - lsp, axis=1)

            nbrs_indx = np.where(dvec <= self.epsilon)[0]
        else:
            return []

        return nbrs_indx

    def accept_point(self, lspoint) -> bool:
        """Returns a boolean whether to accpet or reject a new LS point based on validation information.
           It is better to implement this as an attribute that is selected only once, ow.,
           the if below will be executed every time a point is inquired

        if hasattr(self.interpolator, 'accept_point'):
           return getattr(self.interpolator, 'accept_point')

        Now we use a bernoulli model to assess acceptance criterion for drawing new LS points.
        """
        nbrs = self.get_ls_neighborhood(lspoint)  # this returns a boolean ndarray encoding indexes for feedback_db data structures

        s,n = self._count_success(nbrs)

        if self.log_diagnostics:
            self.sample_size.append(n)

        # if there are no neighbors, accept
        if n == 0:
            return True

        # infer the bernoulli success rate
        _success_probability = float(s) / n

        assert _success_probability <= 1.0

        if self.log_diagnostics:
            self.bernouli_prob.append(_success_probability)

        # use the 95% confidence interval as the acceptance probability for points with only failed neighbors
        if _success_probability == 0.0:
            # rule of threes estimator for confidence interval
            if random.random() <= 3.0 / n:
                return True

        # expected variance reduction: 
        new_var = _success_probability * (1 -_success_probability) / (n*(1+n))

        if self.log_diagnostics:
            self.variance_criterion.append(new_var)

        if new_var <= self.variance_threshold:
            return False

        if random.random() <= _success_probability:
            return True

        return False

    def _feedback_selection(self, k_samples: int, array: npt.NDArray, lsparray: npt.NDArray) -> npt.NDArray:
        """
        Select k_samples from a numpy array based on some feedback information

        TODO: implementation left to @Jeremy
        """

        # pull uniform random 
        _unif_draw = self._random_selection(k_samples, array)

        _accepted = np.array([self.accept_point(lsp) for lsp in lsparray[_unif_draw]],
                             dtype=bool)

        return _unif_draw[_accepted] 

    def select_points(self, k_samples: int, array: npt.NDArray, lsparray: npt.NDArray = None, method: str = "random") -> npt.NDArray:
        """
        Select k_samples from a numpy array
        """
        if method == "random" or lsparray is None:
            return self._random_selection(k_samples, array)
        elif method == "feedback":
            return self._feedback_selection(k_samples, array, lsparray)
        else:
            raise ValueError(f"Invalid method {method} (Accepted are 'random' or 'feedback')")

    def select_samples_from_npz(self, k_samples: int):
        """
        k_samples: int = Number of samples to select

        This method opens the pre-generated DB of structures and select k_samples of it
        The selection method can be random or based on feedback (the method select_points
        is responsible of performing the selection)

        At the end, we update the database of pre-generated structures to update the one
        we selected.
        """

        try:
            load_time = time.time()
            structure_database = np.load(self.structure_db, allow_pickle = True)
            load_time = time.time() - load_time
            status = structure_database["selected"]
            LOGGER.info(f"Loaded {len(status)} structures in {load_time} seconds")
        except Exception as e:
            LOGGER.error(f"While opening {self.structure_db}: {e}")
            raise RuntimeError(e)

        LOGGER.debug(f"We use pre-computed structures from {self.structure_db}")

        if "selected" not in structure_database.files:
            raise RuntimeError(f"Missing keys (selected, status) in {self.structure_db}")

        # Prevent to select more than we
        selected = structure_database["selected"]
        selected_struct = np.where(selected == False)[0]
        vids = np.where(selected == False)[0]

        LOGGER.info(f"Total structures {len(selected)} / structures available for selection : {len(selected_struct)}")

        if len(vids) == 0:
            LOGGER.info(f"No structures found. They have all been selected")
            return np.array([]), np.array([]), np.array([]), np.array([])

        if k_samples > len(vids):
            LOGGER.info(f"Number of samples requested {k_samples} -> {len(vids)}")
            k_samples = len(vids)

        # lsparray = None        
        lsparray = np.squeeze(structure_database["ls_coords"])

        # We select randomly here
        selected_ids = self.select_points(
            k_samples = k_samples,
            array = vids,
            lsparray = lsparray,
            method = "feedback"
        )
        # TODO: @Jeremy here is the API for the feedback
        # selected_ids = self.select_points(
        #     k_samples = k_samples,
        #     array = vids,
        #     method = "feedback"
        # )
        
        np.random.shuffle(selected_ids)

        selected[selected_ids] = np.full(selected_ids.size, True)

        if len(selected_ids) == 0:
            return np.array([]), np.array([]), np.array([]), np.array([])

        v_coords = structure_database["v"][selected_ids]
        v2ucg_coords = structure_database["pred_ucg"][selected_ids]
        pred_cg = structure_database["pred_cg"][selected_ids]

        # We update the pre-generated DB (npz)
        arrays = {key: structure_database[key] for key in structure_database.files}
        for key in structure_database.files:
            arrays[key] = structure_database[key]
        arrays["selected"] = selected
        # File we will be writing (just in case to not corrupt the file in case of interruption)
        timestr = time.strftime("%Y%m%d-%H%M%S")
        db_path = os.path.splitext(self.structure_db)
        input_samples_tmp = f"{db_path[0]}-{timestr}{db_path[1]}"
        # Update NPZ data base
        np.savez_compressed(input_samples_tmp, **arrays)

        # Now we replace the files (atomic all)
        os.replace(input_samples_tmp, self.structure_db)

        return selected_ids, v_coords, v2ucg_coords, pred_cg

    def get_diagnostics(self):
        """Return the diagnostics information."""

        if self.log_diagnostics:
            return self.variance_criterion, self.bernouli_prob, self.sample_size

    def _update_feedback_db_from_file(self) -> None:
        """Re-read the feedback file into memory.
        """
        LOGGER.info(f"Updating feedback database from {self.feedback_frame_path}")

        if self.feedback_exists():
            LOGGER.info(f"Reading feedback sampling db.")
            with np.load(self.feedback_file, allow_pickle=True) as dbf:
                # update the list of validation successes
                self.feedback_db["ls_coords"] = np.squeeze(dbf["ls_coords"])
                self.feedback_db["createsims_status"] = dbf["createsims_status"]
                self.feedback_db["validation_status"] = dbf["validation_status"].astype(bool)

                # if validation fails, we also mark a failed createsims, ignore cs_status = 2 
                self.feedback_db["createsims_status"][np.logical_not(self.feedback_db["validation_status"])] = 0

        if self.feedback_frame_exists():

            # local buffer to accumulate fb information 
            _fb_ls_traj = {}

            # list the split format files
            fb_directory = os.path.dirname(self.feedback_frame_path)
            fbFiles = [dbf for dbf in os.listdir(fb_directory) if re.match("db-feedback-frames_\d{12}.npz", dbf)]
            # just make sure it's ordered right
            fbFilesNames = [f"db-feedback-frames_{fbi:012}.npz" for fbi in range(0, len(fbFiles)-1)]
            fbFilesNames.append("db-feedback-frames.npz")

            for fbn in fbFilesNames:
                LOGGER.info(f"Processing {fbn}.")
                with np.load(os.path.join(fb_directory, fbn), allow_pickle=True) as dbf:
                    # update the interpolator database
                    fb_trajectories = dbf['structures'].tolist()

                # this is a bit goofy... might need to change feedback file formatting later.
                for struct_name in fb_trajectories.keys():

                    #feedback_ls_trj = np.zeros((len(fb_trajectories[struct_name]), fb_trajectories[struct_name][0].shape[1]))
                    #for fid, frame in enumerate(fb_trajectories[struct_name]):
                    #    feedback_ls_trj[fid,:] = frame[0,:]

                    feedback_ls_trj = np.array(fb_trajectories[struct_name])
                    nframes = feedback_ls_trj.shape[0]
                    ls_size = feedback_ls_trj.shape[2]

                    feedback_ls_trj.reshape(nframes, ls_size)

                    if struct_name in _fb_ls_traj.keys():
                        _fb_ls_traj[struct_name] = np.concatenate((_fb_ls_traj[struct_name],
                                                                  feedback_ls_trj.reshape(nframes, ls_size)),
                                                                  axis=0)
                    else:
                        _fb_ls_traj[struct_name] = feedback_ls_trj.reshape(nframes,
                                                                           ls_size)

            # now we need ot update our stopping criteria list and send a stop signal.
            LOGGER.info(f"Checking for stopping criterion.")
            success_flag, fail_flag, stop_flag = Naming.status_flags("cg")
            sims_cg = Naming.dir_root("all-cg")

            for struct_name in _fb_ls_traj.keys():

                sim_dir = os.path.join(sims_cg, struct_name)

                if not os.path.isdir(sim_dir):
                    continue

                # skipp if there are flags and take out of log 
                if self.iointerface.test_signal(sim_dir, success_flag):
                    if struct_name in self.stopping_log.keys():
                        del self.stopping_log[struct_name]
                    continue
                elif self.iointerface.test_signal(sim_dir, fail_flag):
                    if struct_name in self.stopping_log.keys():
                        del self.stopping_log[struct_name]
                    continue
                elif self.iointerface.test_signal(sim_dir, stop_flag):
                    if struct_name in self.stopping_log.keys():
                        del self.stopping_log[struct_name]
                    continue

                nframes = _fb_ls_traj[struct_name].shape[0]

                if struct_name in self.stopping_log.keys():
                    slen = self.stopping_log[struct_name]
                    # if there isn't anything new, move on 
                    if nframes <= slen:
                        continue
                    # update record of what we've seen 
                    elif nframes > slen:
                        self.stopping_log[struct_name] = nframes

                else:
                    self.stopping_log[struct_name] = nframes

                # skip if this one hasn't run long enough
                if len(_fb_ls_traj[struct_name]) * self.dt < self.minimum_sim_length:
                    continue

                stopping = self.check_stopping_criterion(_fb_ls_traj[struct_name][-self.stopping_buffer_size:]) # should take last 10  

                if stopping:
                    if struct_name not in self.stopping_flags:
                        self.stopping_flags.append(struct_name)
                        LOGGER.info(f"Stopping {struct_name} based on feedback")

            # now send stopping signal.
            LOGGER.info(f"Issuing {len(self.stopping_flags)} stop signals out of {len(_fb_ls_traj.keys())} simulations.")
            self._send_stopping_signals()
        else:
            LOGGER.info(f"Skipping updating {self.feedback_frame_path}. DB is empty")


    # Here we rename these to be coherent between C2 and C3 and have a common interface
    def get_new_ls_points(self, npoints: int):
        """
        The method returns npoints from the prebuilt database
        """
        return self.select_samples_from_npz(npoints)

    def check_stopping_criterion(self, trj):
        """
        Return true or false whether this simulation should be stopped.
        """

        acceptance_buffer = np.array([self.accept_point(lsp) for lsp in trj], dtype=bool)
        return acceptance_buffer.sum() == 0


    def _send_stopping_signals(self):
        """
        Read the buffer of simulations to stop and update the buffer and running status.
        """

        for struct_name in self.stopping_flags:
            self._stop_simulation(struct_name)

        return

    def _stop_simulation(self, sim_name: str):
        """
        Stop a simulation bu writing a flag cg_stop in its directory

        sim_name: str simulation name such as "structure_iter00_000000002709"
        """
        success_flag, fail_flag, stop_flag = Naming.status_flags("cg")
        sims_cg = Naming.dir_root("all-cg")
        sim_dir = os.path.join(sims_cg, sim_name)

        if not os.path.isdir(sim_dir):
            LOGGER.error(f"{sim_dir} is not a valid directory")
            return
        if self.iointerface.test_signal(sim_dir, success_flag) or self.iointerface.test_signal(sim_dir, fail_flag):
            LOGGER.info(f"{sim_dir} has already a success or fail flag")
            return
        if not self.iointerface.test_signal(sim_dir, stop_flag):
            self.iointerface.send_signal(sim_dir, stop_flag)
            LOGGER.info(f"Setting stop flag for {sim_dir}")
