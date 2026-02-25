"""
Sampler class for drawing new points from the latent space.
"""

import numpy as np
from typing import Type
import collections
import random
import os
import re
from typing import List, Union

import mummi_core
from mummi_core.utils import Naming
import mummi_ras
from mummi_ras.ml.ls_point import LSPoint

#from .ot_sampler import OTSampler

from logging import getLogger
LOGGER = getLogger(__name__)

class LSSampler_C2:
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

    The LSSampler module relies on another class called an LSInterpolator which encapsulates
    the model used to interpolated between the states A and B. The LSSampler class will propose new
    points by combining the output of LSInterpolator and the validation information it gets from the
    database to define the proposal distribution.
    This allows for some modularity with which the LSSampler builds a proposal distribution.

    NTS: ok, one issue is how we want the sampling to draw points.. if the proposal distribution is
    defined as a weighted draw from teh set of interpolated points, then we might need a tighter
    integration between the interpolator and the sampler classes?? alternatively, we may just have
    the sampler draw points uniformly from the interpolator and accpet / reject points to yield th weighted draw...

    Note: The sample maintains a buffer of new LS points in the new_points_buffer which is a double ended queue structure.
    The "left" side is newest points and the "right" side is the oldest. I.e. the oldest points are pulled off first when asked for points.

    Contents of feedback file:
    1. structure_names[i] contains the name of the structure i
    2. ls_coords[i] contains the coords (32 floats) for structure i
    3. validation_status[i] is True if structure i passed the validation step
    4. createsims_status[i] is equals
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
            states: List[np.ndarray],
            interpolator,
            feedback_file:str = None,
            feedback_frame: str = None,
            epsilon = 1.5,
            variance_threshold = 1e-8,
            proposal_limit = 1e6,
            log_diagnostics = False,
            min_buffer_size = 500,
            stopping_buffer_size = 10,
            minimum_sim_length = 200. # 100 ns for a 0.5 ns framerate.
        ) -> None:
        """Setup a new LS sampler.
        """

        #self.samples = []  # place to store a working list of previously sampled points and validation information
        self.new_points_buffer = collections.deque() # place to store a buffer of new LS points to propose (FIFO queue).
        self.states = states  # the set of points labeled as A
        self.min_buffer_size = min_buffer_size
        # this is a NPZ file that contains informations about the status of createsims runs
        # contains a field named createsims_status where 0 = Fail, 1 = Success and 2 = Unknown (createsim is still running)
        self.feedback_file = feedback_file  # expects a string which is a path
        self.feedback_frame_path = feedback_frame  # Path to DB that contains frames from cganalysis encoded into LS coordinates
        self.feedback_db = {}  # place to store the dictionary read from feedback file
        self.epsilon = epsilon  # this is the neighborhood distance parameter
        self.variance_threshold = variance_threshold  # this parameter determines when to stop sampling this local region.
        self.proposal_limit = proposal_limit
        self.log_diagnostics = log_diagnostics
        self.feedback_db["ls_coords"] = np.array([])
        self.feedback_db["createsims_status"] = np.array([])
        self.feedback_db["validation_status"] = np.array([])
        self.stopping_log = {}
        self.stopping_flags = []  # a list of currently running simulations that need to be stopped
        self.stopping_buffer_size = stopping_buffer_size
        self.minimum_sim_length = minimum_sim_length  # how long each simulation trajector must be before checking stopping 
        self.dt = 0.5 # rate at which a frame is produced in ns (500 ps for MuMMI) @TODO not hardcode that in the code

        if self.log_diagnostics:
            self.variance_criterion = []
            self.bernouli_prob = []
            self.sample_size = []

        self.lsp_index = 0  # this tracks how many latent space points have been created. 

        # initialize the point index
        self._initialize_lsp_index()

        # insert interpolator initialization here.
        self.interpolator = interpolator

        mummi_core.init() # Required to use Naming function (Naming.dir_root)
        # IOinterface to stop simulations
        self.iointerface = mummi_ras.get_io("mummi")

        self._extend_buffer(self.min_buffer_size)

        return

    def _initialize_lsp_index(self):
        """
        Determine initial lsp_index from source and target sets.
        """

        index_set = []
        _max = 0
        for si in self.states:
            index_set.append(np.array([p.get_id() for p in si]))
            if index_set[-1].max() > _max:
                _max = index_set[-1].max()

        self.lsp_index = int(_max + 1)

        return


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

    def get_diagnostics(self):
        """Return the diagnostics information."""

        if self.log_diagnostics:
            return self.variance_criterion, self.bernouli_prob, self.sample_size

        return


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

        nsuccess = np.sum(self.feedback_db["createsims_status"][neighbors] == 1)
        nfailures = np.sum(self.feedback_db["createsims_status"][neighbors] == 0)
        ntrials = nsuccess + nfailures

        return nsuccess, ntrials


    def get_ls_neighborhood(self, lsp) -> list:
        """
        Return an array of neighboring LS points.
        """
        # No feedback at that time
        if self.feedback_db["ls_coords"].shape[0] == 0:
            return []
        dvec = np.linalg.norm(self.feedback_db["ls_coords"] - lsp, axis=1)

        nbrs_indx = np.where(dvec <= self.epsilon)[0]

        return nbrs_indx

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


    def _update_feedback_db_from_file(self) -> None:
        """Re-read the feedback file into memory.
        """
        LOGGER.info(f"Updating feedback database from {self.feedback_frame_path}")

        if self.feedback_exists():
            LOGGER.info(f"Reading feedback sampling db.")
            with np.load(self.feedback_file, allow_pickle=True) as dbf:
                # update the list of validation successes
                self.feedback_db["ls_coords"] = dbf["ls_coords"]
                self.feedback_db["createsims_status"] = dbf["createsims_status"]
                self.feedback_db["validation_status"] = dbf["validation_status"].astype(bool)

                # if validation fails, we also mark a failed createsims, ignore cs_status = 2 
                self.feedback_db["createsims_status"][np.logical_not(self.feedback_db["validation_status"])] = 0

        if self.feedback_frame_exists():
            # rebuild from scratch.  

            LOGGER.info(f"Resetting edge database.")
            self.interpolator.reset_feedback_edges()

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
            # finally pass the database buffer to the interpolator
            for struct_name in _fb_ls_traj.keys():
                self.interpolator.add_new_fb_edge(struct_name,
                                                  _fb_ls_traj[struct_name])

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

    def check_stopping_criterion(self, trj):
        """
        Return true or false whether this simulation should be stopped.
        """

        acceptance_buffer = np.array([self.accept_point(lsp) for lsp in trj], dtype=bool)

        if acceptance_buffer.sum() == 0:
            stop = True
        else:
            stop = False

        return stop

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

    def _extend_buffer(self, npoints: int) -> None:
        """add k new points to the list of new LS points
        """

        if self.feedback_exists():
            # read feedback file
            self._update_feedback_db_from_file()

        accepted_points = []
        counter = 0
        while len(accepted_points) < npoints:

            # this should return a list of nparrays 
            nlsp, lambdas = self.interpolator.draw_k_points(self.min_buffer_size + npoints)

            # filter or weight this proposed set based on the validation information
            accepted_points.extend([pt for pt in zip(nlsp, lambdas) if self.accept_point(pt[0])])

            if counter > self.proposal_limit:
                # NEED TO PASS A MESSGAE OR OTHER THING HERE 
                break

        # init new LSP after accepted. 
        new_buffer_pts = [self.new_lsp_point(pt[0], None, False, initial_lambda=pt[1]) for pt in accepted_points[0:npoints]]

        self.new_points_buffer.extend(new_buffer_pts)

        return


    def new_lsp_point(self, crd, parent, valid, initial_lambda = None):
        """Factory for constructing a new LSPoint instance."""

        lsp = LSPoint(self.lsp_index,
                      parent=parent,
                      ls_coord=crd,
                      valid=valid)

        self.lsp_index += 1

        if initial_lambda is not None:
            lsp.initial_lamda = initial_lambda

        return lsp


    def get_new_ls_points(self, npoints: int):
        """ returns the first npoints and deletes them from the list self.new_ls_points """

        # check to see if we need more points
        if npoints > len(self.new_points_buffer):
            # put enough new points in plus the min_buffer_size
            self._extend_buffer(self.min_buffer_size + npoints)

        _npts = [self.new_points_buffer.popleft() for x in range(npoints)]

        return _npts

