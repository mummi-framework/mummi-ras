#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2023, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------
# @authors
#           Loic Pottier <pottier1@llnl.gov>
# ------------------------------------------------------------------------------

import os
import multiprocessing
import signal
import traceback
import numpy as np
import pathlib
from typing import List, Union, Any
import yaml
import fcntl
import time
import datetime
import tempfile
import pickle
import glob
from filelock import FileLock
from timeit import default_timer as timer

import mummi_core
from mummi_core.brokers import RPCServer
from mummi_core.workflow.flux_env import flux_uri

from mummi_ras import Naming

from mummi_ras.ml.models.autoencoders import FullDenseAutoencoder, FullMiniDenseAutoencoder, HierarchicalCGAutoencoder

from mummi_ras.ml.validators import CGValidator_C2
from mummi_ras.ml import ls_point as lsp
from mummi_ras.ml.samplers import ls_sampler, get_interpolator
from mummi_ras.ml.feedback_frames import FeedbackFrames

from logging import getLogger, Logger
LOGGER = getLogger(__name__)

def write_patches(outpath: str, iteration_id: int, new_positions: np.array) -> Union[List[str], List[Any]]:
    """
    From positions, write a patch as a npz file and 
    return its path.
    """
    iteration_id = "{:02d}".format(iteration_id)
    RPATH = os.path.join(outpath, f"iter{iteration_id}")
    if not os.path.isdir(RPATH):
        os.makedirs(RPATH, exist_ok=True)

    LOGGER.info(f"Write patches in {RPATH}")
    npatches = new_positions.shape[0]
    files = np.empty(npatches, dtype='object')
    positionsArray = []
    namesArray = []
    all_structure_names = []

    # We want to avoid erasing previously generated structures
    offset = _latest_patch_id_generated(RPATH)

    # Create an .npz for each new patch
    for sidx in range(npatches):

        positions = new_positions[sidx]
        structure_id = "{:012d}".format(offset+sidx)

        structure_name = f"structure_iter{iteration_id}_{structure_id}"
        all_structure_names.append(structure_name)
        outfile_positions = os.path.join(RPATH, f"{structure_name}.npz")
        files[sidx] = f"{structure_name}.npz"
        LOGGER.info(f"Writing {outfile_positions}")
        np.savez_compressed(outfile_positions, data = positions)

        positionsArray.append(positions)
        namesArray.append(structure_name)

    return namesArray, positionsArray

def _latest_patch_id_generated(RPATH: str) -> int:
    """
    Find in the directory containing potentially already some structures
    like structure_iter00_X.gro. This function finds the highest X so we can
    write the patch X+1 safely without erasing any structures. 
    """
    nsamples = len(glob.glob(RPATH + '/*.npz'))
    nsamples = nsamples - 2 # the - 2 is to remove master_iterX.npz and table_all_predictions.npz
    return max(nsamples, 0) # to avoid returning a negative number if RPATH was empty


def create_sampling_db(database: str, model_name: str):
    """
    Create an empty sampling feeback DB.
    """
    dirname = os.path.dirname(database)
    ts = datetime.datetime.now().strftime('%d.%m.%Y-%H:%M:%S')
    if dirname != '':
        os.makedirs(dirname, exist_ok=True)
    data = {
        "model_name": model_name,
        "created_at": ts,
        "updated_at": ts,
        "structure_names": np.empty(0, dtype=str),
        "ls_coords": np.empty(0, dtype=float),
        "lambda_values": np.empty(0, dtype=bool),
        "validation_status": np.empty(0, dtype=bool),
        "createsims_status": np.empty(0, dtype=int)
    }
    np.savez_compressed(database, **data)
    LOGGER.info(f"Created new DB at {database}")

def checking_sampling_db(database: str, model_name: str) -> bool:
    """
    Update sampling feeback DB with new samples.
    """
    if not os.path.isfile(database):
        return True
    try:
        sampling_db = np.load(database, allow_pickle=True)
    except Exception as e:
        LOGGER.error(f"{database} {e}")
        return False
    prev_model_name = sampling_db['model_name']
    if model_name != prev_model_name:
        LOGGER.error(f"Old sampling DB cannot be used with different ML model.")
        LOGGER.error(f"This DB {database} has been created with ML model {prev_model_name}")
        LOGGER.error(f"You are currently used Latent Space {model_name}")
        return False
    return True

def update_sampling_db(
    database: str,
    num_new_sample: int,
    model_name: str,
    structure_names: List[str],
    ls_coords: List[List[float]] = [],
    lambda_values: List[float] = [],
    validation_status: List[bool] = [],
    new_createsims_status: List[bool] = []
    ):
    """
    Update sampling feeback DB with new samples.
    """

    # File we will be writing (just in case to not corrupt the file in case of interruption)
    timestr = time.strftime("%Y%m%d-%H%M%S")
    db_path = os.path.splitext(database)
    database_tmp = f"{db_path[0]}-{timestr}{db_path[1]}"

    if not os.path.isfile(database):
        LOGGER.info(f"Sampling DB ({database}) did not exist, creating it.")
        create_sampling_db(database, model_name)

    try:
        with np.load(database, allow_pickle=True) as sampling_db:
            prev_model_name = sampling_db['model_name']
            created_at = sampling_db['created_at']
            prev_struct = sampling_db['structure_names']
            prev_coords = sampling_db["ls_coords"]
            prev_lambda = sampling_db["lambda_values"]
            prev_validation = sampling_db["validation_status"]
            prev_createsims = sampling_db["createsims_status"]
    except Exception as e:
        LOGGER.error(f"{database} {e}")
        return

    new_records = len(structure_names)
    ts = datetime.datetime.now().strftime('%d.%m.%Y-%H:%M:%S')
    if model_name != prev_model_name:
        LOGGER.error(f"Old sampling DB cannot be used with different ML model.")
        LOGGER.error(f"This DB {database} has been created for LS coordinates of LS:")
        LOGGER.error(f" - {prev_model_name}")
        LOGGER.error(f"You are currently used Latent Space {model_name}")
        LOGGER.error(f"Feedback is deactivated.")
        return

    data = {
        "model_name": model_name,
        "created_at": created_at,
        "updated_at": ts
    }

    # Update sampling DB for feedback
    data["structure_names"] = np.append(prev_struct, structure_names)
    if ls_coords != []:
        if len(prev_coords) > 0:
            data["ls_coords"] = np.append(prev_coords, ls_coords, axis=0)
        else:
            data["ls_coords"] = ls_coords
    else:
        data["ls_coords"] = prev_coords
    if lambda_values != []:
        data["lambda_values"] = np.append(prev_lambda, lambda_values)
    else:
        data["lambda_values"] = prev_lambda
    if validation_status != []:
        data["validation_status"] = np.append(prev_validation, validation_status)
    else:
        data["validation_status"] = prev_validation
    # 0 means False, 1 True and 2 is unknown state for createsims
    # As no createsims for these structures are running, they are all unknowns
    if new_createsims_status != []:
        data["createsims_status"] = new_createsims_status
    else:
        data["createsims_status"] = np.append(prev_createsims, np.full(num_new_sample, 2, dtype=int))

    np.savez_compressed(database_tmp, **data)
    # We make sure the new file is not corrupted somehow
    try:
        with np.load(database_tmp, allow_pickle=True) as test:
            LOGGER.debug(f"{database_tmp} is valid {test.files}")
    except Exception as e:
        LOGGER.warning(f"{database_tmp} seems to be corrupted. We keep the old {database} intact")
        return

    os.replace(database_tmp, database)

    LOGGER.debug(f"Added new structures: {structure_names}")
    LOGGER.info(f"Upgraded DB at {database} with {new_records} new records")

def generate_new_samples(**config: dict) -> List[str]:
    """
    Sample k_sample from the latent space, generate the k_samples corresponding
    structures and validate them. Return the structures that passed the validation.
    """
    process_sampler = config["process_generator"]
    # We are waiting for the sampler to be created in a dedicated thread
    if process_sampler.is_alive():
        LOGGER.info(f"Waiting on sampler ({process_sampler}) to be created...")
        process_sampler.join()

    iteration_id = config["iteration_id"]
    k_samples = max(int(config["body"].get("k_samples")), 0)
    strict = bool(config["body"].get("strict"))
    # Get back objects
    sampler = config["manager"]["sampler"]
    generator = config["obj_generator"]
    validator = config["obj_validator"]
    # --------------------------------------------------------------------------
    nMaxSelectedPatchBuffer: int = int(config["nMaxSelectedPatchBuffer"])
    # This factor defines the factor of extra structure we will generate
    # Example: if equals to 2, we will select up to 2 * nMaxSelectedPatchBuffer structure at most
    # It is useful in case a lot of createsim are failing (which the ML server does not know about)
    factor_extra_structures = int(config["sampler"]["factor_extra_structures"])

    # # Setting feedback DB for sampler
    feedback_db_path = config["sampling_db"]
    # We do not need to select more than what is defined in wfmanager.yaml
    num_validated = validator.num_validated(iteration_id)
    if num_validated >= factor_extra_structures * nMaxSelectedPatchBuffer:
        LOGGER.info(f"=> {num_validated} / {factor_extra_structures * nMaxSelectedPatchBuffer} have already been selected")
        return []
    else:
        LOGGER.info(f"=> {num_validated} / {factor_extra_structures * nMaxSelectedPatchBuffer}. We can sample more")

    new_structure_names = []
    new_ls_coords = []
    new_lambda = []
    new_validation_status = []

    total_valid_files = []
    loop_iter = 0
    num_new_sample = 0
    while len(total_valid_files) < k_samples:
        sample_start = time.time()
        ls_coords = sampler.get_new_ls_points(k_samples)
        sample_end = time.time() - sample_start
        LOGGER.info(f"Sampled {k_samples} in {sample_end} seconds")

        for pts in ls_coords:
            new_ls_coords.append(pts.get_coordinates())
            new_lambda.append(pts.get_lamda())
        num_new_sample += len(ls_coords)

        new_positions = generator.decode(ls_coords)
        names_array, positions_array = write_patches(
            outpath = config["generator"]["outpath"],
            iteration_id = iteration_id,
            new_positions = new_positions
        )
        new_structure_names.append(names_array)
        LOGGER.debug(f"generated structures done. new_positions = {new_positions.shape}")

        return_array = validator.validateArrayThreaded(iteration_id, names_array, positions_array)
        # return_array = validator.validateArray(iteration_id, names_array, positions_array)
        if len(return_array) > 0:
            new_validation_status.append(return_array[:,0])

        _, valid_files = validator.write_validation_info(
            iteration_id = iteration_id,
            return_array = return_array,
            all_structure_names = names_array
        )
        valid_files = [os.path.join(validator.current_rpath, f.split('.gro')[0]) for f in valid_files]
        LOGGER.debug(f"ITER = {loop_iter} => valid structures={[os.path.join(validator.current_rpath, f) for f in valid_files]}")
        total_valid_files += list(valid_files)

        with FileLock(config["lock_sampling_db"]) as _:
            # Update sampling DB for feedback
            update_sampling_db(
                database = feedback_db_path,
                model_name = config["encoder_path"],
                num_new_sample = len(ls_coords),
                structure_names = new_structure_names,
                ls_coords = new_ls_coords,
                lambda_values = new_lambda,
                validation_status = new_validation_status
            )

        loop_iter += 1
        num_validated = validator.num_validated(iteration_id)
        if num_validated >= factor_extra_structures * nMaxSelectedPatchBuffer:
            LOGGER.info(f"=> {num_validated} / {factor_extra_structures * nMaxSelectedPatchBuffer} have been selected, returning []")
            break
        if not strict:
            break

    LOGGER.info(f"[{loop_iter}] valid_files({len(total_valid_files)})={total_valid_files}")
    return total_valid_files

def create_sets(training_data_files: List[str]) -> List[List]:
    sets = {}
    assert(len(list(training_data_files)) != 0)

    # detect number of sets
    for fn in training_data_files:
        data = np.load(fn)
        state_labels = sorted(np.unique(data['labels']))
        for label in state_labels:
            bary = data['labels'] == label
            if label in sets:
                sets[label] = np.concatenate((sets[label], data['x_test_encoded'][bary]))
            else:
                sets[label] = data['x_test_encoded'][bary]

    counter = 0
    lsp_sets = []
    for label in sets:        
        lsp_sets.append([lsp.LSPoint(counter + sets[label].shape[0],
                                     parent=str(label),
                                     ls_coord=sets[label][i], 
                                     valid=True) for i in range(sets[label].shape[0])])
        counter += sets[label].shape[0]

    return lsp_sets

class MLServerRPC(object):
    """
    ML server implementation with RPC support 
    """
    def __init__(self, config: dict, lockfile: str, logger: Logger) -> None:
        self.config = config
        self.rpc_server: RPCServer = None
        self.logger = logger
        self.hostname = mummi_core.get_hostname(contract_hostname=False)

        self._exit = multiprocessing.Event()  # Trigger the daemon to exit.
        self._error = multiprocessing.Event() # Trigger the daemon to exit.
        self._setup = multiprocessing.Event() # Wait until the daemon is setup.
        self._manager = multiprocessing.Manager()
        self.feedback = None # Could be process if we choose to activate feedback
        self.lockfile = lockfile or tempfile.mkstemp(prefix="mummi-sampling-", suffix=".lock", dir=None)[1]
        self.waittime = 180 # time between two checks for createsims status (in seconds)

        self.mini_mummi = False

    def is_exit(self):
        return self._exit.is_set()

    def is_error(self):
        return self._error.is_set()

    def exit(self):
        self._exit.set()

    def error(self):
        self._error.set()

    def signal_wrapper(self, name, pid):
        def handler(signum, frame):
            self.logger.warning(f"Received (SIGNUM={signum}) for {name}[pid={pid}], stopping the process...")
            # self.validator.stop()
            self.rpc_server.stop()
            self.rpc_server.purge()
            self.rpc_server.close()
            if self.feedback:
                try:
                    self.feedback.join() # We stop the process updating the feedback file
                    self.feedback_frames.join()
                except AssertionError as e:
                    pass
                try:
                    self.logger.info(f'Removing {self.lockfile}')
                    os.remove(self.lockfile)
                except FileNotFoundError as _:
                    self.logger.warning(f'Lock file not found {self.lockfile}')
                    pass
            self.exit()
            exit(1)
        return handler

    def read_specs(self):
        try:
            # ------------------------------------------------------------------
            # identify scheduler interface
            self.logger.info(f'ML Server launched on host: {self.hostname}')
            self.flux = flux_uri(override_uri_from_file=False)
            self.logger.info(f"  flux  uri: [{self.flux}]")
            flux_uri_local = os.environ.get("FLUX_URI", None)
            self.logger.info(f"  local uri: [{flux_uri_local}]")

            if self.flux is None:
                self.logger.error("Failed to identify Flux URI. Exiting ...")
                # assert False, "Failed to identify flux uri"

            # --------------------------------------------------------------------------
            if Naming.MUMMI_MINI:
                wfmanager_yaml = os.path.join(Naming.MUMMI_SPECS, 'mini-mummi')
            else:
                wfmanager_yaml = os.path.join(Naming.MUMMI_SPECS, 'workflow')
            with open(os.path.join(wfmanager_yaml, 'wfmanager.yaml'), 'r') as data:
                wfmngr = yaml.load(data, Loader=yaml.FullLoader)
            try:
                self.nMaxSelectedPatchBuffer = int(wfmngr["wfmanager"]["config"].get("nMaxSelectedPatchBuffer"))
            except Exception as _:
                self.nMaxSelectedPatchBuffer = 0

            self.config["nMaxSelectedPatchBuffer"] = self.nMaxSelectedPatchBuffer

            try:
                self.iteration_id = int(wfmngr["wfmanager"]["config"].get("mlserver_round_id"))
            except Exception as _:
                self.iteration_id = 0            
            self.config["iteration_id"] = self.iteration_id
            self.iteration_path = "iter{:02d}".format(self.iteration_id)
            self.config["iteration_path"] = self.iteration_path

            self.sampling_db = os.path.join(self.config["sampler"]["outpath"], self.iteration_path)
            self.sampling_db = os.path.join(self.sampling_db, "feedback")
            os.makedirs(self.sampling_db, exist_ok=True)
            self.sampling_db = os.path.join(self.sampling_db, self.config["sampler"]["feedback"]["database"])

            self.feedbackframe_db = os.path.join(self.config["sampler"]["outpath"], self.iteration_path)
            self.feedbackframe_db = os.path.join(self.feedbackframe_db, "feedback")
            os.makedirs(self.feedbackframe_db, exist_ok=True)
            self.feedbackframe_db = os.path.join(self.feedbackframe_db, self.config["sampler"]["feedback"]["frame_database"])

            self.do_feedback = bool(self.config["sampler"]["feedback"]["do_feedback"])

            if "mini_mummi" in self.config["config"] or Naming.MUMMI_MINI:
                self.mini_mummi = bool(self.config["config"]["mini_mummi"])
            
            if self.mini_mummi:
                LOGGER.info("Run in mini MuMMI mode")

            self.config["sampling_db"] = self.sampling_db
            self.config["lock_sampling_db"] = self.lockfile
            self.sampler_interpolator = self.config["sampler"]["interpolator"]
            self.pickle_interpolator = self.config["sampler"].get("pre_computed")
            if self.pickle_interpolator and os.path.isfile(self.pickle_interpolator):
                self.logger.info(f"We will use a pre-computed interpolator: {self.pickle_interpolator}")

            self.encoder_name = self.config["encoder"]["model"]
            self.encoder_path = None

            # We take the first model that exists
            if isinstance(self.config["encoder"]["path"], list):
                for epath in self.config["encoder"]["path"]:
                    self.encoder_path = os.path.join(epath, self.encoder_name)
                    if os.path.exists(self.encoder_path):
                        self.logger.info(f"Found ML model {self.encoder_path}")
                        break
            else:
                self.encoder_path = os.path.join(self.config["encoder"]["path"], self.encoder_name)

            if os.path.exists(self.encoder_path):
                self.logger.info(f"Found ML model {self.encoder_path}")
            else:
                self.logger.error(f"Did not found any valid ML models")

            self.config["encoder_path"] = self.encoder_path
            self.credentials_path = os.path.join(
                self.config["workspace"]["path"],
                self.config["workspace"]["credentials"]
            )   
            self.certificate_path = os.path.join(
                self.config["workspace"]["path"], 
                self.config["workspace"]["certificate"]
            )

            usr = os.environ.get("USER", "mummiusr")
            self.queue = self.config["broker"]["queue"] + '_' + usr
            if usr == "mummiusr":
                self.logger.warning(f"Did not find current user: defaulted to {usr}")
            self.broker_interface = self.config["broker"]["interface"]

            self.logger.info(f"> Initializing MuMMI ML Server")
            self.logger.info(f"  > Server")
            self.logger.info(f"    > Interface                {self.broker_interface}")
            self.logger.info(f"    > Credentials              {self.credentials_path}")
            self.logger.info(f"    > Certificate              {self.certificate_path}")
            self.logger.info(f"    > Queue                    {self.queue}")
            self.logger.info(f"  > Sampler")
            self.logger.info(f"    > Interpolator             {self.sampler_interpolator}")
            self.logger.info(f"    > Run with feedback        {self.do_feedback}")
            self.logger.info(f"    > Createsims Feedback DB   {self.sampling_db}")
            self.logger.info(f"    > CG frames Feedback DB    {self.feedbackframe_db}")
            self.logger.info(f"  > Generator")
            self.logger.info(f"    > Encoder                  {self.encoder_path}")
            self.logger.info(f"  > Validator")
            self.logger.info(f"  > nMaxSelectedPatchBuffer    {self.nMaxSelectedPatchBuffer}")
        except Exception as e:
            traceback.print_exc()
            self.logger.error("ML Server failed during reading specs!")
            self.error()
            raise e

    def _setup_training_sets(self):
        training_dir = os.path.join(self.encoder_path, "training")
        # PC3: round 0: validation data are contained in many different files
        # training_data = pathlib.Path(training_dir).glob("validation_data_r*.npz")
        # PC3: round 1: If Konstantia generates new validation data for MuMMI which are contained in one file
        training_data = list(pathlib.Path(training_dir).glob("validation_data_all.npz"))
        LOGGER.info(f"Training data {training_dir} => {len(training_data)} files")
        states = create_sets(training_data)

        if self.mini_mummi:
            # Only two states in mini MuMMI
            states = states[:2]
        else:
            # states 4 and 5 are not useful for know (A''' and unknown state)
            states = states[:3]

        sub_sample_frac = float(self.config["sampler"]["sub_sample_frac"])
        assert 0.001 < sub_sample_frac <= 1

        # Needed for OTInterpolator because size of sets must be equals
        min_size = len(states[0])
        for s in states:
            if min_size > len(s):
                min_size = len(s)

        mask = np.random.uniform(size=(min_size,)) <= sub_sample_frac
        for i in range(len(states)):
            if len(states[i]) > min_size:
                states[i] = states[i][:min_size]
            # We sub-sample uniformly because POT is too slow for the whole states
            states[i] = list(np.array(states[i])[mask])

        return states

    def setup_sampler(self):
        if self.do_feedback:
            # Setting feedback DB for sampler
            feedback_db_path = self.config["sampling_db"]
            feedbackframe_db = self.feedbackframe_db
            if not os.path.isfile(feedback_db_path):
                create_sampling_db(feedback_db_path, self.encoder_path)

            if not checking_sampling_db(feedback_db_path, self.encoder_path):
                feedback_db_path = None
                feedbackframe_db = None
                LOGGER.warning(f"Feedback {feedback_db_path} is not valid for this ML model.")
                LOGGER.warning(f"All feedback is deactivated for this run.")
            else:
                LOGGER.info(f"Feedback DB for sampling located in {feedback_db_path}")
        else:
            feedback_db_path = None
            feedbackframe_db = None
            LOGGER.info(f"All feedback is deactivated for this run.")

        start = time.time()
        states = self._setup_training_sets()
        try:
            with open(self.pickle_interpolator, 'rb') as file:
                self.interpolator = pickle.load(file)
                end = time.time() - start
                LOGGER.info(f"Loaded pre-computed interpolator {self.pickle_interpolator} in {end:.03f} seconds for {self.interpolator.size()} LS points")
        except Exception as e:
            LOGGER.warning(f"Could not load pre-computed interpolator. Computing interpolator")
            interpolator = get_interpolator(self.sampler_interpolator)
            self.interpolator = interpolator(
                states = states,
                kneigh = int(self.config["sampler"].get("kneighbors", 10)),
                lowerbound = float(self.config["sampler"]["lambda_lowerbound"]),
                upperbound = float(self.config["sampler"]["lambda_upperbound"]),
                num_iter_max = int(self.config["sampler"].get("num_iter_max", 1000000))
            )
            end = time.time() - start
            LOGGER.info(f"Interpolator {self.sampler_interpolator} created in {end:.03f} seconds for {self.interpolator.size()} LS points")

        return ls_sampler.LSSampler(
            states = states,
            interpolator = self.interpolator,
            feedback_file = feedback_db_path,
            feedback_frame = feedbackframe_db
        )

    def setup_generator(self):
        if self.mini_mummi:
            return FullMiniDenseAutoencoder(
                model_path = self.encoder_path
            )
        else:
            return FullDenseAutoencoder(
                model_path = self.encoder_path
            )
            # @TODO add code to swtich between encoders
            # generator = HierarchicalCGAutoencoder(
            #     model_path = model_path
            # )

    def setup_validator(self):
        iteration_path = self.config["validator"]["outpath"]
        resource_name = self.config["validator"]["resources"]
        complex_name = self.config["validator"]["complex"]
        healing = bool(self.config["validator"]["healing"])
        cleanup = bool(self.config["validator"]["cleanup"])
        resource_path = Naming.dir_res(resource_name)

        return CGValidator_C2(
            iteration_path = iteration_path,
            resource_path = resource_path,
            complex_name = complex_name,
            healing = healing,
            cleanup = cleanup,
            mini_mummi = self.mini_mummi
        )

    def start_essential_components(self, results: dict):
        """
        This function just start the sampler. It is supposed to be
        running in a dedicated thread as the sampler can take a 
        long time to get ready.
        """
        try:
            results["sampler"] = self.setup_sampler()
        except Exception as e:
            traceback.print_exc()
            self.logger.error("ML Server failed during setup!")
            self.error()
            raise e

    # @staticmethod
    # def start_feedback_frames(config: dict, database: str, path: str, model_name: str, mini_mummi: bool = False):
    #     # We have to duplicate model because multiprocessing does not support CUDA/Pytorch
    #     # or we should use spawn which is tricky
    #     if mini_mummi:
    #         model = FullMiniDenseAutoencoder(
    #             model_path = model_name,
    #             device = "cpu"
    #         )
    #     else:
    #         model = FullDenseAutoencoder(
    #             model_path = model_name,
    #             device = "cpu"
    #         )

    #     # Process that watches frames being outputed by ddcmd
    #     feedback_frames = FeedbackFrames(
    #         database = database,
    #         path = path,
    #         model = model,
    #         model_name = model_name
    #     )
    #     feedback_frames.start()

    # --------------------------------------------------------------------------
    def setup(self) -> None:
        # if multiprocessing.current_process().name == "MainProcess":
        if self._setup.is_set():
            return
        try:
            self.read_specs()

            self.rpc_server = RPCServer(
                credentials = self.credentials_path,
                cacert = self.certificate_path,
                queue = self.queue,
                logger = self.logger
            )
            self.sampler = None
            self.generator = self.setup_generator()
            self.validator = self.setup_validator()
            # # Will gather result
            self.config["manager"] = self._manager.dict()
            # now we start in another process for the sampler so we do not block the ML server
            process_generator = multiprocessing.Process(target=self.start_essential_components, args=(self.config["manager"],))
            process_generator.start()
            self.config["process_generator"] = process_generator # We save the process to wait on it in the callback.
        except Exception as e:
            traceback.print_exc()
            self.logger.error("ML Server failed during setup!")
            self.error()
            raise e
        
        # Set that the deamon is now setup.
        self._setup.set()

        # if self.do_feedback:
        #     sims_cg = Naming.dir_root("all-cg")
        #     self.feedback = multiprocessing.Process(
        #         target = MLServerRPC.watch_createsims,
        #         args = (self.sampling_db, self.encoder_path, sims_cg, self.filelock, self.waittime, self.logger,)
        #     )

        #     # Process that watches frames being outputed by ddcmd
        #     tmp_config = self.config.copy()
        #     feed_path = Naming.dir_root("feedback-cg")
        #     self.feedback_frames = multiprocessing.Process(
        #         target = MLServerRPC.start_feedback_frames,
        #         args = (tmp_config, self.feedbackframe_db, feed_path, self.encoder_path, self.mini_mummi,)
        #     )

    def run(self) -> None:
        p = multiprocessing.current_process()
        signal.signal(signal.SIGTERM, self.signal_wrapper("mlserver", p.pid))
        signal.signal(signal.SIGINT, self.signal_wrapper("mlserver", p.pid))
        # We start listening for commands (RPC) from the Workflow Manager
        input_config = self.config.copy()
        input_config["obj_generator"] = self.generator
        input_config["obj_validator"] = self.validator

        # if self.do_feedback:
        #     self.feedback.start()
        #     self.logger.info(f"Started process {self.feedback} to watch createsims")

        #     self.feedback_frames.start()
        #     self.logger.info(f"Started process to process feedback frames in {self.feedbackframe_db}")

        self.rpc_server.connect()
        self.rpc_server.run(
            callback = generate_new_samples,
            args = input_config
        )

    # @staticmethod
    # def watch_createsims(sampling_db: str, encoder_path: str, createsims_dir: str, filelock: str, wait: int, logger: Logger):
    #     """
    #     Check all directories in createsims_dir to see if they contains createsims_success or createsims_failure.
    #     then update sampling_db file. This function will wait for "wait" seconds between each check.
    #     """
    #     createsims_status = {}
    #     file_patterns = ["createsims_success", "createsims_failure"]

    #     if not os.path.isfile(sampling_db):
    #         create_sampling_db(sampling_db, encoder_path)
        
    #     logger.info(f"Start monitoring for feedback {sampling_db}")

    #     # File we will be writing (just in case to not corrupt the file in case of interruption)
    #     timestr = time.strftime("%Y%m%d-%H%M%S")
    #     db_path = os.path.splitext(sampling_db)
    #     database_tmp = f"{db_path[0]}-{timestr}{db_path[1]}"

    #     while True:
    #         with FileLock(filelock):
    #             start = timer()
    #             for folder in os.scandir(createsims_dir):
    #                 logger.debug(f"Checking {folder.path}")
    #                 struct_name = os.path.basename(folder.path)
    #                 success = os.path.exists(os.path.join(folder.path, file_patterns[0]))
    #                 failure = os.path.exists(os.path.join(folder.path, file_patterns[1]))
    #                 if success:
    #                     createsims_status[struct_name] = 1 # True
    #                     logger.debug(f"Found success for {struct_name}")
    #                 if failure:
    #                     createsims_status[struct_name] = 0 # False
    #                     logger.debug(f"Found failure for {struct_name}")
    #                 if success and failure:
    #                     logger.warning(f"Found success and failure for {struct_name}")
    #             try:
    #                 with np.load(sampling_db, allow_pickle=True) as data:
    #                     structure_names = data['structure_names']
    #                     previous_createsims_status = data['createsims_status']
    #                     model_name = data['model_name']
    #                     created_at = data['created_at']
    #                     ls_coords = data['ls_coords']
    #                     lambda_values = data["lambda_values"]
    #                     validation_status = data['validation_status']
    #             except Exception as e:
    #                 logger.error(f"{sampling_db} {e}")
    #                 break

    #             # Update existing createsims status
    #             for updated_struct, status in createsims_status.items():
    #                 logger.debug(f"Updating {updated_struct} status={status}")
    #                 index = np.where(structure_names == updated_struct)[0]
    #                 if len(index) == 0:
    #                     logger.warning(
    #                         f"{updated_struct} is running as a createsims but "
    #                         f"has not been sampled by that sampler. It should not "
    #                         f"happened unless you mixed feedback DBs from two different runs. Index = {index}")
    #                 elif len(index) == 1:
    #                     logger.debug(f"Updated {updated_struct} from {previous_createsims_status[index]} => {status}")
    #                     previous_createsims_status[index] = status
    #                 else:
    #                     logger.error(f"{updated_struct} has duplicate in {sampling_db}. The DB is likely corrupted. (Index={index})")

    #             np.savez_compressed(
    #                 database_tmp,
    #                 model_name = model_name,
    #                 created_at = created_at,
    #                 updated_at = datetime.datetime.now().strftime('%d.%m.%Y-%H:%M:%S'),
    #                 structure_names = structure_names,
    #                 ls_coords = ls_coords,
    #                 lambda_values = lambda_values,
    #                 validation_status = validation_status,
    #                 createsims_status = previous_createsims_status
    #             )
    #             # We make sure the new file is not corrupted somehow
    #             try:
    #                 with np.load(database_tmp, allow_pickle=True) as test:
    #                     logger.debug(f"{database_tmp} is valid {test.files}")
    #             except Exception as e:
    #                 logger.warning(f"{database_tmp} seems to be corrupted. We keep the old {sampling_db} intact")
    #                 return
    #             os.replace(database_tmp, sampling_db)
    #             end = timer() - start
    #             if len(createsims_status) > 0:
    #                 logger.info(f"Updated {len(createsims_status)} structures in {end:.3f} seconds ({len(createsims_status)/end} struct/sec)")
    #         time.sleep(wait)
