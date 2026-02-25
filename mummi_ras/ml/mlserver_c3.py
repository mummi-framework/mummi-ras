#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2025, Lawrence Livermore National Security, LLC. All rights
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
import pathlib
from typing import List, Union, Any, Tuple
import yaml
import fcntl
import time
import datetime
import tempfile
import pickle
import glob
from filelock import FileLock
from timeit import default_timer as timer

import torch
import numpy as np
import numpy.typing as npt

import mummi_core
from mummi_core.brokers import RPCServer
from mummi_core.workflow.flux_env import flux_uri

from mummi_ras import Naming
from mummi_ras.ml.models import V_Generator, UCG2CG_Generator, V2UCG_Generator, write_patches
from mummi_ras.ml.validators import CGValidator_C3 as CGValidator

from mummi_ras.feedback import LSSamplingDB
from mummi_ras.ml.samplers import LSSampler

from mummi_ras.ml.models.autoencoders import FullDenseAutoencoderUCG as DenseAutoEncoder

from logging import getLogger, Logger
LOGGER = getLogger(__name__)

def generate_new_samples(**config: dict) -> List[str]:
    """
    Sample k_sample from the latent space, generate the k_samples corresponding
    structures and validate them. Return the structures that passed the validation.
    """

    iteration_path = config["validator"]["outpath"]
    iteration_id = config["iteration_id"]
    k_samples = max(int(config["body"].get("k_samples")), 0)
    strict = bool(config["body"].get("strict"))
    # Get back objects
    v2ucg_generator = config["obj_v2ucg_generator"]
    ucg2cg_generator = config["obj_ucg2cg_generator"]
    structure_database_path = config["structure_database"]
    validator = config["obj_validator"]
    lssampling_db = config["obj_sampling_db"]

    sampler = LSSampler(
        structures_db = structure_database_path,
        feedback_file = config["sampling_db"],
        feedback_frame = config["feedbackframe_db"],
        stopping_buffer_size = 10
    )

    # --------------------------------------------------------------------------
    nMaxSelectedPatchBuffer = int(config["nMaxSelectedPatchBuffer"])
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

    checkpoint_v2ucg_paths, checkpoint_ucg2cg_paths = config["gen_models"]

    new_structure_names = []
    new_validation_status = []

    total_valid_files = []
    loop_iter = 0
    num_new_sample = 0

    while len(total_valid_files) < k_samples:
        gen_time = time.time()
        LOGGER.info(f"Selecting {k_samples} structures")
        selected_ids, v_coords, v2ucg_coords, pred_cg = sampler.get_new_ls_points(k_samples)
        gen_time = time.time() - gen_time

        if len(selected_ids) == 0:
            break

        names_array, _ = write_patches(
            outpath = iteration_path,
            iteration_id = iteration_id,
            new_positions = pred_cg
        )

        LOGGER.info(f"Selection: sec/sample: {gen_time/k_samples}")

        # Dict to match the v_coords to the name of the correspondong structure
        match_v_coords = {}
        match_v2ucg_coords = {}
        match_pred_cg = {}

        for struct_id, v, v2ucg, ucg2cg in zip(names_array, v_coords, v2ucg_coords, pred_cg):
            match_v_coords[struct_id] = v
            match_v2ucg_coords[struct_id] = v2ucg
            match_pred_cg[struct_id] = ucg2cg

        return_array = validator.validateArrayThreaded(iteration_id, names_array, pred_cg)
        if len(return_array) > 0:
            new_validation_status.append(validator.get_status(return_array))

        _, valid_files = validator.write_validation_info(
            iteration_id = iteration_id,
            return_array = return_array,
            all_structure_names = names_array
        )

        new_structure_names = validator.get_structure_names(return_array)
        LOGGER.info(f"new_structure_names = {new_structure_names}")

        expended_v_coords = []
        expended_v2ucg_coords = []
        expended_pred_cg = []

        # Because one structure_0000X can generate multiple structure_000X_FF
        # We need to expand the v's, v2ucg and pred_cg to match the different Force-Fields
        for structure in new_structure_names:
            strip_structure = structure.split("_")
            if len(strip_structure) == 4:
                strip_structure = "_".join(structure.split("_")[:-1]) # We remove the force fields
            else:
                # No force field (pre C3 structures)
                strip_structure = structure
            expended_v_coords.append(match_v_coords[strip_structure])
            expended_v2ucg_coords.append(match_v2ucg_coords[strip_structure])
            expended_pred_cg.append(match_pred_cg[strip_structure])

        np.random.shuffle(valid_files)
        valid_files = [os.path.join(validator.current_rpath, f) for f in valid_files]
        LOGGER.debug(f"ITER = {loop_iter} => valid structures={valid_files}")
        total_valid_files += list(valid_files)

        # Update sampling DB for feedback
        lssampling_db.update_db(
            num_new_sample = len(new_structure_names),
            structure_names = new_structure_names.tolist(),
            v_coords = expended_v_coords,
            v2ucg_coords = expended_v2ucg_coords,
            ucg2cg_coords = expended_pred_cg,
            validation_status = new_validation_status
        )
 
        LOGGER.debug(f"lssampling_db = {lssampling_db.print_db()}")

        loop_iter += 1
        num_validated = validator.num_validated(iteration_id)
        if num_validated >= factor_extra_structures * nMaxSelectedPatchBuffer:
            LOGGER.info(f"=> {num_validated} / {factor_extra_structures * nMaxSelectedPatchBuffer} have been selected, returning []")
            break
        if not strict:
            break

    LOGGER.info(f"[{loop_iter}] valid_files({len(total_valid_files)})={total_valid_files}")
    return total_valid_files

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
                    os.remove(self.lockfile)
                except FileNotFoundError as _:
                    pass
            self.exit()
            exit(1)
        return handler

    def select_mlmodels(self, lists_paths):
        """
        The YAML config can provided multiple paths for the same model.
        THis happens when we run on different machines. Models are located on
        different paths. We select the first valid one.
        """
        if not isinstance(lists_paths, list):
            lists_paths = [lists_paths]
        if isinstance(lists_paths, list):
            for epath in lists_paths:
                if os.path.exists(epath):
                    self.logger.info(f"Found valid path {epath}")
                    return epath
        

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

            self.config["feedbackframe_db"] = self.feedbackframe_db
            self.config["sampling_db"] = self.sampling_db

            self.do_feedback = bool(self.config["sampler"]["feedback"]["do_feedback"])

            if "mini_mummi" in self.config["config"]:
                self.mini_mummi = bool(self.config["config"]["mini_mummi"])
            
            if self.mini_mummi:
                raise RuntimeError("ML server for C3 cannot be used with mini-MuMMI")

            self.config["sampling_db"] = self.sampling_db
            self.config["lock_sampling_db"] = self.lockfile
            self.sampler_interpolator = self.config["sampler"]["interpolator"]
            self.pickle_interpolator = self.config["sampler"].get("pre_computed")
            if self.pickle_interpolator and os.path.isfile(self.pickle_interpolator):
                self.logger.info(f"We will use a pre-computed interpolator: {self.pickle_interpolator}")


            # UCG autoencoder for Feedback
            if self.do_feedback:
                encoder_ucg_path = self.select_mlmodels(self.config["models"]["ucg_autoencoder"]["paths"])
                encoder_path = os.path.join(encoder_ucg_path, self.config["models"]["ucg_autoencoder"]["model"])
                indice_cluster = os.path.join(self.config["models"]["ucg_autoencoder"]["resources"], "indices_per_cluster_dimer.npz")
                self.config["obj_fb_autoencoder"] = DenseAutoEncoder(model_path = encoder_path, device = "cpu", indice_per_cluster = indice_cluster)
            else:
                self.logger.info(f"No feedback as do_feedback = False. Deactivating UCG autoencoder")
                self.config["obj_fb_autoencoder"] = None

            # We take the first model that exists
            self.checkpoint_v2ucg_paths = self.select_mlmodels(self.config["models"]["generative"]["v2ucg"]["paths"])
            self.checkpoint_ucg2cg_paths = self.select_mlmodels(self.config["models"]["generative"]["ucg2cg"]["paths"])
            self.config["gen_models"] = [self.checkpoint_v2ucg_paths, self.checkpoint_ucg2cg_paths]

            # We precompute structure in C3
            self.structure_database = self.select_mlmodels(self.config["models"]["generative"]["database"]["paths"])
            self.config["structure_database"] = self.structure_database

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
            self.logger.info(f"    > V2UCG                    {self.checkpoint_v2ucg_paths}")
            self.logger.info(f"    > UCG2CG                   {self.checkpoint_ucg2cg_paths}")
            self.logger.info(f"    > Pre-computed structures  {self.structure_database}")
            self.logger.info(f"  > Validator")
            self.logger.info(f"  > nMaxSelectedPatchBuffer    {self.nMaxSelectedPatchBuffer}")
        except Exception as e:
            traceback.print_exc()
            self.logger.error("ML Server failed during reading specs!")
            self.error()
            raise e

    def setup_sampling_db(self):
        # We declare the LSSampling DB
        return LSSamplingDB(
            database_path = self.sampling_db,
            model_name = self.config["gen_models"],
            ls_encoder = self.config["obj_fb_autoencoder"],
            lockfile = self.lockfile,
            logger = self.logger
        )

    def setup_generators(self):
        """
        This method create the objects needed to generate structures (inference call).
        It returns both ML models
        """

        indice_path = self.config["models"]["generative"]["resources"]

        v2ucg_generator = V2UCG_Generator.load_from_checkpoint(self.checkpoint_v2ucg_paths)
        ucg2cg_generator = UCG2CG_Generator.load_from_checkpoint(self.checkpoint_ucg2cg_paths, ucg_index_file = os.path.join(indice_path, "dimers-ucg-idx.txt"))

        v2ucg_generator.setup('predict')
        ucg2cg_generator.setup('predict')

        return v2ucg_generator, ucg2cg_generator

    def setup_validator(self):
        iteration_path = self.config["validator"]["outpath"]
        resource_name = self.config["validator"]["resources"]
        complex_name = self.config["validator"]["complex"]
        healing = bool(self.config["validator"]["healing"])
        cleanup = bool(self.config["validator"]["cleanup"])
        nthreads = int(self.config["validator"]["nthreads"])
        resource_path = Naming.dir_res(resource_name)

        return CGValidator(
            iteration_path = iteration_path,
            resource_path = resource_path,
            complex_name = complex_name,
            healing = healing,
            cleanup = cleanup,
            nthreads = nthreads,
            mini_mummi = self.mini_mummi
        )

    # --------------------------------------------------------------------------
    def setup(self) -> None:
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
            self.generators = self.setup_generators()
            self.validator = self.setup_validator()
            self.lssampling_db = self.setup_sampling_db()
            # # Will gather result
            self.config["manager"] = self._manager.dict()
            # now we start in another process for the sampler so we do not block the ML server
            # process_generator = multiprocessing.Process(target=self.start_essential_components, args=(self.config["manager"],))
            # process_generator.start()
            # self.config["process_generator"] = process_generator # We save the process to wait on it in the callback.
        except Exception as e:
            traceback.print_exc()
            self.logger.error("ML Server failed during setup!")
            self.error()
            raise e

        # Set that the deamon is now setup.
        self._setup.set()

    def run(self) -> None:
        p = multiprocessing.current_process()
        signal.signal(signal.SIGTERM, self.signal_wrapper("mlserver", p.pid))
        signal.signal(signal.SIGINT, self.signal_wrapper("mlserver", p.pid))
        # We start listening for commands (RPC) from the Workflow Manager
        input_config = self.config.copy()
        input_config["obj_v2ucg_generator"] = self.generators[0]
        input_config["obj_ucg2cg_generator"] = self.generators[1]
        input_config["obj_validator"] = self.validator
        input_config["obj_sampling_db"] = self.lssampling_db

        self.rpc_server.connect()
        self.rpc_server.run(
            callback = generate_new_samples,
            args = input_config
        )
