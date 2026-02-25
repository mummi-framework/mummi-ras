#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------

import os
import multiprocessing
import signal
import traceback
import numpy as np
from typing import List, Union
import time
import yaml

import mummi_core
from mummi_core.brokers import RPCServer
from mummi_core.workflow.flux_env import flux_uri
from mummi_core.utils import MummiFilterLogger

from mummi_ras import Naming
from mummi_ras.ucg.generators import UCGGenerator
from mummi_ras.ucg.validators import UCGValidator
from mummi_ras.ucg.selectors import UCGSelector, UCGSelectorType

import logging
from logging import Logger
LOGGER = logging.getLogger(__name__)

# filter = MummiFilterLogger(str_to_filter=["MDAnalysis", "mummi_ras.dynim", "pika"])
# LOGGER.addFilter(filter)

def generate_new_samples(**config: dict) -> List[str]:
    """
    Generate k_samples from UCG structures, validate and select them.
    Return the structures that passed the validation. The input
    parameters are the following:
        {'iteration_id':int ,'k_samples': int, 'strict': bool}
    where:
        - iteration_id: represents the iteration ID of the 
                        current simulation campaign.
        - k_samples:    Number of requested samples.
        - strict:       if True, the returned List is of length k_samples.
                        If False, len(List) <= k_samples.

    :param config: Hashmap containing the various input parameters
    :type fname: dict

    :return: List of valid and selected structures
    :rtype:  List[str]
    """
    start_init = time.time()

    mummi_core.init()
    model_path = os.path.join(config["encoder"]["path"], config["encoder"]["model"])
    outpath = config["generator"]["outpath"]

    iteration_path = config["validator"]["outpath"]
    resource_path = Naming.dir_res("ucg-validator")
    data_path = Naming.dir_res("ucg")
    trajectory_path =  os.path.join(config["generator"]["outpath"], config["generator"]["trajectories"])
    # --------------------------------------------------------------------------

    # Number of valid samples we need before doing a selection
    # Example: if equals to 10, then to select k sampels we need 2*10 valid candidates
    factor_down_selection = 10
    # This factor defines the factor of extra structure we will generate
    # Example: if equals to 2, we will select up to 2 * nMaxSelectedPatchBuffer structure at most
    # It is useful in case a lot of createsim are failing (which the UCG server does not know about)
    factor_extra_structures = 4

    # Internal parameters
    iteration_id = int(config["body"].get("iteration_id"))
    k_samples = max(int(config["body"].get("k_samples")), 0)

    LOGGER.info(f"Starting generating {k_samples} samples for round {iteration_id}")

    generator = UCGGenerator(
        iteration_id = iteration_id,
        trajectory_path = trajectory_path,
        data_path = data_path,
        model_path = model_path,
        outpath = outpath
    )
    validator = UCGValidator(
        iteration_path = iteration_path,
        resource_path = resource_path,
        cleanup = True
    )
    # nb_cores = 12

    # Selection process
    iteration_id_str = "{:02d}".format(iteration_id)
    selection_path = os.path.join(config["generator"]["outpath"], f"iter{iteration_id_str}/selection")
    ucgconfig = Naming.ucg()
    assert ucgconfig['selection'] == 'importance'

    selector = UCGSelector(UCGSelectorType.Manager, 'manager', selection_path, ucgconfig)
    selector.restore()
    selector.test()
    min_cands_b4_sel = int(ucgconfig['min_cands_b4_sel'])
    num_candidates = selector.num_candidates()
    end_init = time.time() - start_init
    LOGGER.info(f"Initialization of sampling method in {end_init:.3f} secs")

    # --------------------------------------------------------------------------
    if Naming.MUMMI_MINI:
        wfmanager_yaml = os.path.join(Naming.MUMMI_SPECS, 'mini-mummi')
    else:
        wfmanager_yaml = os.path.join(Naming.MUMMI_SPECS, 'workflow')
    with open(os.path.join(wfmanager_yaml, 'wfmanager.yaml'), 'r') as data:
        wfmngr = yaml.load(data, Loader=yaml.FullLoader)
    try:
        nMaxSelectedPatchBuffer = int(wfmngr["wfmanager"]["config"].get("nMaxSelectedPatchBuffer"))
    except Exception as _:
        nMaxSelectedPatchBuffer = 0

    # We do not need to select more than what is defined in wfmanager.yaml
    if selector.num_selected() >= factor_extra_structures * nMaxSelectedPatchBuffer:
        LOGGER.info(f"{selector.num_selected()}/{factor_extra_structures * nMaxSelectedPatchBuffer} have been selected. returning []")
        return []

    # We generate factor_down_selection * K valid samples from UCG trajectories and then select K structures
    total_valid_files = []
    loop_iter = 0
    while len(total_valid_files) <= factor_down_selection * k_samples:
        start_loop = time.time()
        # Load the trajectories file and load at most factor_down_selection*k_samples
        trajectories_ids, raw_data = generator.load_trajectories(max_traj_load = factor_down_selection*k_samples)
        ucg_coords = generator.process_trajectories(raw_data)
        new_positions = generator.new_positions(coords = ucg_coords)

        names_array, positions_array = generator.write_patches(
            iteration_id = iteration_id,
            new_positions = new_positions,
            trajectories_ids = trajectories_ids
        )
        LOGGER.debug(f"generated new_positions = {new_positions.shape}")
        # --------------------------------------------------------------------------

        LOGGER.debug(f"[{loop_iter}] Patches = {names_array}")
        LOGGER.info(f"[{loop_iter}] Checking {len(set(names_array))} / {len(names_array)} structures.")

        for i in range(len(positions_array)):
            # res is under the form => [valid_status, gromacs_file, mf, pe]
            res = validator.validate(iteration_id, names_array[i], positions_array[i])
            _, _ = validator.write_validation_info(
                    iteration_id = iteration_id,
                    return_array = np.array([res]),
                    all_structure_names = names_array[i]
            )
            if res[0]:
                valid_file = os.path.join(validator.current_rpath, res[1])
                total_valid_files.append(valid_file)

                selector.add_snapshot(structure_file = valid_file, topology_file = os.path.splitext(valid_file)[0]+'.tpr', cleanup = True)
                num_candidates = selector.num_candidates()

                LOGGER.debug(f"Candidates ({num_candidates}) in selector "
                        f"(we need >= {factor_down_selection * k_samples}) "
                        f", min_cands_b4_sel = {min_cands_b4_sel} "
                        f"#valid patches in iter {loop_iter} => {total_valid_files}")

                # If enought in selector we can select right now
                if num_candidates >= min_cands_b4_sel and num_candidates >= factor_down_selection * k_samples:
                    LOGGER.info(f"Enough candidates ({num_candidates}) in selector to generate {k_samples} samples.")

                    selections = [s.id for s in selector.select(k_samples)]
                    end_total = time.time() - start_init

                    LOGGER.info(f"Total time to generate {len(selections)} samples => {end_total:.3f} secs")
                    LOGGER.info(f"[{loop_iter}] selections({len(selections)})={selections}")
                    return selections

        end_loop = time.time() - start_loop
        LOGGER.info(f"Iteration {loop_iter} selector: number of candidates = {selector.num_candidates()} in {end_loop:.3f} secs")
        loop_iter += 1

    # --------------------------------------------------------------------------
    LOGGER.info(f"[{loop_iter}] {len(total_valid_files)} valid_files")

    selections = [s.id for s in selector.select(k_samples)]
    end_total = time.time() - start_init
    LOGGER.info(f"Total time to generate {len(selections)} samples => {end_total:.3f} secs")

    LOGGER.info(f"[{loop_iter}] selections({len(selections)})={selections}")
    return selections

class UCGServerRPC(object):
    """
    UCG server implementation with RPC support 
    """
    def __init__(self, config: dict, logger: Logger) -> None:
        self.config = config
        self.rpc_server: RPCServer = None
        self.logger = logger
        self.hostname = mummi_core.get_hostname(contract_hostname=False)

        self._exit = multiprocessing.Event()  # Trigger the daemon to exit.
        self._error = multiprocessing.Event() # Trigger the daemon to exit.
        self._setup = multiprocessing.Event() # Wait until the daemon is setup.

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
            self.rpc_server.stop()
            self.rpc_server.purge()
            self.rpc_server.close()
            exit(1)
        return handler

    # --------------------------------------------------------------------------
    def setup(self) -> None:
        if self._setup.is_set():
            return

        try:
            # ------------------------------------------------------------------
            # identify scheduler interface
            self.logger.info(f'UCG Server launched on host: {self.hostname}')

            self.flux = flux_uri(override_uri_from_file=False)
            self.logger.info(f"  flux  uri: [{self.flux}]")
            flux_uri_local = os.environ.get("FLUX_URI", None)
            self.logger.info(f"  local uri: [{flux_uri_local}]")

            if self.flux is None:
                self.logger.error("Failed to identify Flux URI. Exiting ...")
                assert False, "Failed to identify flux uri"

            # ------------------------------------------------------------------
            # Message broker initialization (e.g., RabbitMQ)
            # ------------------------------------------------------------------
            credentials_path = os.path.join(
                self.config["workspace"]["path"],
                self.config["workspace"]["credentials"]
            )

            certificate_path = os.path.join(
                self.config["workspace"]["path"], 
                self.config["workspace"]["certificate"]
            )

            self.encoder_name = self.config["encoder"]["model"]

            self.encoder_path = os.path.join(
                self.config["encoder"]["path"], 
                self.encoder_name
            )
            self.queue = self.config["broker"]["queue"] + '_' + os.getlogin()

            self.broker_interface = self.config["broker"]["interface"]
            # interface = get_broker(self.broker_interface)
            self.rpc_server = RPCServer(
                credentials = credentials_path,
                cacert = certificate_path,
                queue = self.queue,
                logger = self.logger
            )

            self.logger.info(f"> Initializing MuMMI UCG Server")
            self.logger.info(f"  > Server")
            self.logger.info(f"    > Interface      {self.broker_interface}")
            self.logger.info(f"    > Credentials    {credentials_path}")
            self.logger.info(f"    > Certificate    {certificate_path}")
            self.logger.info(f"    > Queue          {self.queue}")
            self.logger.info(f"  > Sampler")
            self.logger.info(f"    > Encoder        {self.encoder_path}")

        except Exception as e:
            traceback.print_exc()
            self.logger.error("UCG Server failed during setup!")
            self.error()
            raise e

        # Set that the deamon is now setup.
        self._setup.set()

    def run(self) -> None:
        p = multiprocessing.current_process()
        signal.signal(signal.SIGTERM, self.signal_wrapper("ucgserver", p.pid))
        signal.signal(signal.SIGINT, self.signal_wrapper("ucgserver", p.pid))

        # We start listening for commands (RPC) from the Workflow Manager
        input_config = self.config.copy()
        self.rpc_server.connect()
        self.rpc_server.run(
            callback = generate_new_samples,
            args = input_config
        )
