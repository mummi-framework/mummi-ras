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

import argparse
import os
import signal
import sys
import yaml
import timeit
import traceback
import multiprocessing
from multiprocessing.managers import SyncManager
import pika

import mummi_core
from mummi_core.brokers import RPCClient

import mummi_ras
from mummi_ras import Naming

from logging import getLogger
LOGGER = getLogger(__name__)

def main(nPatches: int):
    mummi_core.init()
    mummi_core.create_root()

    if Naming.MUMMI_MINI:
        spath = os.path.join(mummi_ras.Naming.MUMMI_SPECS, "mini-mummi")
    else:
        spath = os.path.join(mummi_ras.Naming.MUMMI_SPECS, "workflow")

    config = mummi_core.read_specs(os.path.join(spath, "mlserver.yaml"))
    mummi_core.init_logger(config = config["config"])

    wspace = config["workspace"]["path"]
    credentials_path = os.path.join(wspace,
        config["workspace"]["credentials"]
    )   
    certificate_path = os.path.join(wspace,
        config["workspace"]["certificate"]
    )

    usr = os.environ.get("USER", "mummiusr")
    routing_key = config["broker"]["queue"] + '_' + usr
    if usr == "mummiusr":
        LOGGER.warning(f"Did not find current user: defaulted to {usr}")
    broker_interface = config["broker"]["interface"]
    try:
        # interface = get_broker(broker_interface)
        rpc_client = RPCClient(
            credentials = credentials_path,
            cacert = certificate_path,
            routing_key = routing_key,
        )
    except Exception as e:
        traceback.print_exc()
        LOGGER.error("ML Server failed during setup!")
        raise e

    LOGGER.info(f"> Initializing RPC client to interact with ML Server")
    LOGGER.info(f"  > Client")
    LOGGER.info(f"    > Interface   {broker_interface}")
    LOGGER.info(f"    > Credentials {credentials_path}")
    LOGGER.info(f"    > Certificate {certificate_path}")
    LOGGER.info(f"    > Routing key {routing_key}")

    iterCounterMLServer = 0
    hostname = mummi_core.get_hostname()

    # select new patches
    p = multiprocessing.current_process()
    LOGGER.debug(f"> Running RPC client {p.name} on [{hostname}/{p.pid}]")
    LOGGER.info(f" [x] Requesting generate_new_samples(iter={iterCounterMLServer}, k={nPatches})")
    # If strict equal True then call() will ALWAYS return the requested number of VALID patches (but it will take longer to run)
    selections = []
    try:
        selections = rpc_client.call({'k_samples': nPatches, 'iteration_id': iterCounterMLServer, 'strict': False})
        if isinstance(selections, str):
            LOGGER.error(f"[.] RPC server sent back an exception: {selections}")
            return -1

        selections = [x.split('/')[-1] for x in selections]
        LOGGER.info(f" [.] Got generate_new_samples(iter={iterCounterMLServer}, k={nPatches}) = {selections}")

    except pika.exceptions.AMQPError as e:
        # we try to restart the connection
        LOGGER.error(f"We got {e}")
        return 0


# -----------------------------------------------------------------------------
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script to test the ML server, it sends a request to a ML server running via RabbitMQ. Useful for testing")
    parser.add_argument("-n", "--nsamples", help="Number of samples to generate", required=True, type=int)
    args = parser.parse_args()

    d = main(
        nPatches = args.nsamples
    )
    
