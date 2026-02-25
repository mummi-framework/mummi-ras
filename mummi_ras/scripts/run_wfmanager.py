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

import os
import yaml
import sys
import signal
import traceback
import numpy as np
from os.path import join

import mummi_core
from mummi_core.utils import Naming
from mummi_ras.workflow.wfmanager import WFManager

from logging import getLogger
LOGGER = getLogger(__name__)

# recursive flux scheduling is requiring reset of the TMPDIR
# https://github.com/flux-framework/flux-core/issues/4328
# https://github.com/flux-framework/flux-core/issues/4329

# TODO: this should be set from the script that launches the job
tmpdir = '/tmp'
os.makedirs(tmpdir, exist_ok=True)
os.environ["TMPDIR"] = tmpdir

# ------------------------------------------------------------------------------
def read_specs(spath):
    # read the specs
    LOGGER.info(f'Reading specs from ({spath})')
    with open(os.path.join(spath, 'wfmanager.yaml'), 'r') as data:
        wfmngr = yaml.load(data, Loader=yaml.FullLoader)

    with open(os.path.join(spath, 'maestro.yaml'), 'r') as data:
        maestro = yaml.load(data, Loader=yaml.FullLoader)

    with open(os.path.join(spath, 'pf2pfpatches.yaml'), 'r') as data:
        pf2pfpatches = yaml.load(data, Loader=yaml.FullLoader)

    config = wfmngr
    config['wfmanager']['batch'] = maestro['maestro']['batch']

    # macro gc sims now run within the createsim
    config['macro_patch_creator'] = pf2pfpatches['macro_patch_creator']

    print(config)

    with open(os.path.join(spath, 'jobs_createsim.yaml'), 'r') as data:
        config['createsim'] = yaml.load(data, Loader=yaml.FullLoader)
        if config["wfmanager"]["config"].get("do_ucgserver") and bool(config["wfmanager"]["config"]["do_ucgserver"]):
            config["createsim"]["variables"]["inpath"]["eval"] = "mummi_ras.Naming.dir_root('ucg')"
        elif config["wfmanager"]["config"].get("do_mlserver") and bool(config["wfmanager"]["config"]["do_mlserver"]):
            config["createsim"]["variables"]["inpath"]["eval"] = "mummi_ras.Naming.dir_root('ml')"

    with open(os.path.join(spath, 'jobs_cg.yaml'), 'r') as data:
        config['cg'] = yaml.load(data, Loader=yaml.FullLoader)

    # # For GROMACS, only use for campaign 1 because ddcMD is not working!!
    # with open(os.path.join(spath, 'jobs_gmx.yaml'), 'r') as data:
    #     config['cg'] = yaml.load(data, Loader=yaml.FullLoader)

    with open(os.path.join(spath, 'jobs_backmapping.yaml'), 'r') as data:
        config['backmapping'] = yaml.load(data, Loader=yaml.FullLoader)

    with open(os.path.join(spath, 'jobs_aa.yaml'), 'r') as data:
        config['aa'] = yaml.load(data, Loader=yaml.FullLoader)

    # --------------------------------------------------------------------------
    # read list of macro gc sims and override the config
    # macro gc sims now run within the createsim
    simlist = []
    try:
        simlistfile = os.path.join(Naming.dir_root('macro'), 'simlist.spec')
        simlist = np.genfromtxt(simlistfile, delimiter=',', dtype=str)
        simlist = [os.path.splitext(s)[0] for s in simlist[:,0]]
        simlist = [s for s in simlist if s != 'large']
        print(f'> Read ({len(simlist)}) simulations from ({simlistfile})')
    except Exception as e:
        print(f'> Failed to read macro simulation list: {e}')

    config['macro_patch_creator']['params']['sim2frame'] = {s: 1 for s in simlist}
    print(f'> Updated sim2frame: {config["macro_patch_creator"]["params"]}')

    # --------------------------------------------------------------------------
    return config


# ------------------------------------------------------------------------------
def term_wrapper(job):
    def signal_handler(signum, frame):
        LOGGER.info(f"Signal caught! {signum} -- Terminating")
        job.stop()
    return signal_handler

# ------------------------------------------------------------------------------
def main():
    # --------------------------------------------------------------------------
    hostname = mummi_core.get_hostname(contract_hostname=False)
    print(f"> Launching workflow manager on ({hostname})")

    mummi_core.init()
    mummi_core.create_root()

    # --------------------------------------------------------------------------
    if Naming.MUMMI_MINI:
        config = read_specs(os.path.join(Naming.MUMMI_SPECS, "mini-mummi"))
    else:
        config = read_specs(os.path.join(Naming.MUMMI_SPECS, "workflow"))
    # --------------------------------------------------------------------------
    # create logger
    mummi_core.init_logger(config = config['wfmanager']['config'])

    # --------------------------------------------------------------------------
    # manager = Manager()
    # exceptions = manager.Queue()
    job = WFManager('wfmanager') # , exceptions)

    # --------------------------------------------------------------------------
    # start the job
    try:
        env = None
        job.setup(env, config)
        job.start()

        # Add signals for termination
        # signal.signal(signal.SIGTERM, term_wrapper(job))
        signal.signal(signal.SIGINT, term_wrapper(job))

        # Join threads after start
        job.join()

    except Exception as e:
        LOGGER.error(f"Exiting due to error ({e})")
        traceback.print_exc()
        # job.stop()
    finally:
        sys.exit(1)


# ------------------------------------------------------------------------------
if __name__ == '__main__':
    main()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
