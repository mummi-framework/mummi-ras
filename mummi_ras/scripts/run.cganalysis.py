# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
# start a persistent job splash app to analyze a running CG simulation
import os
import sys
import yaml

from mummi_core import Naming
import mummi_ras


if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('Usage:', sys.argv[0], ' 0/1')
        print('\t 0: run locally')
        print('\t 1: run in pbatch')
        exit()

    rtype = int(sys.argv[1])
    if rtype != 0 and rtype != 1:
        print('Usage:', sys.argv[0], ' 0/1')
        print('\t 0: run locally')
        print('\t 1: run in pbatch')
        exit()

    mummi_ras.init()

    # configuration for pfpatches
    rpath = os.path.join(Naming.MUMMI_RESOURCES, 'wfspecs/wfmanager.yaml')
    with open(rpath, 'r') as data:
        pyaml = yaml.load(data)

    env = pyaml['cganalysis']['env']
    config = pyaml['cganalysis']['config']
    params = pyaml['cganalysis']['params']
    batch = pyaml['wfmanager']['batch']

    # empty params will be read from common configuration
    if len(params['nscgpatches']) == 0:
        params['nscgpatches'] = Naming.dbr('cgpatches')

    if len(params['nscg2bkilled']) == 0:
        params['nscg2bkilled'] = Naming.dbr('cg2bkilled')

    if len(config['workspace']) == 0:
        config['workspace'] = Naming.dir_root('workspace')

    # hardcoded parameters
    # params['simname'] = 'pfpatch_000000000000'
    # params['inpath'] = '/p/lscratchd/bhatia4/pilot2_data/test'

    params['simname'] = 'pfpatch_000000000048'
    params['inpath'] = \
        '/p/lscratchd/bhatia4/splash-run-180313/sims/' \
        'pfpatch_000000000048/cg-prod/'

    if rtype == 0:

        cmd = config["jobbin"]
        for var, value in list(params.items()):

            if var == 'inpath':
                p = value
                p = os.path.expandvars(p)    # expands environment variables
                if p[0] == '~':
                    p = os.path.expanduser(p)       # expands home directory
                p = os.path.abspath(p)
                value = p

            cmd = cmd + ' --' + var + ' ' + str(value)

        print(cmd)
        os.system(cmd)

    else:

        from queue import Queue
        from core.workflow.modules import MLPersistentJob

        exceptions = Queue()
        job = MLPersistentJob('cganalysis', batch, exceptions)
        job.setup(env, config, params)
        job.start()
