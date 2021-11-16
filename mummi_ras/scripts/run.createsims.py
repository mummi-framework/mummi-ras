# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
# start a persistent job splash app to analyze a running CG simulation
import os
import yaml

from mummi_core import Naming
import mummi_ras
from maestrowf.abstracts.enums import SubmissionCode
from maestrowf.datastructures.core import StudyStep
from maestrowf.interfaces import ScriptAdapterFactory
if __name__ == '__main__':

    '''
    if len(sys.argv) != 2:
        print 'Usage:', sys.argv[0], ' 0/1'
        print '\t 0: run locally'
        print '\t 1: run in pbatch'
        exit()

    rtype = int(sys.argv[1])
    if rtype != 0 and rtype != 1:
        print 'Usage:', sys.argv[0], ' 0/1'
        print '\t 0: run locally'
        print '\t 1: run in pbatch'
        exit()
    '''
    mummi_ras.init()

    # configuration for pfpatches
    # rpath = os.path.join(Naming.RESOURCES, 'wfspecs/wfmanager.yaml')
    rpath = os.path.join(Naming.MUMMI_SPECS, 'workflow/wfmanager.yaml')
    with open(rpath, 'r') as data:
        pyaml = yaml.load(data)

    env = pyaml['createsim']['env']
    conf = pyaml['createsim']['config']
    params = pyaml['createsim']['params']
    batch = pyaml['wfmanager']['batch']

    nprocs = int(conf['nprocs'])
    nnodes = int(conf['nnodes'])

    adapter = \
        ScriptAdapterFactory.get_adapter(batch['type'])(**batch)

    outdir = Naming.dir_root('workspace')

    cmd = conf['jobbin']
    for var, value in list(params.items()):

        if var == 'mpi':
            value = '\"$(LAUNCHER)[{}n, {}p] ' + value + '\"'
            # value = '\"srun -N{} -n{} ' + value + '\"'
        if var == 'logpath' and len(value) == 0:
            value = Naming.dir_root('workspace')
        if var == 'inpath' and len(value) == 0:
            value = Naming.dbr('pfpatches')

        cmd = cmd + ' --' + var + ' ' + str(value)

    cmd = cmd.format(nnodes, nprocs)

    if len(conf['redirect']) > 0:
        cmd = cmd + ' >> ' + conf['redirect'] + ' 2>&1'

    print('\ncreatesim command [[', cmd, ']]')

    patch = params['patch']

    step = StudyStep()
    step.name = conf['jobname'] + '-' + patch
    step.description = conf['jobdesc'].format(patch)
    step.run["cmd"] = cmd
    step.run["procs"] = nprocs
    step.run["nodes"] = nnodes
    step.run["walltime"] = conf['walltime']

    to_be_scheduled, cmd_script, restart_script = \
        adapter.write_script(outdir, step)

    print('wrote script to', outdir)

    exit()
    print('now submitting job!')
    try:
        retcode, jobid = adapter.submit(step, cmd_script, outdir)
        if retcode != SubmissionCode.OK:
            raise Exception("Failed to submit cganalysis to queue.")

    except Exception as e:
        # self._exception_queue.put(sys.exc_info())
        raise e

    print(retcode, jobid)
