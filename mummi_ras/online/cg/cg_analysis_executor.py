"""
Python script for continuously monitoring and running analysis
"""

from mpi4py import MPI
import pathlib
import time
import argparse
import os
import sys

import mummi_core
from mummi_core.utils import Naming
from  mummi_ras.online.cg import cganalysis as cga

from logging import getLogger
LOGGER = getLogger(__name__)

# figure out list of simulations to run
sims_cg_dir = pathlib.Path(Naming.dir_sim('cg'))


if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nprocs = comm.Get_size()

    # Set up the necessary base data structures to begin study set up.
    parser = cga.setup_argparser()
    args = parser.parse_args()

    ## Sets up MuMMI and the logger
    mummi_core.init()
    #mummi_core.create_root()
    mummi_core.init_logger(argparser=args)

    # run a finite number of loops through the scan
    for i in range(2000):
        sims = list(sims_cg_dir.glob("structure*"))
        for sim in sims[rank::nprocs]:

            # detect how much analysis was done
            analysis_pylst = os.path.join(sim, "analysis.tar.pylst")
            if os.path.isfile(analysis_pylst):
                with open(analysis_pylst, "r") as fh:
                    lines = fh.readlines()
                    fnum = lines[-1][10:21]

            else:
                fnum = 0

            # edit args for this run
            args.simname = sim.name
            args.path = sim
            args.fcount = int(fnum)

            try:
                cga.run_cganalysis(args)
            except Exception as e:
                LOGGER.error(e)
                #traceback.print_exc()
                pass
                #raise e

        comm.Barrier()
        # wait two minutes before doing another pass
        time.sleep(120)

    MPI.Finalize()
