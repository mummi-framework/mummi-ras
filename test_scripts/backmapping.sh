#!/bin/bash

#spack environment setup.
source ~/.mummi/config.mummi.sh
source $MUMMI_APP/setup/setup.env.sh

spack unload py-mdanalysis-mummi@mda0.19.2_ddcmd_lessmem
spack load -r py-mdanalysis-mummi@mda_1.0.1_ddcmd
spack load py-tqdm
spack load -r py-parmed

spack load -r gromacs@2019.6.mummifix +double
spack load -r gromacs@2019.6.mummifix ~double

# the test input is in
# /g/g92/zhang30/W/MuMMI/scratch/mummi_root_20200802/sims-cg/pfpatch_000000384620
# require input files: topol.tpr system.top positions.tar

# ML module need to provide patch id (384620) and frame id(100000)
cd $MUMMI_ROOT
python $MUMMI_APP/mummi/transformations/backmapping/backmapping.py -p pfpatch_000000384620 -f 100000

# If you want to run with specified outlocal path use "--outlocal" option
# python backmapping.py -p 384620 -f 100000 --outlocal /tmp/mummi/backmapping/pfpatch_000000384620_20200820-233210
