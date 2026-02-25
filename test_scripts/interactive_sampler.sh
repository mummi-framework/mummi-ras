#!/usr/bin/env bash
#BSUB -nnodes 1
#BSUB -W 30
#BSUB -q pdebug
#BSUB -G cancer
#BSUB -o test_samplers.%J.log
#BSUB -e test_samplers.%J.err
#BSUB -J test_samplers

MAX_MEM=$((256*1024)) #Max memory per node on Lassen 256GB
MEM=$((250*1024))

source ~/.mummi/config.mummi.sh
source $MUMMI_APP/setup/setup.env.sh

spack load gromacs~mpi+mummi-c2+cuda arch=$(spack arch)

#spack load gromacs@2022.06~cuda arch=$(spack arch)
# spack load gromacs@2022.6~cuda arch=$(spack arch)
python3 $MUMMI_APP/tests/test-sampler-validation.py -n "$1"
