#!/bin/bash

# ------------------------------------------------------------------------------
# likely to be set for each run
# ------------------------------------------------------------------------------

# path where you want to run mummi (workspace)
export MUMMI_ROOT=~/tst_root_apr29

export MUMMI_REDIS_NNODES=1
export MUMMI_MACRO_NNODES=1
export MUMMI_FAKECG_NNODES=0

# ------------------------------------------------------------------------------
# all these are set up automatically for the group (setup/envs/env.host.sh)
# but they can be overwritten per user
# ------------------------------------------------------------------------------

# virtual environment name of your choice
# no need to change (unless you want to do a complete reinstall)
#export MUMMI_VENV=mummi-ras

# path where you clone mummi repositories
# no need to change (unless you move your repositories)
#export MUMMI_RESOURCES=/gpfs/alpine/lrn005/proj-shared/pilot2/mummi_resources
#export MUMMI_CORE=/gpfs/alpine/lrn005/proj-shared/pilot2/mummi-core
#export MUMMI_APP=/gpfs/alpine/lrn005/proj-shared/pilot2/mummi-ras

# ------------------------------------------------------------------------------
# likey not to be changed
# ------------------------------------------------------------------------------

# group root and spack root
#export MUMMI_GRP_ROOT=/autofs/nccs-svm1_proj/lrn005
#export MUMMI_SPACK_ROOT=/autofs/nccs-svm1_proj/lrn005/spack

# ------------------------------------------------------------------------------
