#!/usr/bin/env bash
# ------------------------------------------------------------------------------

echo "(`hostname`: `date`) --> Loading Frontier environment"

# Frontier modules (May 2024)
# need to export this because we need to load this again later!
export MUMMI_MPI_MODULE="cray-mpich/8.1.31"

module load gcc-native/12.3
module load rocm/6.3.1
module load cray-fftw/3.3.10

# ------------------------------------------------------------------------------
# Set host properties
export MUMMI_HOST="frontier"
export MUMMI_GROUP="lrn005"

# Setting permisions so every users of $MUMMI_GROUP can access file from the workflow
# newgrp $MUMMI_GROUP
# umask 002

# Explicit paths
#TODO: should be in Spack
export MUMMI_AMBER_PATH="/autofs/nccs-svm1_proj/lrn005/amber18"
export MUMMI_AUTOBIND_PATH="/autofs/nccs-svm1_proj/lrn005/bin:/lustre/orion/lrn005/proj-shared/spack/bin/"

# locale messes up how flux job ids are reported
# https://github.com/flux-framework/flux-docs/issues/61
export LC_ALL="C"

# For ML server (RabbitMQ credentials)
export MUMMI_CRED="/lustre/orion/lrn005/proj-shared/spack/rabbitmq/"

# For flux job
export MUMMI_MPI_FLUX="simple"
# export MUMMI_MPI_FLUX="cray-pals"

# ------------------------------------------------------------------------------
# spack packages
# ------------------------------------------------------------------------------
export MUMMI_SPACK_VER="0.23"
export MUMMI_SPACK_ENV="$MUMMI_HOST"
#export MUMMI_SPACK_ROOT="/lustre/orion/lrn005/proj-shared/spack/mummi-spack/"

# This spack contains a more recent Flux and a specific 
#   maestrowf: https://github.com/lpottier/maestrowf/tree/feature/async-flux
#   it's a modified version of maestro that can submit job asynchronously and must faster than normal maestro
export MUMMI_SPACK_ROOT="/lustre/orion/lrn005/proj-shared/spack/mummi-dev-spack"

if [ ! -d "$MUMMI_SPACK_ROOT" ] ;then
  echo "(`hostname`: `date`) --> ERROR: [$MUMMI_SPACK_ROOT] is not a valid Spack installation." >&2
fi

echo "(`hostname`: `date`) --> Loading Spack installation ($MUMMI_SPACK_ROOT)"
source $MUMMI_SPACK_ROOT/$MUMMI_SPACK_VER/spack/share/spack/setup-env.sh

export MUMMI_SPACK_TARGET=$(spack arch --target)
echo "(`hostname`: `date`) --> Loading Spack target $MUMMI_SPACK_TARGET"

spack load flux-sched target=$MUMMI_SPACK_TARGET
if [ $? -ne 0 ]; then
  echo "Failed to load Spack for unknown reasons. aborting."
fi

spack load mummi-core target=$MUMMI_SPACK_TARGET
spack load mummi-ras target=$MUMMI_SPACK_TARGET

# spack load hwloc target=$MUMMI_SPACK_TARGET
# spack load py-maestrowf target=$MUMMI_SPACK_TARGET
# spack unload py-pip target=$MUMMI_SPACK_TARGET
spack load py-psutil target=$MUMMI_SPACK_TARGET
spack load py-pika target=$MUMMI_SPACK_TARGET
spack load py-keras target=$MUMMI_SPACK_TARGET
spack load py-keras-applications target=$MUMMI_SPACK_TARGET
spack load py-keras-preprocessing target=$MUMMI_SPACK_TARGET
# spack load gridsim2d target=$MUMMI_SPACK_TARGET
spack load py-torch target=$MUMMI_SPACK_TARGET
spack load py-lightning target=$MUMMI_SPACK_TARGET
spack load py-torch-sparse target=$MUMMI_SPACK_TARGET
spack load py-mpi4py target=$MUMMI_SPACK_TARGET

# Dirty solution to remove mummi-core/ras from PYTHONPATH
# spack unload does not work and unload all packages + dependencies
export PYTHONPATH=$(python3 $MUMMI_APP/setup/utils/clean-pythonpath.py)

# Specific to ROCm to avoid: RuntimeError: miopenStatusInternalError
# See: https://github.com/pytorch/pytorch/issues/60477
export MIOPEN_USER_DB_PATH="/tmp/mummi-miopen-cache"
export MIOPEN_CUSTOM_CACHE_DIR=${MIOPEN_USER_DB_PATH}
rm -rf ${MIOPEN_USER_DB_PATH}
mkdir -p ${MIOPEN_USER_DB_PATH}

spack load py-maestrowf target=$MUMMI_SPACK_TARGET
spack load py-parmed target=$MUMMI_SPACK_TARGET
spack load py-pyyaml target=$MUMMI_SPACK_TARGET

# ------------------------------------------------------------------------------
