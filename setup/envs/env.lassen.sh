#!/usr/bin/env bash
# ------------------------------------------------------------------------------

# To use module command we must source this file
source /etc/profile.d/z00_lmod.sh

echo "(`hostname`: `date`) --> Loading lassen environment"

# lassen modules (Jun 2022)
# need to export this because we need to load this again later!
#export MUMMI_MPI_MODULE="spectrum-mpi/2019.06.24"
# export MUMMI_MPI_MODULE="spectrum-mpi/2020.08.19"
export MUMMI_MPI_MODULE="spectrum-mpi/rolling-release"

# module load gcc/8.3.1
module load gcc/11.2.1 # For Spack 0.21
module load cmake/3.23.1
module load cuda/11.6.1
module load $MUMMI_MPI_MODULE
module load fftw/3.3.9

# ------------------------------------------------------------------------------
# Set host properties
export MUMMI_HOST="lassen"
export MUMMI_GROUP="mummiusr"

# Explicit paths
export MUMMI_AMBER_PATH="/usr/gapps/mummi/amber18"
# export MUMMI_AUTOBIND_PATH="/usr/gapps/mummi/bin"
#export MUMMI_AUTOBIND_PATH="/usr/gapps/mummi/bin:/usr/gapps/mummi/lassen/macro-gc/bin"
# export MUMMI_AUTOBIND_PATH="/usr/gapps/mummi/bin:/usr/gapps/mummi/lassen/macro-gc/pc3-iso/"
export MUMMI_AUTOBIND_PATH="/usr/gapps/mummi/bin:/usr/gapps/mummi/lassen/macro-gc/c1-iso:/usr/gapps/mummi/lassen/ddcmd/c1/"

# ------------------------------------------------------------------------------
# Flux shims
#source /etc/profile.d/z00_lmod.sh
export MUMMI_FLUX_MODULE_FILE="/usr/global/tools/flux/blueos_3_ppc64le_ib/modulefiles"
export MUMMI_FLUX_SHIM_MODULE="pmi-shim"
export MUMMI_FLUX_MPI_MODULE="spectrum-mpi/2019.06.24-flux"

# locale messes up how flux job ids are reported
# https://github.com/flux-framework/flux-docs/issues/61
export LC_ALL="C"

# For ML server (RabbitMQ credentials)
export MUMMI_CRED="/usr/workspace/mummiusr/cz-ws-mummi/rabbitmq/"

export MUMMI_MPI_FLUX="spectrum"

# ------------------------------------------------------------------------------
# spack packages
# ------------------------------------------------------------------------------
export MUMMI_SPACK_VER="0.21"
export MUMMI_SPACK_ENV="$MUMMI_HOST"
export MUMMI_SPACK_ROOT="/usr/workspace/mummiusr/mummi-spack/"

if [ ! -d "$MUMMI_SPACK_ROOT" ] ;then
  echo "(`hostname`: `date`) --> ERROR: [$MUMMI_SPACK_ROOT] is not a valid Spack installation." >&2
  return 1
fi

echo "(`hostname`: `date`) --> Loading Spack installation ($MUMMI_SPACK_ROOT)"
source $MUMMI_SPACK_ROOT/spack/$MUMMI_SPACK_VER/share/spack/setup-env.sh

export MUMMI_SPACK_TARGET=$(spack arch --target)
echo "(`hostname`: `date`) --> Loading Spack target $MUMMI_SPACK_TARGET"

spack load flux-sched target=$MUMMI_SPACK_TARGET

if [ $? -ne 0 ]; then
  echo "Failed to load Spack for unknown reasons. aborting."
  return 1
fi
spack load mummi-core target=$MUMMI_SPACK_TARGET
spack load mummi-ras target=$MUMMI_SPACK_TARGET

spack load hwloc target=$MUMMI_SPACK_TARGET
spack load py-psutil target=$MUMMI_SPACK_TARGET
spack load py-pika target=$MUMMI_SPACK_TARGET
spack load py-keras-applications target=$MUMMI_SPACK_TARGET
spack load py-keras-preprocessing target=$MUMMI_SPACK_TARGET
spack load py-torch target=$MUMMI_SPACK_TARGET
spack load py-lightning target=$MUMMI_SPACK_TARGET
spack load py-torch-sparse target=$MUMMI_SPACK_TARGET
spack load py-torch-scatter target=$MUMMI_SPACK_TARGET
spack load py-pot target=$MUMMI_SPACK_TARGET
spack load py-mpi4py target=$MUMMI_SPACK_TARGET

# Dirty solution to remove mummi-core/ras from PYTHONPATH
# spack unload does not work and unload all packages + dependencies
export PYTHONPATH=$(python3 $MUMMI_APP/setup/utils/clean-pythonpath.py)

# ------------------------------------------------------------------------------
