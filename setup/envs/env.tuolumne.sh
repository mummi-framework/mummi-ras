#!/usr/bin/env bash
# ------------------------------------------------------------------------------

# To use module command we must source this file
source /etc/profile.d/z00_lmod.sh

echo "(`hostname`: `date`) --> Loading tuolumne environment"

export MUMMI_MPI_MODULE=cray-mpich/8.1.31

#module load gcc/12.2.1
module load cmake/3.23.1
module load $MUMMI_MPI_MODULE
module load cray-fftw/3.3.10.9
module load rocm/6.3.1


# ------------------------------------------------------------------------------
# Set host properties
export MUMMI_HOST="tuolumne"
export MUMMI_GROUP="mummiusr"

# Explicit paths
export MUMMI_AMBER_PATH="/usr/gapps/mummi/amber18"
# export MUMMI_AUTOBIND_PATH="/usr/gapps/mummi/bin"
#export MUMMI_AUTOBIND_PATH="/usr/gapps/mummi/bin:/usr/gapps/mummi/lassen/macro-gc/bin"
# export MUMMI_AUTOBIND_PATH="/usr/gapps/mummi/bin:/usr/gapps/mummi/lassen/macro-gc/pc3-iso/"
export MUMMI_AUTOBIND_PATH="/usr/gapps/mummi/bin:/usr/gapps/mummi/lassen/macro-gc/c1-iso:/usr/gapps/mummi/lassen/ddcmd/c1/"

# locale messes up how flux job ids are reported
# https://github.com/flux-framework/flux-docs/issues/61
export LC_ALL="C"

# For ML server (RabbitMQ credentials)
export MUMMI_CRED="/usr/workspace/mummiusr/cz-ws-mummi/rabbitmq/"

export MUMMI_MPI_FLUX="simple" #???

# ------------------------------------------------------------------------------
# spack packages
# ------------------------------------------------------------------------------
export MUMMI_SPACK_VER="0.23"
export MUMMI_SPACK_ENV="$MUMMI_HOST"
export MUMMI_SPACK_ROOT="/usr/workspace/mummiusr/mummi-spack/"

if [ ! -d "$MUMMI_SPACK_ROOT" ] ;then
  echo "(`hostname`: `date`) --> ERROR: [$MUMMI_SPACK_ROOT] is not a valid Spack installation." >&2
  return 1
fi

echo "(`hostname`: `date`) --> Loading Spack installation ($MUMMI_SPACK_ROOT)"
source $MUMMI_SPACK_ROOT/spack/$MUMMI_SPACK_VER/share/spack/setup-env.sh

export MUMMI_SPACK_TARGET="linux-rhel8-zen3"
echo "(`hostname`: `date`) --> Loading Spack arch $MUMMI_SPACK_TARGET"

spack load flux-sched arch=$MUMMI_SPACK_TARGET

if [ $? -ne 0 ]; then
  echo "Failed to load Spack for unknown reasons. aborting."
  return 1
fi
spack load mummi-core arch=$MUMMI_SPACK_TARGET
spack load mummi-ras arch=$MUMMI_SPACK_TARGET

# Dirty solution to remove mummi-core/ras from PYTHONPATH
# spack unload does not work and unload all packages + dependencies
export PYTHONPATH=$(python3 $MUMMI_APP/setup/utils/clean-pythonpath.py)

spack load hwloc arch=$MUMMI_SPACK_TARGET
spack load py-psutil arch=$MUMMI_SPACK_TARGET
spack load py-pika arch=$MUMMI_SPACK_TARGET
spack load py-keras-applications arch=$MUMMI_SPACK_TARGET
spack load py-keras-preprocessing arch=$MUMMI_SPACK_TARGET
spack load py-torch arch=$MUMMI_SPACK_TARGET
spack load py-lightning arch=$MUMMI_SPACK_TARGET
spack load py-torch-sparse arch=$MUMMI_SPACK_TARGET
spack load py-pot arch=$MUMMI_SPACK_TARGET
spack load py-mpi4py arch=$MUMMI_SPACK_TARGET

# ------------------------------------------------------------------------------
