#!/usr/bin/env bash
# ------------------------------------------------------------------------------

# To use module command we must source this file
source /etc/profile.d/z00_lmod.sh

echo "(`hostname`: `date`) --> Loading pascal environment"

# pascal modules (Feb 2023)
# need to export this because we need to load this again later!
export MUMMI_MPI_MODULE="mvapich2/2.3"

module load gcc/8.3.1
module load cmake/3.23.1
module load cuda/11.6.1
module load $MUMMI_MPI_MODULE
#module load fftw/3.3.9

# ------------------------------------------------------------------------------
# Set host properties
export MUMMI_HOST="pascal"
export MUMMI_GROUP="mummiusr"

# Explicit paths
export MUMMI_AMBER_PATH="/usr/gapps/mummi/pascal/amber20/"
export MUMMI_AUTOBIND_PATH="/usr/gapps/mummi/bin"

# ------------------------------------------------------------------------------
# Flux shims
#source /etc/profile.d/z00_lmod.sh
export MUMMI_FLUX_MODULE_FILE="/usr/global/tools/flux/toss_3_x86_64_ib/modulefiles"
export MUMMI_FLUX_SHIM_MODULE="pmi-shim"
#export MUMMI_FLUX_MPI_MODULE="spectrum-mpi/2019.06.24-flux" # doesnt exist on Pascal

# locale messes up how flux job ids are reported
# https://github.com/flux-framework/flux-docs/issues/61
export LC_ALL="C"


# ------------------------------------------------------------------------------
# spack packages
# ------------------------------------------------------------------------------
#export MUMMI_SPACK_ROOT="/usr/gapps/mummi/spack"
#export MUMMI_SPACK_ROOT="/usr/workspace/mummidev/spack2"
export MUMMI_SPACK_VER="0.19"
export MUMMI_SPACK_ENV="$MUMMI_HOST"
export MUMMI_SPACK_ROOT="/usr/workspace/mummiusr/mummi-spack/"

if [ ! -d "$MUMMI_SPACK_ROOT" ] ;then
  echo "(`hostname`: `date`) --> ERROR: [$MUMMI_SPACK_ROOT] is not a valid Spack installation." >&2
  exit 1
fi

echo "(`hostname`: `date`) --> Loading spack installation ($MUMMI_SPACK_ROOT)"
source $MUMMI_SPACK_ROOT/spack/$MUMMI_SPACK_VER/share/spack/setup-env.sh
export MUMMI_SPACK_TARGET=$(spack arch --target)

spack load flux-sched target=$MUMMI_SPACK_TARGET
if [ $? -ne 0 ]; then
  echo "Failed to load Spack for unknown reasons. aborting."
  exit 1
fi
spack load mummi-ras target=$MUMMI_SPACK_TARGET
spack load py-maestrowf target=$MUMMI_SPACK_TARGET
spack load py-pip target=$MUMMI_SPACK_TARGET
spack load py-keras-applications target=$MUMMI_SPACK_TARGET
spack load py-keras-preprocessing target=$MUMMI_SPACK_TARGET
spack load gridsim2d target=$MUMMI_SPACK_TARGET
spack load py-torch target=$MUMMI_SPACK_TARGET
spack load py-pot target=$MUMMI_SPACK_TARGET

# unload mummi core and mummi ras
spack unload mummi-core mummi-ras target=$MUMMI_SPACK_TARGET

# ------------------------------------------------------------------------------
#export MUMMI_DDCMDCONVERTER_VERSION=ddcmdconverter/x7z6hqq

