#!/usr/bin/env bash
# ------------------------------------------------------------------------------

echo "(`hostname`: `date`) --> Loading lassen environment"

# lassen modules (Jun 2022)
# need to export this because we need to load this again later!
export MUMMI_MPI_MODULE="spectrum-mpi/2019.06.24"

module load gcc/8.3.1
module load cmake/3.20.2
module load cuda/10.1.105
module load $MUMMI_MPI_MODULE
module load fftw/3.3.8

# ------------------------------------------------------------------------------
# Set host properties
export MUMMI_HOST="lassen"
export MUMMI_GROUP="mummiusr"

# Explicit paths
export MUMMI_AMBER_PATH="/usr/gapps/mummi/amber18"
export MUMMI_AUTOBIND_PATH="/usr/gapps/mummi/bin"

# ------------------------------------------------------------------------------
# Flux shims
source /etc/profile.d/z00_lmod.sh
export MUMMI_FLUX_MODULE_FILE="/usr/global/tools/flux/blueos_3_ppc64le_ib/modulefiles"
export MUMMI_FLUX_SHIM_MODULE="pmi-shim"
export MUMMI_FLUX_MPI_MODULE="spectrum-mpi/2019.06.24-flux"

# locale messes up how flux job ids are reported
# https://github.com/flux-framework/flux-docs/issues/61
export LC_ALL="C"


# ------------------------------------------------------------------------------
# spack packages
# ------------------------------------------------------------------------------
export MUMMI_SPACK_ROOT="/usr/gapps/mummi/spack"

echo "(`hostname`: `date`) --> Loading spack environment ($MUMMI_SPACK_ROOT)"
source $MUMMI_SPACK_ROOT/share/spack/setup-env.sh

spack load -r mummi

spack unload mummi
spack unload hwloc~cuda

spack load -r hwloc+cuda
spack load -r py-maestrowf@1.1.9
spack load -r flux-sched@0.16.0

spack load -r py-pytaridx
spack load -r py-h5py
spack load -r py-keras

spack load swig@4.0.2/lk5op3c

spack load py-virtualenvwrapper
source `which virtualenvwrapper.sh`

# ------------------------------------------------------------------------------
