#!/usr/bin/env bash
# ------------------------------------------------------------------------------

# need to export this because we need to load this again later!
export MUMMI_MPI_MODULE="spectrum-mpi/2019.06.24"

# lassen modules (Jun 2022)
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

# Flux shims
source /etc/profile.d/z00_lmod.sh
export MUMMI_FLUX_MODULE_FILE="/usr/global/tools/flux/blueos_3_ppc64le_ib/modulefiles"
export MUMMI_FLUX_SHIM_MODULE="pmi-shim"
export MUMMI_FLUX_MPI_MODULE="spectrum-mpi/2019.06.24-flux"


# locale messes up how flux job ids are reported
# https://github.com/flux-framework/flux-docs/issues/61
export LC_ALL="C"

# ------------------------------------------------------------------------------
