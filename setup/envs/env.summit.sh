#!/usr/bin/env bash
# ------------------------------------------------------------------------------

# need to export this because we need to load this again later!
export MUMMI_MPI_MODULE="spectrum-mpi/10.4.0.3-20210112"

# summit modules (May 2022)
module load gcc/9.1.0          # > /dev/null 2>&1
module load cuda/11.0.3        # > /dev/null 2>&1
module load cmake/3.23.1       # > /dev/null 2>&1
module load rust/1.51.0        # > /dev/null 2>&1
module load sqlite/3.36.0      # > /dev/null 2>&1
module load $MUMMI_MPI_MODULE  # > /dev/null 2>&1

# ------------------------------------------------------------------------------
# Set host properties
export MUMMI_HOST="summit"
export MUMMI_GROUP="lrn005"

# Explicit paths
export MUMMI_SPACK_ROOT="/autofs/nccs-svm1_proj/lrn005/spack2022"
export MUMMI_AMBER_PATH="/autofs/nccs-svm1_proj/lrn005/amber18"
export MUMMI_AUTOBIND_PATH="/autofs/nccs-svm1_proj/lrn005/bin"

# Flux shims
export MUMMI_FLUX_MODULE_FILE="/sw/summit/modulefiles/ums/gen007flux/Core"
export MUMMI_FLUX_SHIM_MODULE="pmi-shim"
export MUMMI_FLUX_MPI_MODULE="spectrum-mpi/10.3.0.1-20190611-flux"

# locale messes up how flux job ids are reported
# https://github.com/flux-framework/flux-docs/issues/61
export LC_ALL="C"

# ------------------------------------------------------------------------------
# two packages (theano and octave) create temp files in home directory
# but summit doesnt allow accessing the home directory
# so we hack them and use a symlink to a write-able directory
gpath="$MUMMI_GRP_ROOT/mummi_run_temp"
if [ ! -d ~/.theano ]; then
   tname="theano-$USER"
   mkdir -p $gpath/theano/$tname
   ln -s $gpath/theano/$tname ~/.theano 2>/dev/null
fi

if [ ! -f ~/.octave_hist ]; then
   tname="octave_hist-$USER"
   oct_hist=$gpath/octave_hist/$tname
   mkdir -p `dirname $oct_hist`
   touch $oct_hist
   ln -s $oct_hist ~/.octave_hist 2>/dev/null
fi

# ------------------------------------------------------------------------------
