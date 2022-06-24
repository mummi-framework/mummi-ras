#!/usr/bin/env bash
# ------------------------------------------------------------------------------

echo "(`hostname`: `date`) --> Loading summit environment"

# summit modules (May 2022)
# need to export this because we need to load this again later!
export MUMMI_MPI_MODULE="spectrum-mpi/10.4.0.3-20210112"

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
export MUMMI_AMBER_PATH="/autofs/nccs-svm1_proj/lrn005/amber18"
export MUMMI_AUTOBIND_PATH="/autofs/nccs-svm1_proj/lrn005/bin"

# ------------------------------------------------------------------------------
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
# spack packages
# ------------------------------------------------------------------------------
export MUMMI_SPACK_ROOT="/autofs/nccs-svm1_proj/lrn005/spack2022"

echo "(`hostname`: `date`) --> Loading spack environment ($MUMMI_SPACK_ROOT)"
source $MUMMI_SPACK_ROOT/share/spack/setup-env.sh

spack load py-virtualenvwrapper
source `which virtualenvwrapper.sh`

spack load py-pytest
spack load py-psutil
spack load py-filelock
spack load py-pyyaml
spack load py-numpy
spack load py-pytaridx
spack load py-redis
spack load py-maestrowf
spack load flux-sched
spack load py-cryptography@36.0.1
spack load py-scipy
spack load py-parmed
spack load py-matplotlib
spack load py-mdanalysis-mummi
spack load faiss
spack load py-keras
spack load py-h5py@2.8.0~mpi
spack load swig

spack load py-pip
# this is preventing proper functioning of ssh (on summit)
spack unload openssl

# ------------------------------------------------------------------------------
