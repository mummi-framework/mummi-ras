#!/usr/bin/env bash
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# set up script for mummi software stack
# ------------------------------------------------------------------------------
# TODO: this one is for bash only. handle zsh!
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "(`hostname`: `date`) --> ERROR: Please source this script!" >&2
  exit
fi

# ------------------------------------------------------------------------------
# identify the machine
# ------------------------------------------------------------------------------
host=`hostname --long`      # hostname

if [[ $host == *summit* ]]; then
  host=summit
#elif [[ $host == lassen* ]]; then
#  host=lassen
else
  echo "(`hostname`: `date`) --> ERROR: Unidentified host $host" >&2
  return
fi

# ------------------------------------------------------------------------------
if [ -n "$VIRTUAL_ENV" ] ; then
  echo "(`hostname`: `date`) --> WARNING: Deactivating the virtualenv ($VIRTUAL_ENV)"
  deactivate 2>/dev/null
fi

# ------------------------------------------------------------------------------
# find the path of this script
#if [ -n "$BASH" ] ;then
#   SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
#elif [ -n "$ZSH_NAME" ] ;then
#   SCRIPT_DIR="$(realpath `dirname ${0}`)"
#fi

# ------------------------------------------------------------------------------
# config script
mummi_config=$HOME/.mummi/config.mummi.sh
if [ ! -f $mummi_config ]; then
    echo "(`hostname`: `date`) --> ERROR: Could not find $mummi_config" >&2
    return
fi

echo "(`hostname`: `date`) --> Loading mummi configuration ($mummi_config)"
source $mummi_config


# host config
host_config=$MUMMI_APP/setup/envs/env.$host.sh
if [ ! -f $host_config ]; then
    echo "(`hostname`: `date`) --> ERROR: Could not find $host_config" >&2
    return
fi

echo "(`hostname`: `date`) --> Loading host environment ($host_config)"
source $host_config


# ------------------------------------------------------------------------------
# set up the spack installed for kras
# ------------------------------------------------------------------------------
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


spack load py-pip
# this is preventing proper functioning of ssh (on summit)
spack unload openssl

# ------------------------------------------------------------------------------
# load the virtual environment
# ------------------------------------------------------------------------------
# todo: is there a better way of checking if virtualenv exists?
# if the virtual env exists
if [ -d $VIRTUALENVWRAPPER_HOOK_DIR/$MUMMI_VENV ]; then
    echo "(`hostname`: `date`) --> Loading virtual environment (venv = $MUMMI_VENV)"
    workon $MUMMI_VENV
else
    echo "(`hostname`: `date`) --> Creating virtual environment (venv = $MUMMI_VENV)"
    mkvirtualenv $MUMMI_VENV
fi

#TODO: we should not have to add this path everytime..
# the venv should automatically do it for us!
export PYTHONPATH=$VIRTUAL_ENV/lib/python3.9/site-packages:$PYTHONPATH
export PATH=$MUMMI_AUTOBIND_PATH:$PATH

# ------------------------------------------------------------------------------
# utility functions
# ------------------------------------------------------------------------------
#export MUMMI_UTILS=$MUMMI_APP/setup/utils
source $MUMMI_APP/setup/utils/misc_functions.sh

# ------------------------------------------------------------------------------
# Export explicit compiler binaries
# ------------------------------------------------------------------------------
echo "(`hostname`: `date`) --> Setting up compilers etc."
export CC=`which mpicc`
export CXX=`which mpic++`
export F90=mpif90
export F77=mpif77
export FC=mpif90

if false; then
echo ''
gcc --version | head -n1
which gcc
g++ --version | head -n1
which g++
python --version | head -n1
which python
fi

# ------------------------------------------------------------------------------
# other necessary environment variables
# ------------------------------------------------------------------------------
export KERAS_BACKEND='theano'
export OMP_NUM_THREADS=4
umask 007

echo "(`hostname`: `date`) --> Setup to run MuMMI in ($MUMMI_ROOT)"

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
