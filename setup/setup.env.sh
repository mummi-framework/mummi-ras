#!/usr/bin/env bash -i
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# set up script for mummi software stack
# ------------------------------------------------------------------------------
# if [[ "${BASH_SOURCE[0]}" == "${0}" || ! $ZSH_EVAL_CONTEXT =~ :file$  ]]; then
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "(`hostname`: `date`) --> ERROR: Please source this script!" >&2
  exit
fi

# ------------------------------------------------------------------------------
# identify the machine
# ------------------------------------------------------------------------------
host=$(hostname --long)

if [[ $host == *summit* ]]; then
  host=summit
elif [[ $host == lassen* ]]; then
  host=lassen
elif [[ $host == *frontier* ]]; then
  host=frontier
elif [[ $host == *ruby* || $host == *dane* ]]; then
  host=dane
elif [[ $host == *tuolumne* ]]; then
  host=tuolumne
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
if [ -z "$MUMMI_CONF" ]; then
  mummi_config=$HOME/.mummi/config.mummi.sh
  echo "(`hostname`: `date`) --> \$MUMMI_CONF was not set, defaulting to config in \$HOME directory: $mummi_config" >&2
else
  mummi_config="$MUMMI_CONF"
  echo "(`hostname`: `date`) --> using \$MUMMI_CONF=$MUMMI_CONF for config " >&2
fi


if [ ! -f $mummi_config ]; then
  echo "(`hostname`: `date`) --> ERROR: Could not find $mummi_config" >&2
  return
fi

echo "(`hostname`: `date`) --> Loading mummi configuration ($mummi_config)"
source $mummi_config

# We don't want interferences from local Spack
export SPACK_DISABLE_LOCAL_CONFIG=1
# Can make Spack faster
export SPACK_SKIP_MODULES=1

# ------------------------------------------------------------------------------
# set up custom config for this host (available on lassen and summit)
# ------------------------------------------------------------------------------
host_config=$MUMMI_APP/setup/envs/env.$host.sh
if [ ! -f $host_config ]; then
  echo "(`hostname`: `date`) --> ERROR: Could not find $host_config" >&2
  return
fi

# Prevent Spack to use local configs
export SPACK_DISABLE_LOCAL_CONFIG=true
echo "(`hostname`: `date`) --> Loading host ($host) environment ($host_config)"
source $host_config


# ------------------------------------------------------------------------------
# load the virtual environment
# ------------------------------------------------------------------------------
# todo: is there a better way of checking if virtualenv exists?
# if the virtual env exists
if [ -f "$MUMMI_ROOT/$MUMMI_VENV/bin/activate" ]; then
  echo "(`hostname`: `date`) --> Loading virtual environment (venv = $MUMMI_VENV)"
  source $MUMMI_ROOT/$MUMMI_VENV/bin/activate
  echo "python:    " $(command -v python)
  echo "pip:       "  $(command -v pip)
else
  echo "(`hostname`: `date`) --> Creating virtual environment (venv = $MUMMI_VENV)"
  python3 -m venv --system-site-packages $MUMMI_ROOT/$MUMMI_VENV
  source $MUMMI_ROOT/$MUMMI_VENV/bin/activate
  echo "python:    " $(command -v python)
  echo "pip:       "  $(command -v pip)
  $MUMMI_APP/setup/install/install_mummi_editable.sh
fi

#TODO: we should not have to add this path everytime..
# the venv should automatically do it for us!
# Fix by pottier1: one must use absolute paths when installing packages inside
# virtual envs
#python_version=`python -c 'import sys; print(str(sys.version_info[0])+"."+str(sys.version_info[1]))'`
#export PYTHONPATH=$VIRTUAL_ENV/lib/python$python_version/site-packages:$PYTHONPATH

export PATH=$MUMMI_AUTOBIND_PATH:$PATH

# ------------------------------------------------------------------------------
# check the paths. these should be pointing to the virtual env
echo "mummi_core:" `python -c "import mummi_core; print (mummi_core.__path__)"`
echo "mummi_ras:"  `python -c "import mummi_ras; print (mummi_ras.__path__)"`

# ------------------------------------------------------------------------------
# utility functions
# ------------------------------------------------------------------------------
#export MUMMI_UTILS=$MUMMI_APP/setup/utils
source $MUMMI_APP/setup/utils/misc_functions.sh

# ------------------------------------------------------------------------------
# Export explicit compiler binaries
# ------------------------------------------------------------------------------
echo "(`hostname`: `date`) --> Setting up compilers etc."
export CC=$(command -v mpicc)
export CXX=$(command -v mpic++)
export F90=mpif90
export F77=mpif77
export FC=mpif90

if false; then
  echo ''
  gcc --version | head -n1
  command -v gcc
  g++ --version | head -n1
  command -v g++
  python --version | head -n1
  command -v python
fi


# ------------------------------------------------------------------------------
# other necessary environment variables
# ------------------------------------------------------------------------------
export KERAS_BACKEND='theano'
export OMP_NUM_THREADS=4
umask 007

echo "(`hostname`: `date`) --> Setup to run MuMMI in ($MUMMI_ROOT)"
