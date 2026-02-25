#!/usr/bin/env bash
# ------------------------------------------------------------------------------
if [ -f "~/.bash_profile" ]; then
  source ~/.bash_profile
fi
USRNAME=`whoami`

echo ""
echo "(`hostname`: `date`) --> Launching MuMMI workflow manager only for (`whoami`)"

# ------------------------------------------------------------------------------
# input parameters
# ------------------------------------------------------------------------------
#mummi_nnodes=$1
mummi_nnodes=$1
re='^[0-9]+$'
if ! [[ $mummi_nnodes =~ $re ]] ; then
  echo "(`hostname`: `date`) --> ERROR: mummi_nnodes=$mummi_nnodes is not a number" >&2
  exit 1
fi

# ------------------------------------------------------------------------------
# load the mummi config from ($HOME/.mummi/config.mummi.sh)
# ------------------------------------------------------------------------------
mummi_config=$HOME/.mummi/config.mummi.sh
if [ ! -f $mummi_config ]; then
  echo "(`hostname`: `date`) --> ERROR: Could not find ($mummi_config)" >&2
  exit 1
fi
source $mummi_config

# step 4: launch workflow
export MUMMI_NNODES=$mummi_nnodes
export MUMMI_REDIS_NNODES=0
export MUMMI_MACRO_NODES=0
export MUMMI_MLSERVER_NNODES=1
export MUMMI_UCGSERVER_NNODES=0
export MUMMI_FB_NNODES=1

# ------------------------------------------------------------------------------
# load the mummi environment
# ------------------------------------------------------------------------------
source $MUMMI_APP/setup/setup.env.sh
source $MUMMI_ROOT/flux/flux_client.env

# Just a safety measure
wfm=$(flux jobs --count=80000 | grep wfmanager)
if [[ "$?" -eq 0 ]]; then
  echo "(`hostname`: `date`) --> A workflow manager already runs. Abort"
  echo "(`hostname`: `date`) $wfm"
  exit 1
fi

# need access to mummi. but paths getting overridden? this is a temporary fix!
spack load py-pip target=$(spack arch --target)
# spack unload openssl target=$(spack arch --target)

echo "(`hostname`: `date`) --> Launching MuMMI workflow ($mummi_nnodes nodes)"
pushd $MUMMI_ROOT/workspace > /dev/null 2>&1

#if in python env, we make sure we activate it
if [[ ! -z "${VIRTUAL_ENV}" ]]; then
  echo "(`hostname`: `date`) --> Sourcing $VIRTUAL_ENV/bin/activate"
  source $VIRTUAL_ENV/bin/activate
fi


NUM_THREADS=4
FLUX_CMD="flux"

FLUX_VERSION=$(version $(flux version | awk '/^commands/ {print $2}'))
MIN_VER_FLUX=0.55 # Below 0.55 we can use flux mini but >= that we cannot
MIN_VER_FLUX_LONG=$(version ${MIN_VER_FLUX})
# Remove leading 0 because they are interpreted as octal numbers in bash
if [[ "${FLUX_VERSION#00}" -lt "${MIN_VER_FLUX_LONG#00}" ]]; then
  FLUX_CMD="flux mini"
  jobid=$(${FLUX_CMD} submit -n 1 -c ${NUM_THREADS} -o mpi=${MUMMI_MPI_FLUX} \
    --env=OMP_NUM_THREADS=${NUM_THREADS} \
    --job-name=wfmanager. \
    --output=wfmanager-{{id}}.err \
    sh -c "python3 $MUMMI_APP/mummi_ras/scripts/run_wfmanager.py")
else
  jobid=$(${FLUX_CMD} submit -n 1 -c ${NUM_THREADS} -o mpi=${MUMMI_MPI_FLUX} \
    --env=OMP_NUM_THREADS=${NUM_THREADS} \
    --requires=rank:0 \
    --job-name=wfmanager. \
    --output=wfmanager-{{id}}.err \
    sh -c "python3 $MUMMI_APP/mummi_ras/scripts/run_wfmanager.py")
fi

echo "(`hostname`: `date`) --> Launched workflow on ${NUM_THREADS} cores (job ID: $jobid / MPI=${MUMMI_MPI_FLUX})"

