#!/usr/bin/env bash
# ------------------------------------------------------------------------------

# the script needs the number of mummi nodes
mummi_nnodes=$1
re='^[0-9]+$'
if ! [[ $mummi_nnodes =~ $re ]] ; then
   echo "(`hostname`: `date`) --> ERROR: Nodes is not a number" >&2
   exit 1
fi

if [ -z "$MUMMI_ROOT" ] ;then
  echo "(`hostname`: `date`) --> ERROR: Not launching workflow because MUMMI_ROOT is not set" >&2
  exit 1
fi

# need access to mummi. but paths getting overridden? this is a temporary fix!
spack load py-pip
spack unload openssl

# should we need this again?
#export FLUX_URI=`cat $MUMMI_ROOT/flux/flux_server.info`

echo "(`hostname`: `date`) --> Launching MuMMI workflow ($mummi_nnodes nodes)"
pushd $MUMMI_ROOT/workspace > /dev/null 2>&1

## from pilot/mummi-ras (old) repo!
## flux mini submit -n 1 -c 24 -N 1 --env=OMP_NUM_THREADS=4 -o mpi=spectrum $SCRIPT_DIR/wfmanager.sh

flux mini submit -n 1 -c 24 -N 1 --env=OMP_NUM_THREADS=4 \
	sh -c "autobind-24 python3 $MUMMI_APP/mummi_ras/scripts/run_wfmanager.py"

#"python $MUMMI_APP/mummi_ras/scripts/run_wfmanager.py >> wfmanager.log 2>&1"
#autobind-24 mummi_workflow >> wfmanager.log 2>&1

popd > /dev/null 2>&1

# ------------------------------------------------------------------------------
