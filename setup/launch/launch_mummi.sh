#!/usr/bin/env bash

# These wars need to be exporetd beforand (x3 first by calling script, LSB_JOBIN by LSF)
echo 'Running splash launcher'
echo '    MUMMI_ROOT            :' $MUMMI_ROOT
echo '    MUMMI_APP             :' $MUMMI_APP
echo '    MUMMI_CORE            :' $MUMMI_CORE
echo '    MUMMI_RESTART         :' $MUMMI_RESTART
echo '    LSB_JOBID             :' $LSB_JOBID
echo ''

if [ -n "$BASH" ] ;then
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
elif [ -n "$ZSH_NAME" ] ;then
    SCRIPT_DIR="$(realpath `dirname ${0}`)"
fi

if [ -z "$MUMMI_NNODES" ] ;then
    echo "ERROR: MUMMI_NNODES is not set";
    exit 1;
fi

export FLUX_URI=`cat $MUMMI_ROOT/flux/flux.info`
flux mini submit -n 1 -c 24 -N 1 --env=OMP_NUM_THREADS=4 -o mpi=spectrum $SCRIPT_DIR/wfmanager.sh
echo 'Returning from launch_script'
