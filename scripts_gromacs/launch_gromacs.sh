#!/usr/bin/env bash

MUMMI_CONF=~/.mummi/config.mummi.sh
NODES="$1"
TIMEOUT="$2"
#QUEUE="pdebug"
#QUEUE="pbatch"
QUEUE="expedite"
DEP="$3"

CDIR="sims-cg"
SIMTIME="12" #In hours

EXE="./run_cgsims_gromacs.sh $CDIR $NODES $SIMTIME"

source $MUMMI_CONF
if [ -z "${DEP}" ]; then
    bsub -G cancer -nnodes $NODES -W $TIMEOUT -q $QUEUE $EXE
else    
    bsub -G cancer -nnodes $NODES -W $TIMEOUT -q $QUEUE -w "ended($DEP)" $EXE
fi
