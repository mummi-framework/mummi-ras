#!/usr/bin/env bash

ts=`date "+%Y.%m.%d-%H.%M.%S"`
log_file=macro.${ts}.log
err_file=macro.${ts}.err
#macro_bin="/autofs/nccs-svm1_proj/lrn005/gridcorr2d/gridsim2dras"
macro_bin=`which gridsim2dras`

MACRO_OUT_FREQ=2500
MACRO_LOAD_FREQ=100000
MACRO_NITER=100000000
MACRO_CFG="-multi simlist.spec"
MACRO_OCT=$MUMMI_HOST
PARAMS="iter0=last continue"

# node_counts.sh removed by JM (11/09/21)
# source $MUMMI_APP/setup/node_counts.sh
SETUP_CONFIG=$HOME/.mummi/config.mummi.sh
source $SETUP_CONFIG

cd $MUMMI_ROOT/macro

echo "------------------------------------------------"
echo "MUMMI_MACRO_NNODES: $MUMMI_MACRO_NNODES"
echo "Macro execurable:   $macro_bin"
echo "Working directory:  `pwd`"
echo "Octave Path:        `echo $MUMMI_HOST`"
echo "Macro Params:       $PARAMS"
echo "------------------------------------------------"

# TODO: check if we need to launch this via subbroker
#flux mini run -n 1 -N $MUMMI_MACRO_NNODES -c $((MUMMI_MACRO_NNODES * 24)) --output=$log_file --error=$err_file -o cpu-affinity=per-task -o mpi=spectrum \
MACRO_CMD="$macro_bin $MACRO_CFG octave_path=$MACRO_OCT rdf_path=$MUMMI_ROOT/feedback-cg2macro reload_freq=$MACRO_LOAD_FREQ niter=$MACRO_NITER $PARAMS"
echo "Macro Command: $MACRO_CMD"
echo "Macro Command: $MACRO_CMD" >> $log_file

flux mini run -n $((MUMMI_MACRO_NNODES * 24)) -N $MUMMI_MACRO_NNODES -c 1 --output=$log_file --error=$err_file -o cpu-affinity=per-task -o mpi=spectrum \
$MACRO_CMD &
