#!/usr/bin/env bash
# ------------------------------------------------------------------------------

function version() {
  echo "$@" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }';
}

mummi_macro_nnodes=$1
re='^[0-9]+$'
if ! [[ $mummi_macro_nnodes =~ $re ]] ; then
  echo "(`hostname`: `date`) --> ERROR: $mummi_macro_nnodes is not a number" >&2
  exit 1
fi

if [ $mummi_macro_nnodes -eq 0 ]; then
  echo "(`hostname`: `date`) --> WARNING: Not Launching Macro because mummi_macro_nnodes = $mummi_macro_nnodes" >&2
  exit 1
fi

# ------------------------------------------------------------------------------
#spack load gridsim2dras target=$(spack arch --target)
macro_exe=`command -v gridsim2dras`

ts=`date "+%Y.%m.%d-%H.%M.%S"`
log_file=macro.${ts}.log
err_file=macro.${ts}.err

macro_out_freq=2500
macro_load_freq=100000
macro_niter=100000000
macro_cfg="-multi simlist.spec"
macro_params="iter0=last continue"

macro_oct=$MUMMI_HOST
macro_rdfpath=$MUMMI_ROOT/feedback-cg2macro

macro_cmd="$macro_cfg octave_path=$macro_oct rdf_path=$macro_rdfpath reload_freq=$macro_load_freq niter=$macro_niter $macro_params"

# ------------------------------------------------------------------------------
pushd $MUMMI_ROOT/macro

echo "(`hostname`: `date`) --> Launching Macro with $mummi_macro_nnodes nodes"
echo "   Working directory:  `pwd`"
echo "   Macro executable:   $macro_exe"
echo "   Macro Command:      $macro_cmd"

# also write the command in the log file
echo "----> Macro Command: $macro_exe $macro_cmd" >> $log_file

# TODO: check if we need to launch this via subbroker
#flux mini run -n 1 -N $MUMMI_MACRO_NNODES -c $((MUMMI_MACRO_NNODES * 24)) --output=$log_file --error=$err_file -o cpu-affinity=per-task -o mpi=spectrum \
# -o mpi=spectrum \

FLUX_CMD="flux"

FLUX_VERSION=$(version $(flux version | awk '/^commands/ {print $2}'))
MIN_VER_FLUX=0.55 # Below 0.55 we can use flux mini but >= that we cannot
MIN_VER_FLUX_LONG=$(version ${MIN_VER_FLUX})
# Remove leading 0 because they are interpreted as octal numbers in bash
if [[ "${FLUX_VERSION#00}" -lt "${MIN_VER_FLUX_LONG#00}" ]]; then
  FLUX_CMD="flux mini"
fi

${FLUX_CMD} run -n $((mummi_macro_nnodes * 24)) -N $mummi_macro_nnodes -c 1 \
 --output=$log_file --error=$err_file -o cpu-affinity=per-task \
 $macro_exe $macro_cmd &

popd
# ------------------------------------------------------------------------------
