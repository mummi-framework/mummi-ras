#!/usr/bin/env bash
# ------------------------------------------------------------------------------

function version() {
  echo "$@" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }';
}

mummi_mlserver_nnodes=$1
re='^[0-9]+$'
if ! [[ $mummi_mlserver_nnodes =~ $re ]] ; then
  echo "(`hostname`: `date`) --> ERROR: $mummi_mlserver_nnodes is not a number" >&2
  exit 1
fi

if [ $mummi_mlserver_nnodes -eq 0 ]; then
  echo "(`hostname`: `date`) --> WARNING: Not Launching ML server because mummi_mlserver_nnodes = $mummi_mlserver_nnodes" >&2
  exit 1
fi

# ------------------------------------------------------------------------------
ts=`date "+%Y.%m.%d-%H.%M.%S"`
log_file=mlserver.${ts}.log
err_file=mlserver.${ts}.err

mlserver_exe="python3 $MUMMI_APP/mummi_ras/scripts/run_mlserver.py"
mlserver_cmd=""

mkdir -p $MUMMI_ROOT/mlserver

if [[ -z "${MUMMI_CRED}" ]]; then
  echo "(`hostname`: `date`) --> MUMMI_CRED is empty \"${MUMMI_CRED}\". Cannot find credentials for RabbitMQ"
  exit 1
fi

cp ${MUMMI_CRED}/rabbitmq-credentials.* $MUMMI_ROOT/mlserver

# ------------------------------------------------------------------------------
pushd $MUMMI_ROOT/mlserver

# HERE: the number of cores NUM_THREADS + NUM_THREADS from launch_workflow.sh must be less than
# total number of cores on one node
NUM_THREADS=52 # BE CAREFUL TO NOT ALLOCATE MORE CORES THAN AVAILABLE

echo "(`hostname`: `date`) --> Launching ML Server with $mummi_mlserver_nnodes nodes"
echo "   Working directory:      `pwd`"
echo "   ML Server executable:   $mlserver_exe"

spack load gromacs@2024.01~mpi~mummi-c2 target=$MUMMI_SPACK_TARGET

#if in python env, we make sure we activate it
if [[ ! -z "${VIRTUAL_ENV}" ]]; then
  echo "(`hostname`: `date`) --> Sourcing $VIRTUAL_ENV/bin/activate"
  source $VIRTUAL_ENV/bin/activate
fi

FLUX_CMD="flux"

FLUX_VERSION=$(version $(flux version | awk '/^commands/ {print $2}'))
MIN_VER_FLUX=0.55 # Below 0.55 we can use flux mini but >= that we cannot
MIN_VER_FLUX_LONG=$(version ${MIN_VER_FLUX})
# Remove leading 0 because they are interpreted as octal numbers in bash
if [[ "${FLUX_VERSION#00}" -lt "${MIN_VER_FLUX_LONG#00}" ]]; then
  FLUX_CMD="flux mini"
  jobid=$(${FLUX_CMD} submit -N "$mummi_mlserver_nnodes" -n 1 -c $NUM_THREADS \
    --env=OMP_NUM_THREADS=$NUM_THREADS \
    --job-name=mlserver. \
    --output=mlserver-{{id}}.log \
    sh -c "$mlserver_exe $mlserver_cmd")
else
  jobid=$(${FLUX_CMD} submit -N "$mummi_mlserver_nnodes" -n 1 -c $NUM_THREADS \
    --env=OMP_NUM_THREADS=$NUM_THREADS \
    --requires=rank:0 \
    --job-name=mlserver. \
    --output=mlserver-{{id}}.log \
    sh -c "$mlserver_exe $mlserver_cmd")
fi

echo "(`hostname`: `date`) --> Launched ML Server (job ID: $jobid)"

popd > /dev/null 2>&1
# ------------------------------------------------------------------------------
