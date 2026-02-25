#!/usr/bin/env bash
# ------------------------------------------------------------------------------

function version() {
  echo "$@" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }';
}

# the script needs the number of mummi nodes
mummi_cga_nnodes=$1
re='^[0-9]+$'
if ! [[ $mummi_cga_nnodes =~ $re ]] ; then
   echo "(`hostname`: `date`) --> ERROR: Nodes is not a number" >&2
   exit 1
fi

if [ $mummi_cga_nnodes -eq 0 ]; then
  echo "(`hostname`: `date`) --> WARNING: Not Launching CG analysis because node = 0" >&2
  exit 1
fi

if [ -z "$MUMMI_ROOT" ] ;then
  echo "(`hostname`: `date`) --> ERROR: Not launching workflow because MUMMI_ROOT is not set" >&2
  exit 1
fi

echo "(`hostname`: `date`) --> Launching CG analysis ($mummi_cga_nnodes nodes)"
pushd $MUMMI_ROOT/workspace > /dev/null 2>&1

# load mummi_environment (already loaded by launch_mummi.sh)
# source $HOME/.mummi/config.mummi.sh
# source $MUMMI_APP/setup/setup.env.sh

#if in python env, we make sure we activate it
if [[ ! -z "${VIRTUAL_ENV}" ]]; then
  echo "(`hostname`: `date`) --> Sourcing $VIRTUAL_ENV/bin/activate"
  source $VIRTUAL_ENV/bin/activate
fi

NCORES=80

FLUX_VERSION=$(version $(flux version | awk '/^commands/ {print $2}'))
MIN_VER_FLUX=0.55 # Below 0.55 we can use flux mini but >= that we cannot
MIN_VER_FLUX_LONG=$(version ${MIN_VER_FLUX})
# Remove leading 0 because they are interpreted as octal numbers in bash
if [[ "${FLUX_VERSION#00}" -lt "${MIN_VER_FLUX_LONG#00}" ]]; then
  FLUX_CMD="flux mini"
fi

jobid=$(${FLUX_CMD} submit  -N ${mummi_cga_nnodes} -n ${NCORES} -c 1 --job-name analysis_cg \
            -o mpi=spectrum -o cpu-affinity=per-task -o gpu-affinity=per-task \
            python3 $MUMMI_APP/mummi_ras/online/cg/cg_analysis_executor.py --simname "" \
            --noninteractive \
            --path "" \
            --fstype "mummi" \
            --fbio "mummi" \
            --frameProcessBatchSize 100 \
            --siminputs "" \
            --simbin "" \
            --fcount 0 \
            --useGromacs)

echo "(`hostname`: `date`) --> Launched CG analysis (job ID: $jobid)"

popd > /dev/null 2>&1

# ------------------------------------------------------------------------------
