#!/usr/bin/env bash
# ------------------------------------------------------------------------------

function version() {
  echo "$@" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }';
}

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
# spack load py-pip target=$(spack arch --target)
# spack unload openssl target=$(spack arch --target)

# should we need this again?
#export FLUX_URI=`cat $MUMMI_ROOT/flux/flux_server.info`

echo "(`hostname`: `date`) --> Launching MuMMI workflow ($mummi_nnodes nodes)"
pushd $MUMMI_ROOT/workspace > /dev/null 2>&1

TMP_DIR="tmp_gen_chk"
# check if wfmanager.chk exists
if [ ! -f "wfmanager.chk" ]; then
    echo "(`hostname`: `date`) --> [Warning] File wfmanager.chk doesn't exist. Generate a new wfmanager.chk..."
    rm -rf $TMP_DIR
    mkdir $TMP_DIR
    python3 $MUMMI_APP/mummi_ras/scripts/update_wf_chkpt.py create $TMP_DIR
    cat $TMP_DIR/wfmanager_*.chk.yml
    cp $TMP_DIR/wfmanager_*.chk.yml wfmanager.chk
    rm -rf $TMP_DIR
    if [ ! -f "wfmanager.chk" ]; then
      echo "(`hostname`: `date`) --> [Error] Fail to generate wfmanager.chk. Exiting..."
      exit 1
    fi
fi

chk_file_size=`wc -c wfmanager.chk | awk '{print $1}'`
if [ $chk_file_size -lt 400 ]; then
  echo "(`hostname`: `date`) --> Checkpoint file wfmanager.chk is corrupted with size of $chk_file_size"
  chk_file=`ls -tr wfmanager.chk.bak.* |tail -n 1`
  cp $chk_file wfmanager.chk
  cp $chk_file wfmanager_test.chk
  echo "(`hostname`: `date`) --> Replace wfmanager.chk with backup $chk_file"
else
  echo "(`hostname`: `date`) --> Checkpoint file wfmanager.chk is fine size of $chk_file_size"
fi

sh $MUMMI_APP/setup/utils/backup_wspace.sh

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

popd > /dev/null 2>&1

# ------------------------------------------------------------------------------
