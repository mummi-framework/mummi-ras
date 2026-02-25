#!/usr/bin/env bash
# ------------------------------------------------------------------------------

function version() {
  echo "$@" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }';
}

mummi_ucg_nnodes=$1

nucg_sims=20              # Number of UCG sims to launch
nb_threads_sim=1          # Number of OpenMP threads per UCG simulation
nb_threads_ucgserver=12   # Number of threads for the UCG server
nb_gpus_ucgserver=4       # Number of GPUs per node

re='^[0-9]+$'
if ! [[ $mummi_ucg_nnodes =~ $re ]] ; then
  echo "(`hostname`: `date`) --> ERROR: $mummi_ucg_nnodes is not a number"
  exit 1
fi

if [ $mummi_ucg_nnodes -eq 0 ]; then
  echo "(`hostname`: `date`) --> WARNING: Not Launching UCG server because mummi_ucg_nnodes = $mummi_ucg_nnodes"
  exit 1
fi

if [[ -z "${MUMMI_ROOT+x}" ]]; then
  echo "(`hostname`: `date`) --> MUMMI_ROOT is not set or is empty, please source your MuMMI config."
  exit 1
fi

export FLUX_SSH=ssh
flux_server="${MUMMI_ROOT}/flux/flux_server.info"

if [ -z ${FLUX_URI+x} ]; then
  echo "(`hostname`: `date`) --> FLUX_URI is not set or is empty, trying to recover it from flux/"
  if [[ -f ${flux_server} ]]; then
    export FLUX_URI=$(cat ${flux_server})
  else
    echo "(`hostname`: `date`) --> ERROR: Flux has not started, did not find ${flux_server} .. Aborting"
    exit 1
  fi
fi

# ------------------------------------------------------------------------------
echo "(`hostname`: `date`) --> Loading LAMMPS"
spack load lammps target=$MUMMI_SPACK_TARGET
if [ $? -ne 0 ]; then
  echo "Failed to load Spack and/or LAMMPS for unknown reasons. Aborting."
  exit 1
fi
# Needed to use LAMMPS Python package, This should be done by Spack automatically
# export PYTHONPATH="${PYTHONPATH}:$(spack location -i lammps)/lib/python3.9/site-packages/"

echo "(`hostname`: `date`) --> Loading GROMACS"
spack load gromacs@2019.6.classic target=$MUMMI_SPACK_TARGET
if [ $? -ne 0 ]; then
  echo "Failed to load GROMACS for unknown reasons. Aborting."
  exit 1
fi

ucg_exe=$(command -v lmp)
ucg_input="lammps.in"
ucg_traj="ucg.traj.xyz"
ucg_seed="seed.txt"

# ------------------------------------------------------------------------------
mkdir -p $MUMMI_ROOT/ucg
pushd $MUMMI_ROOT/ucg > /dev/null 2>&1

list_trajectories="$(pwd)/ucg_trajectories.txt"
rm -rf $list_trajectories
touch $list_trajectories

echo "(`hostname`: `date`) --> Launching ${nucg_sims} UCG simulations on ${mummi_ucg_nnodes} nodes with ${nb_threads_sim} threads"

#if in python env, we make sure we activate it
if [[ ! -z "${VIRTUAL_ENV}" ]]; then
  echo "(`hostname`: `date`) --> Sourcing $VIRTUAL_ENV/bin/activate"
  source $VIRTUAL_ENV/bin/activate
fi

max_cores_per_node=$(flux kvs get resource.R | flux R decode --include=0 --count=core)
ncores=$(($nucg_sims * $nb_threads_sim + $nb_threads_ucgserver))
echo "(`hostname`: `date`) --> total number of cores requested: $ncores (max cores/node: $max_cores_per_node)"

if [[ "$ncores" -gt "$max_cores_per_node" ]]; then
  echo "(`hostname`: `date`) --> ncores=$ncores cannot exceed $max_cores_per_node"
  echo "(`hostname`: `date`)     Too many UCG simulations, request more nodes to run the UCG server."
  exit 1
else
  # We the use the max number of cores/node
  ncores=$max_cores_per_node
fi
echo "(`hostname`: `date`) --> total number of cores granted: $ncores"

ucg_jid=$(
  flux mini batch --wrap --nodes=${mummi_ucg_nnodes} \
  --nslots=1 --cores-per-slot=${ncores} -g ${nb_gpus_ucgserver} \
  --output=ucgserver-{{id}}.log --job-name=ucgserver. \
  sleep inf)

ucg_flux_uri=$(flux uri --remote $ucg_jid)

FLUX_CMD="flux"

FLUX_VERSION=$(version $(flux version | awk '/^commands/ {print $2}'))
MIN_VER_FLUX=0.55 # Below 0.55 we can use flux mini but >= that we cannot
MIN_VER_FLUX_LONG=$(version ${MIN_VER_FLUX})
# Remove leading 0 because they are interpreted as octal numbers in bash
if [[ "${FLUX_VERSION#00}" -lt "${MIN_VER_FLUX_LONG#00}" ]]; then
  FLUX_CMD="flux mini"
fi

echo "(`hostname`: `date`) --> Launch batch for UCG $ucg_jid -> UCG FLUX_URI: $ucg_flux_uri"
# resource_info=$(flux proxy $ucg_jid flux hwloc info)
# echo "(`hostname`: `date`)     $resource_info"

echo "> UCG executable    : $ucg_exe"
echo "> UCG Command       : $ucg_cmd"

# --------------------- Launch all UCG simulations -----------------------------
for i in $(seq 1 $nucg_sims); do
  output="ucg-lammps-${i}"
  mkdir -p "${output}"
  pushd ${output} > /dev/null 2>&1
  cp $MUMMI_RESOURCES/ucg/* $(pwd)

  ts=$(date "+%Y%m%d%H%M")
  log_file=ucg.${i}.${ts}.log
  SEED=$RANDOM
  # Note that the part with -v must be before -in, otherwise LAMMPS will fail
  ucg_cmd="-v seed $SEED -in $ucg_input"

  echo "(`hostname`: `date`) --> Starting Ultra-Coarse Grained LAMMPS ${i}"

  if [[ -f "$ucg_traj" ]]; then
    echo "(`hostname`: `date`) --> Found trajectory file $ucg_traj for simulation ${i}. Restarting LAMMPS"
    ucg_input="lammps_restart.in"
    # We try to find the seed to restart the simulation with the same seed
    if [[ -f "$ucg_seed" ]]; then
      SEED=$(cat $ucg_seed)
    fi
    ucg_cmd="-v seed $SEED -in $ucg_input"
  fi

  echo "    > Working directory : $(pwd)"
  echo "    > Seed              : $SEED"

  # BUG: FLux proxy makes LAMMPS exit after 1s: flux proxy $ucg_jid 
  ucgsim_jid=$(FLUX_URI=$ucg_flux_uri ${FLUX_CMD} submit -n 1 -N 1 -c $nb_threads_sim \
    -o mpi=spectrum -o cpu-affinity=per-task --requires=rank:0 \
    --job-name=ucg.${i} --output=$log_file \
    sh -c "OMP_NUM_THREADS=$nb_threads_sim $ucg_exe $ucg_cmd")
  
  # Writing path for each trajectory so we can easily read them later (in UCG server for example)
  # This file is structured for K uCG simulations as follows:
  #
  #   name_ucgsim_1 seed_1 last_pos_read_1 traj_path_1
  #   ...
  #   name_ucgsim_K seed_K last_pos_read_K traj_path_K
  # 
  # last_pos_read_K represents the number of bytes read from file traj_path_K

  echo "ucg.${i} $SEED 0 $(pwd)/ucg.traj.xyz" >> $list_trajectories
  # Writing the seed we used for that simulation for easy access
  # TODO: replace with cut $list_trajectories
  echo "$SEED" > $ucg_seed

  popd > /dev/null 2>&1
done
popd > /dev/null 2>&1

# ----------------------------- Launch UCG Server ------------------------------
ucgserver_exe="python3 $MUMMI_APP/mummi_ras/scripts/run_ucgserver.py"
ucgserver_cmd=""

# We copy credentials to connect to RabbitMQ broker
echo "(`hostname`: `date`) --> Copying credentials for RabbitMQ server"
mkdir -p $MUMMI_ROOT/ucgserver
cp -v /usr/workspace/mummiusr/cz-ws-mummi/rabbitmq/rabbitmq-credentials.* $MUMMI_ROOT/ucgserver

pushd $MUMMI_ROOT/ucgserver > /dev/null 2>&1

echo "(`hostname`: `date`) --> Launching UCG Server with $mummi_ucg_nnodes nodes ($nb_gpus_ucgserver GPUs)"
echo "   Working directory      : $(pwd)"
echo "   UCG Server executable  : $ucgserver_exe"
echo "   UCG Server command     : $ucgserver_cmd"

# also write the command in the log file
echo "----> UCG Server Command: $ucgserver_exe $ucgserver_cmd"

# flux proxy $ucg_jid
ucgserver_jid=$(FLUX_URI=$ucg_flux_uri ${FLUX_CMD} submit -n 1 \
  -g $nb_gpus_ucgserver -N "$mummi_ucg_nnodes" --requires=rank:0 \
  --job-name=server. --output=server-{{id}}.log \
  sh -c "OMP_NUM_THREADS=$nb_threads_ucgserver $ucgserver_exe $ucgserver_cmd")

popd > /dev/null 2>&1

sleep 30s
# Safety check
nb_jobs_running=$(FLUX_URI=$ucg_flux_uri flux jobs -n | wc -l)
total_nb_jobs=$(($nucg_sims + 1))
if [[ $nb_jobs_running -ne $total_nb_jobs ]]; then
  echo "(`hostname`: `date`) --> ERROR: Number of jobs running is too low ($nb_jobs_running). $total_nb_jobs jobs should be running!"
  flux jobs -aR
  exit 1
fi

# ------------------------------------------------------------------------------