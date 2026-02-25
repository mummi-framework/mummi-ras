#!/usr/bin/env bash
# Copyright 2021-2023 Lawrence Livermore National Security, LLC and other
#
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

usage="Usage: $(basename "$0") [#NODES] -- Script that bootstrap Flux on NNODES."

function version() {
  echo "$@" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }';
}

# Check if the allocation has the right size
# args:
#   - $1 : number of nodes requested for Flux
function check_main_allocation() {
  if [[ "$(flux getattr size)" -eq "$1" ]]; then
    echo "[$(date +'%m%d%Y-%T')@$(hostname)] Flux launch successful with $(flux getattr size) nodes"
  else
    echo "[$(date +'%m%d%Y-%T')@$(hostname)] Error: Requested nodes=$1 but Flux allocation size=$(flux getattr size)"
    exit 1
  fi
}

# Check if the 3 inputs are integers
function check_input_integers() {
  re='^[0-9]+$'
  if ! [[ $1 =~ $re && $2 =~ $re && $3 =~ $re ]] ; then
    echo "[$(date +'%m%d%Y-%T')@$(hostname)] Error: number of nodes is not an integer ($1, $2, $3)"
    exit 1
  fi
}

# Wait for a file to be created
#   - $1 : the file
#   - $2 : Max number of retry (one retry every 5 seconds)
function wait_for_file() {
  local FLUX_SERVER="$1"
  local EXIT_COUNTER=0
  local MAX_COUNTER="$2"
  while [ ! -f $FLUX_SERVER ]; do
    sleep 5s
    echo "[$(date +'%m%d%Y-%T')@$(hostname)] $FLUX_SERVER does not exist yet."
    exit_counter=$((EXIT_COUNTER + 1))
    if [ "$EXIT_COUNTER" -eq "$MAX_COUNTER" ]; then
      echo "[$(date +'%m%d%Y-%T')@$(hostname)] Timeout: Failed to find file (${FLUX_SERVER})."
      exit 1
    fi
  done
}

# ------------------------------------------------------------------------------
source /etc/profile.d/z00_lmod.sh
# ------------------------------------------------------------------------------
# spack packages
# ------------------------------------------------------------------------------
export MUMMI_SPACK_VER="0.19"
export MUMMI_SPACK_ROOT="/usr/workspace/mummiusr/mummi-spack/"

if [ ! -d "$MUMMI_SPACK_ROOT" ] ;then
  echo "[$(date +'%m%d%Y-%T')@$(hostname)] [$MUMMI_SPACK_ROOT] is not a valid Spack installation." >&2
  exit 1
fi
echo "[$(date +'%m%d%Y-%T')@$(hostname)] Loading Spack installation ($MUMMI_SPACK_ROOT)"
source $MUMMI_SPACK_ROOT/spack/$MUMMI_SPACK_VER/share/spack/setup-env.sh
export MUMMI_SPACK_TARGET=$(spack arch --target)
echo "[$(date +'%m%d%Y-%T')@$(hostname)] Loading Spack target $MUMMI_SPACK_TARGET"
spack load flux-sched target=$MUMMI_SPACK_TARGET

# the script needs the number of nodes for flux
FLUX_NODES="$1"
FLUX_SERVER="mummi-uri.log"
FLUX_LOG="mummi-flux.log"
# Flux-core Minimum version required
[[ -z ${MIN_VER_FLUX+z} ]] && MIN_VER_FLUX="0.45.0"

re='^[0-9]+$'
if ! [[ $FLUX_NODES =~ $re ]] ; then
  echo "[$(date +'%m%d%Y-%T')@$(hostname)] ERROR: '$FLUX_NODES' is not a number."
  echo $usage  
  exit 1
fi

echo "[$(date +'%m%d%Y-%T')@$(hostname)] Launching Flux with $FLUX_NODES nodes"

unset FLUX_URI
export LC_ALL="C"
export FLUX_F58_FORCE_ASCII=1
export FLUX_SSH="ssh"
# Cleanup from previous runs
rm -f $FLUX_SERVER $FLUX_LOG

if ! [[ -x "$(command -v flux)" ]]; then
  echo "[$(date +'%m%d%Y-%T')@$(hostname)] Error: flux is not installed."
  exit 1
fi
echo "[$(date +'%m%d%Y-%T')@$(hostname)] flux = $(which flux)"

flux_version=$(version $(flux version | awk '/^commands/ {print $2}'))
MIN_VER_FLUX_LONG=$(version ${MIN_VER_FLUX})
# We need to remove leading 0 because they are interpreted as octal numbers in bash
if [[ "${flux_version#00}" -lt "${MIN_VER_FLUX_LONG#00}" ]]; then
  echo "[$(date +'%m%d%Y-%T')@$(hostname)] Error Flux $(flux version | awk '/^commands/ {print $2}') is not supported.\
  Bootstrap requires flux>=${MIN_VER_FLUX}"
  exit 1
fi

echo "[$(date +'%m%d%Y-%T')@$(hostname)] flux version"
flux version

# We create a Flux wrapper around sleep on the fly to get the main Flux URI
FLUX_SLEEP_WRAPPER="./$(mktemp flux-wrapper.XXXX.sh)"
cat << 'EOF' > $FLUX_SLEEP_WRAPPER
#!/usr/bin/env bash
echo "ssh://$(hostname)$(flux getattr local-uri | sed -e 's!local://!!')" > "$1"
sleep inf
EOF
chmod u+x $FLUX_SLEEP_WRAPPER

MACHINE=$(echo $HOSTNAME | sed -e 's/[0-9]*$//')
if [[ "$MACHINE" == "lassen" ]] ; then
  # To use module command we must source this file
  # Those options are needed on IBM machines (CORAL)
  # Documented: https://flux-framework.readthedocs.io/en/latest/tutorials/lab/coral.html
  source /etc/profile.d/z00_lmod.sh
  module use /usr/tce/modulefiles/Core
  module use /usr/global/tools/flux/blueos_3_ppc64le_ib/modulefiles
  module load pmi-shim

  PMIX_MCA_gds="^ds12,ds21" \
    jsrun -a 1 -c ALL_CPUS -g ALL_GPUS -n ${FLUX_NODES} \
      --bind=none --smpiargs="-disable_gpu_hooks" \
      flux start -o,-S,log-filename=$FLUX_LOG -v $FLUX_SLEEP_WRAPPER $FLUX_SERVER &
elif [[ "$MACHINE" == "pascal" || "$MACHINE" == "ruby" ]] ; then
    srun -n ${FLUX_NODES} -N ${FLUX_NODES} --pty --mpi=none --mpibind=off \
      flux start -o,-S,log-filename=$FLUX_LOG -v $FLUX_SLEEP_WRAPPER $FLUX_SERVER &
else
  echo "[$(date +'%m%d%Y-%T')@$(hostname)] machine $MACHINE is not supported at the moment."
  exit 1
fi

echo ""
# now, wait for the flux info file
# we retry 20 times (one retry every 5 seconds)
wait_for_file $FLUX_SERVER 20
export FLUX_URI=$(cat $FLUX_SERVER)
echo "[$(date +'%m%d%Y-%T')@$(hostname)] Run: export FLUX_URI=$(cat $FLUX_SERVER)"
check_main_allocation ${FLUX_NODES}
