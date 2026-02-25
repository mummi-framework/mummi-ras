#!/usr/bin/env bash
# Script to run GROMACS simulations 

export PATH=$PATH:/usr/gapps/mummi/bin # to get the executable jsrun_get_gpu
get_pinoffset() {
    data=$(OMP_NUM_THREADS=1 jsrun -n1 -r1 -a1 -c1 -g1 /usr/gapps/mummi/bin/jsrun_get_gpu | awk -F'=' '{print $3}' | awk 'NF')
    # jsrun_get_gpu return if we reqest 3 GPUs
    #   OMP_thread=01:GPU_id=G1 G2 G3
    #   OMP_thread=XX:GPU_id=G1 G2 G3
    #   etc for X all threads requested
    # Print first line and first value:
    #   - all line should have the same GPU IDs
    #   - and there should only be only one GPUs
    gpu_id=$(echo $data | awk 'NR==1{print $1}') 
    local offset=8
    case $gpu_id in
        0) offset=8;;
        1) offset=48;;
        2) offset=96;;
        3) offset=136;;
        *) echo "GPUs ID unknown error.";; 
    esac
    echo $offset
}

cdir=$1
nodecount=$2
simtime=$3

if [[ $# != 3 ]]; then
  echo "error 3 arguments expected: $0 cdir nodecount simtime"
  exit 1
fi

if [[ ! -e $cdir ]]; then
  echo "error: directory $cdir does not exist."
  exit 1
fi

maxsims=$(($nodecount * 4))

echo "Running CG simulations with GROMACS"
echo " sim dir: $cdir"
echo " Maximum number of sims possible: $maxsims"
echo " Each sim run for $simtime hours"

#source /etc/profile.d/z00_lmod.sh
source /usr/workspace/mummiusr/mummi-spack/spack/0.19/share/spack/setup-env.sh
export MUMMI_SPACK_TARGET=$(spack arch --target)
echo " Loading Gromacs"
spack load gromacs@2022.06.gpu target=$MUMMI_SPACK_TARGET
export OMP_NUM_THREADS=8

if ! [ -x "$(command -v gmx)" ]; then
  echo "error: gmx is not installed."
  exit 1
fi

echo " Loaded Gromacs: $(which gmx)"

currentDir=$(pwd)
ts=$(date +"%Y%m%d%H%M")
# file where we store the list of sims-cg to process
sim_run_file="${currentDir}/gmx_to_run-${ts}.txt"
rm -rf $sim_run_file
touch $sim_run_file

cd $cdir

for cSim in structure_iter00_*
do 
  if [ ! -f $cSim/createsims_success ]; then
    echo "$cSim - NOT createsims_success"
  elif [ -f $cSim/cg_stop ]; then
    echo "$cSim - NOT cg_stop"
  elif [ -f $cSim/cg_success ]; then
    echo "$cSim - NOT cg_success"
  elif [ -f $cSim/cg_failure ]; then
    echo "$cSim - NOT cg_failure"
  else
    echo "$cSim - YES extend sim"
    echo "$(realpath $cSim)" >> $sim_run_file
  fi
done

echo "Found $(cat $sim_run_file | wc -l) sims to run"

i=0
for sim in $(cat $sim_run_file); do
    cd $sim
    # We get the offset for Gromacs that depends on the GPU IDs
    offset=$(get_pinoffset)
    echo "Run $sim with pinoffset $offset maxh $simtime "
    sleep 1 # we sleep 1 seconds for safety between the jsrun call done in get_pinoffset and the one we are doing now
    jsrun -n1 -r1 -a1 -c1 -g1 --smpiargs=none gmx mdrun -cpi -nb gpu -ntmpi 1 -ntomp $OMP_NUM_THREADS -pin on -pinstride 2 -pinoffset "$offset" -maxh $simtime >> mdrun.log 2>&1 &
    echo "Submitted $sim for $simtime hours"
    cd $currentDir
    i=$((i+1))
    if [[ "$i" -eq "$maxsims" ]]; then
        echo "We have reached the max number of simulations to schedule."
        break
    fi
done

jslist -r

jswait all
echo "All sims done"

