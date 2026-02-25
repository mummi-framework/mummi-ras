#!/bin/bash
#SBATCH -N 1
#SBATCH -t 1:00:00
#SBATCH -p batch
#SBATCH -A lrn005
#SBATCH -J mummi-c3-gen
#SBATCH -o %x-%j.out


module load rocm/6.3.1
module load craype-accel-amd-gfx90a

# Environment setup
# TODO: replace with your own config.mummi.sh
source /ccs/home/lpottier/mummi-framework-c3/config.mummi.sh

source $MUMMI_APP/setup/setup.env.sh
deactivate
spack load gromacs@2024.01~mummi-c2~mpi+rocm target=$(spack arch --target)
source $MUMMI_ROOT/$MUMMI_VENV/bin/activate

# # Distributed setup
env | grep *ROCM*
env | grep *HIP*
env | grep GPU_DEVICE_ORDINAL

env > env.txt

# export WORLD_SIZE=$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))
# echo "WORLD_SIZE = "$WORLD_SIZE

# master_addr=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n 1)
# export MASTER_ADDR=$master_addr
# export MASTER_PORT=$(expr 10000 + $(echo -n $SLURM_JOBID | tail -c 4) + $SLURM_ARRAY_TASK_ID)
# echo "MASTER_ADDR="$MASTER_ADDR
# echo "MASTER_PORT="$MASTER_PORT

# Based on https://docs.olcf.ornl.gov/software/analytics/pytorch_frontier.html#example
# Get address of head node
export MASTER_ADDR=$(hostname -i)
export NCCL_SOCKET_IFNAME=hsn0
export MASTER_PORT=3442

# Needed to bypass MIOpen, Disk I/O Errors
export MIOPEN_USER_DB_PATH="/tmp/my-miopen-cache"
export MIOPEN_CUSTOM_CACHE_DIR=${MIOPEN_USER_DB_PATH}
rm -rf ${MIOPEN_USER_DB_PATH}
mkdir -p ${MIOPEN_USER_DB_PATH}

# export "CUDA_VISIBLE_DEVICES" = "0,1,2,3"
# export "ROCR_VISIBLE_DEVICES" = "0,1,2,3"

## Launch multi-GPU generation
# python3 -m torch.distributed.launch --nnodes $SLURM_NTASKS_PER_NODE --use-env test-generators.py -n 1000

srun -N1 -n8 -c7 --gpus-per-task=1 --gpu-bind=closest python3 -W ignore -u $MUMMI_APP/tests/test-generators.py -n 10