#!/bin/bash
#SBATCH -N 2
#SBATCH -t 8:00:00
#SBATCH -p erl
#SBATCH -J v2cg-generate
#SBATCH -o %x-%j.log
#SBATCH --exclusive          # get all gpus and cpu cores per node
#SBATCH --ip-isolate yes     # for compliance with security policy

# Load environment
module load cuda/12.6
source /p/vast1/$USER/miniconda3-matrix/bin/activate
conda activate admirral_torch271

# Distributed setup
export MASTER_ADDR=$(scontrol show hostname ${SLURM_NODELIST} | head -n1)
export MASTER_PORT=23456
export OMP_NUM_THREADS=16 # change this depending on what your machine can provide
export NUM_GPUS_PER_NODE=4 # change this to match your machine

# Launch multi-node, multi-GPU generation
SRUN="srun \
    --unbuffered \
    --mpibind=off \
    --gpu-bind=none"
TORCHRUN="torchrun \
    --nnodes=$SLURM_NNODES \
    --nproc_per_node=$NUM_GPUS_PER_NODE \
    --rdzv_id=5487 \
    --rdzv_backend=c10d \
    --rdzv_endpoint=$MASTER_ADDR:$MASTER_PORT"
$SRUN $TORCHRUN \
    generate.py \
    --v-file /usr/workspace/$USER/projects/admirral/v2cg/sampled/true-v_C.20250429.npy \
    --first-n-samples 1000 \
    --output-dir ./sampled/C/ \
    --ucg-generator /usr/workspace/$USER/projects/admirral/v2ucg/lit_logs/BCR1/version_2/checkpoints/epoch=99-step=1000000.ckpt \
    --cg-generator /usr/workspace/$USER/projects/admirral/ucg2cg/lit_logs/BCR1/version_2/checkpoints/epoch=99-step=1000000.ckpt