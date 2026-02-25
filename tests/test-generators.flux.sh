#!/bin/bash
#flux: -N 4
#flux: -t 8h
#flux: -q pbatch
#flux: --job-name=test-generators
#flux: --output={{name}}-{{id}}.log

# Environment setup
module load rocm/6.3
# source /usr/workspace/hsu16/miniconda/bin/activate
# conda activate torch270_rocm

# Distributed setup
export MASTER_ADDR=$(flux getattr hostlist | /bin/hostlist -n 1)
export MASTER_PORT=23456
export NUM_NODES=$(flux getattr size)

flux run -N $NUM_NODES -n $NUM_NODES --unbuffered --exclusive \
    -ompibind=omp_proc_bind,omp_places \
    torchrun \
    --nnodes=$NUM_NODES \
    --nproc_per_node=4 \
    --rdzv_id=5487 \
    --rdzv_backend=c10d \
    --rdzv_endpoint=$MASTER_ADDR:$MASTER_PORT \
    test-generators.py \
    --num-samples 1000 \
    --batch-size 32 \
    --steps 256