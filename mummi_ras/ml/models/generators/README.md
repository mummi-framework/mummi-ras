# RasDoubleDiff

RasDoubleDiff is a two-stage, protein descriptor-to-structure generative model, with an emphasis on generative interpolation for protein pathway elucidation. RasDoubleDiff has been tested on RAS-RAF-14-3-3 protein complexes.
- RasDoubleDiff consists of two models: a descriptor-to-UCG generator and a UCG-to-CG generator. The terms UCG (ultra-coarse-grained) and CG (coarse-grained) refer to different levels of resolution.
- For both generators, the architecture is a modified version of a standard transformer (pre-norm, rotary positional embeddings, SiLU activation) for encoding/decoding protein structures.
- The models were trained from scratch based on the diffusion model formalism, namely the denoising score matching objective.
- The application of generative interpolation is achieved by performing descriptor-to-structure generation from a descriptor value between two end points of a pathway.


## Installation

Dependencies:
- PyTorch
- PyTorch-Lightning: for facilitating model training and logging without the need to write manual training loops. It's also used for defining custom model inference behaviors.
- PyTorch-Geometric: used for optimized scatter-reduce operation for CG-to-UCG aggregation.


### Example installation steps

Refer to the PyTorch website for installing `torch`. Then, clone this repo and cd to the repo folder. Lastly, run:
```
pip install ".[rasdiff]"
```


## How to use

The `generate.py` python script loads the generative model checkpoints and generates RAS-RAF-14-3-3 structures based on a given numpy .npy file containing an array of descriptors (also called `v`), in the shape of `[N, 4]`. Effectively, the generation process is v → UCG → CG. This means `N` structures will be generated corresponding to the number of descriptors in `v`. Also, the script assumes a distributed setup so it should be launched with `torchrun`. In the case of multi-node jobs, an additional command for launching parallel jobs is typically required, e.g., `srun` for the SLURM scheduler. Overall, a shell script for launching a multi-node, multi-GPU generation process may look like this:

```bash
# Load python environment (conda)
source /path/to/conda-dir/bin/activate
conda activate your_torch_env

# Distributed setup (assuming SLURM is used)
export MASTER_ADDR=$(scontrol show hostname ${SLURM_NODELIST} | head -n1)
export MASTER_PORT=23456
export OMP_NUM_THREADS=16 # change this depending on what your machine can provide
export NUM_GPUS_PER_NODE=4 # change this to match your machine

# SLURM srun command
SRUN="srun \
    --unbuffered \
    --mpibind=off \
    --gpu-bind=none"

# torchrun command
TORCHRUN="torchrun \
    --nnodes=$SLURM_NNODES \
    --nproc_per_node=$NUM_GPUS_PER_NODE \
    --rdzv_id=5487 \
    --rdzv_backend=c10d \
    --rdzv_endpoint=$MASTER_ADDR:$MASTER_PORT"

# Launch multi-node, multi-GPU generation
$SRUN $TORCHRUN \
    generate.py \
    --v-file /path/to/npy_file/containing/v_array \
    --first-n-samples 1000 \ # limit to the first n samples of the v array for debugging purposes
    --output-dir /path/to/output_dir \
    --ucg-generator /path/to/v-to-ucg-generator/checkpoint-file.ckpt \
    --cg-generator /path/to/ucg-to-cg-generator/checkpoint-file.ckpt
```

The shell script `generate.sh` is provided as an example. Note that it assumes the use of Nvidia cuda and the SLURM scheduler.


### Model checkpoint files and v examples

An example .npy file for `v`, a model checkpoint for the v-to-UCG generator, and a model checkpoint for the UCG-to-CG generator are provided at ... 


### Converting raw generation output to protein data format

The as-generated outputs are saved in .npy files, which contain the UCG/CG bead coordinates in the shape of `[N, num_beads, 3]`. To convert to file formats for ease of parsing and visualization, refer to the following steps:
...
