#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2023, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------
# @authors
#           Loic Pottier <pottier1@llnl.gov>
# ------------------------------------------------------------------------------

# LC hack: workaround so that "import torch" will not change CPU affinity
import os
if 'OMP_PLACES'    in os.environ: del os.environ['OMP_PLACES']
if 'OMP_PROC_BIND' in os.environ: del os.environ['OMP_PROC_BIND']

# LC hack for distributed setup
if 'OMPI_COMM_WORLD_RANK'       in os.environ: os.environ["RANK"]       = os.environ['OMPI_COMM_WORLD_RANK']
if 'OMPI_COMM_WORLD_SIZE'       in os.environ: os.environ["WORLD_SIZE"] = os.environ['OMPI_COMM_WORLD_SIZE']
if 'OMPI_COMM_WORLD_LOCAL_RANK' in os.environ: os.environ["LOCAL_RANK"] = os.environ['OMPI_COMM_WORLD_LOCAL_RANK']

# Slurm
if 'SLURM_PROCID' in os.environ: os.environ["LOCAL_RANK"] = os.environ["SLURM_PROCID"]

import torch
import numpy as np
import argparse
from tqdm import tqdm
from logging import getLogger
from torch.utils.data import DataLoader
import time
import sys

import mummi_core
import mummi_ras
from mummi_ras import Naming
from mummi_ras.ml.models import V_Generator, UCG2CG_Generator, V2UCG_Generator

from mummi_ras.ml.models.autoencoders import FullDenseAutoencoderUCG

LOGGER = getLogger(__name__)

def main():
    # CLI arguments
    parser = argparse.ArgumentParser(description="Script to test sampler, generator and validator routines.")
    parser.add_argument("-n", "--num-samples", type=int, required=True,               help="Number of samples requested")
    parser.add_argument("-b", "--batch-size",  type=int, required=False, default=32,  help="Batch size (per-GPU)")
    parser.add_argument("-s", "--steps",       type=int, required=False, default=128, help="Number of generation/denoising steps")
    parser.add_argument("-o", "--output",      type=str, required=False, default="table_all_generations.npz", help="Output file (npz)")
    args = parser.parse_args()

    # DDP setup
    torch.cuda.set_device(int(os.environ["LOCAL_RANK"]))
    torch.distributed.init_process_group(backend='nccl')
    rank = torch.distributed.get_rank()
    world_size = torch.distributed.get_world_size()
    device_id = rank % torch.cuda.device_count()
    device = torch.device(f'cuda:{device_id}')
    print(f'{world_size=:}, {rank=:}, {device=:}', flush=True)

    # MUMMI initialization and reading mlserver config
    mummi_core.init()
    mummi_core.create_root()
    spath = os.path.join(Naming.MUMMI_SPECS, "workflow")
    config = mummi_core.read_specs(os.path.join(spath, "mlserver.yaml"))
    if rank == 0:
        mummi_core.init_logger(config=config["config"])
        validation_time = time.time()

    # Determine model checkpoint file paths depending on the machine (e.g., lassen vs. frontier) 
    v2ucg_ckpt = None
    for path in config["models"]["generative"]["v2ucg"]["paths"]:
        if os.path.isfile(path):
            v2ucg_ckpt = path
    assert v2ucg_ckpt is not None

    ucg2cg_ckpt = None
    for path in config["models"]["generative"]["ucg2cg"]["paths"]:
        if os.path.isfile(path):
            ucg2cg_ckpt = path
    assert ucg2cg_ckpt is not None


    ucg_model_dir = None
    for path in config["models"]["ucg_autoencoder"]["paths"]:
        print(f"path= {path} => {os.path.isdir(path)}")
        ucg_model_dir = path if os.path.isdir(path) else None

    if ucg_model_dir is None:
        print(f"Error: did not find any valid paths for the UCG model")
        sys.exit(-1)

    ucg_model_path = os.path.join(ucg_model_dir, config["models"]["ucg_autoencoder"]["model"])

    indice_cluster = os.path.join(config["models"]["ucg_autoencoder"]["resources"], "indices_per_cluster_dimer.npz")

    ucg_encoder = FullDenseAutoencoderUCG(
        model_path = ucg_model_path,
        indice_per_cluster = indice_cluster
    )

    ##### Generation ######

    # Generate v samples
    v_generator = V_Generator(seed=0)
    v = v_generator.generate(10*2, args.num_samples//10)  # this should result in more than num_samples of v's

    # More than num_samples of v's were generated intentionally. Trim down to num_samples
    random_indices = torch.randperm(v.size(0))
    v = v[random_indices]
    v = v[:args.num_samples]

    # Split v's into num_gpu chunks, each on a different GPU
    assert len(v) == args.num_samples, 'After filtering, size of v is no longer the specified number'
    loader = DataLoader(torch.tensor_split(v, world_size)[rank], batch_size=args.batch_size)

    # Load generators
    model_config_path = config["models"]["generative"]["resources"]
    v2ucg_generator = V2UCG_Generator.load_from_checkpoint(v2ucg_ckpt).to(device)
    ucg2cg_generator = UCG2CG_Generator.load_from_checkpoint(
        ucg2cg_ckpt, ucg_index_file=os.path.join(model_config_path, "dimers-ucg-idx.txt")
    ).to(device)
    v2ucg_generator.eval()
    ucg2cg_generator.eval()

    # Generate UCG and CG structures
    if rank == 0:
        LOGGER.info(f"Generating {len(v)} structures")
    pred_ucg, pred_cg = [], []
    with torch.inference_mode():
        for batch in tqdm(loader, position=rank, desc=f'GPU {rank}'):
            pred_ucg_i = v2ucg_generator.generate(batch.to('cuda'), steps=args.steps)
            pred_ucg_i *= 20
            pred_cg_i = ucg2cg_generator.generate(pred_ucg_i, steps=args.steps)
            pred_ucg.append(pred_ucg_i.clone())
            pred_cg.append(pred_cg_i.clone())
    pred_ucg = torch.cat(pred_ucg, dim=0)
    pred_cg = torch.cat(pred_cg, dim=0)

    # Gather up all outputs from all GPUs
    gather_list = [torch.zeros_like(pred_ucg) for _ in range(world_size)]
    torch.distributed.all_gather(gather_list, pred_ucg)
    pred_ucg = torch.cat(gather_list).cpu().numpy()
    if rank == 0:
        LOGGER.info(f'Generated UCG structures: {pred_ucg.shape}')

    gather_list = [torch.zeros_like(pred_cg) for _ in range(world_size)]
    torch.distributed.all_gather(gather_list, pred_cg)
    pred_cg = torch.cat(gather_list).cpu().numpy()
    if rank == 0:
        validation_time = time.time() - validation_time
        LOGGER.info(f'Generated CG structures: {pred_cg.shape} in {validation_time:.2f} secs ({pred_cg.shape[0]/validation_time} sample/sec)')

    # Encode all structures into UCG latent space
    ucg_encoding_time = time.time()
    ls_coords = []
    # Could be parralel?
    for cg_position in pred_cg:
        ls_coords.append(ucg_encoder.encode(cg_position))

    ls_coords = np.array(ls_coords)
    ucg_encoding_time = time.time() - ucg_encoding_time
    LOGGER.info(f'Encoded UCG structures into LS: {pred_cg.shape} in {ucg_encoding_time:.2f} secs ({pred_cg.shape[0]/ucg_encoding_time} sample/sec)')

    # Save generation outputs
    if rank == 0:
        out_fname = args.output
        if os.path.isfile(out_fname):
            data_dict = np.load(out_fname, allow_pickle=True)
            v        = np.append(data_dict['v'],        v,        axis=0)
            pred_ucg = np.append(data_dict['pred_ucg'], pred_ucg, axis=0)
            pred_cg  = np.append(data_dict['pred_cg'],  pred_cg,  axis=0)
            ls_coords  = np.append(data_dict['ls_coords'],  ls_coords,  axis=0)
        np.savez_compressed(
            out_fname,
            model_v2ucg=[v2ucg_ckpt],
            model_ucg2cg=[ucg2cg_ckpt],
            v=v,
            pred_ucg=pred_ucg,
            pred_cg=pred_cg,
            ls_coords=ls_coords,
            selected=np.full(len(pred_cg), False)
        )
        LOGGER.info(f'Generations saved to {out_fname}')

if __name__ == "__main__":
    main()