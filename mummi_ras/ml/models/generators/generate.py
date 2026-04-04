# LC hack: workaround so that "import torch" will not change CPU affinity
import os
if 'OMP_PLACES'    in os.environ: del os.environ['OMP_PLACES']
if 'OMP_PROC_BIND' in os.environ: del os.environ['OMP_PROC_BIND']

import torch
import numpy as np
import argparse
from tqdm import tqdm
from torch.utils.data import DataLoader
from datetime import datetime, timedelta
from c3.ucg_generator import V2UCG_Generator
from c3.cg_generator import UCG2CG_Generator


def main():
    # CLI arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--v-file', type=str, help='path to numpy .npy file containing v arrays in the shape of [N, 4]')
    parser.add_argument('--first-n-samples', type=int, help='only use the first n samples from v-file', default=None)
    parser.add_argument('--output-dir', type=str, help='path to output directory')
    parser.add_argument('--ucg-generator', type=str, help='path to v-to-UCG generator checkpoint file')
    parser.add_argument('--cg-generator', type=str, help='path to UCG-to-CG generator checkpoint file')
    parser.add_argument('--ucg-index-file', type=str, help='plain text file describing UCG-CG mapping', default='./configs/ucg-idx-dimer.txt')
    args = parser.parse_args()

    # DDP setup
    torch.distributed.init_process_group(
        backend='nccl',
        timeout=timedelta(hours=4), # use long timeout duration because some GPUs will finish inference earlier than others
        init_method=f'env://', # use environment variables 
        rank=int(os.environ.get('RANK', 0)),
        world_size=int(os.environ.get('WORLD_SIZE', 1)),
    )
    rank = torch.distributed.get_rank()
    world_size = torch.distributed.get_world_size()
    device_id = rank % torch.cuda.device_count()
    device = torch.device(f'cuda:{device_id}')
    print(f'{world_size=}, {rank=}, {device=}', flush=True)

    # Data setup
    v = np.load(args.v_file)
    if args.first_n_samples is not None:
        v = v[:args.first_n_samples]
    if rank == 0:
        print(f'Generating from {v.shape[0]} descriptors (v)', flush=True)
    # logv = torch.tensor(v, dtype=torch.float).log10()  # the v-to-UCG generator takes log(v) values instead of v
    v = np.array_split(v, world_size)[rank]
    loader = DataLoader(v, batch_size=32)

    # Model setup
    v2ucg_generator = V2UCG_Generator.load_from_checkpoint(
        args.ucg_generator,
        map_location='cpu',
    ).to(device)
    ucg2cg_generator = UCG2CG_Generator.load_from_checkpoint(
        args.cg_generator,
        ucg_index_file=args.ucg_index_file,
        map_location='cpu',
    ).to(device)
    v2ucg_generator.eval()
    ucg2cg_generator.eval()

    # Start generations
    pred_ucg, pred_cg = [], []
    with torch.inference_mode():
        for batch in tqdm(loader, position=rank, desc=f'GPU {rank}'):
            pred_ucg_i = v2ucg_generator.generate(batch.to(device), steps=256)
            # The v-to-UCG generator and the UCG-to-CG generator operate at different scales
            # For the v-to-UCG generator, a scale of 1/20 is applied to the protein UCG bead coordinates
            pred_ucg_i *= 20 
            pred_cg_i = ucg2cg_generator.generate(pred_ucg_i, steps=256)
            pred_ucg.append(pred_ucg_i.clone())
            pred_cg.append(pred_cg_i.clone())

    # Concatenate the generation outputs per GPU process
    pred_ucg = torch.cat(pred_ucg)
    pred_cg = torch.cat(pred_cg)

    # Gather up all outputs from all GPUs
    gather_list = [torch.zeros_like(pred_ucg) for _ in range(world_size)]
    torch.distributed.all_gather(gather_list, pred_ucg)
    all_pred_ucg = torch.cat(gather_list).cpu().numpy()

    gather_list = [torch.zeros_like(pred_cg) for _ in range(world_size)]
    torch.distributed.all_gather(gather_list, pred_cg)
    all_pred_cg = torch.cat(gather_list).cpu().numpy()

    # Save generation outputs
    if rank == 0:
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)
        timestamp = datetime.now().strftime('%Y%m%d')
        np.save(os.path.join(args.output_dir, f'pred-ucg.{timestamp}.npy'), all_pred_ucg)
        np.save(os.path.join(args.output_dir, f'pred-cg.{timestamp}.npy'), all_pred_cg)

    torch.distributed.destroy_process_group()


if __name__ == '__main__':
    main()