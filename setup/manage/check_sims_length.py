#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Copyright (c) 2023, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------

import argparse
import logging
import os
from pathlib import Path
from struct import unpack
import traceback
import sys
import time
import glob
from multiprocessing.pool import ThreadPool as Pool
from typing import List

import numpy as np

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)

def plot(data: np.array, output: str, max_ns: int = 600, timestr: str = None) -> None:
    try:
        import matplotlib.mlab as mlab
        import matplotlib.pyplot as plt
    except Exception as e:
        logging.error(f"Did not find Matplolib")
        return

    font_size=14
    valid_x = data[data >= max_ns]

    plt.figure(figsize=(12,6))
    plt.style.use('default')
    # plt.style.use('seaborn-whitegrid')
    plt.title(f'Simulated length ({len(data)} #simulations | {len(valid_x)} >= {max_ns} ns | {np.sum(data) / (1000*1000)} ms)', fontsize=font_size)
    plt.xlabel('Simulation length (ns)', fontsize=font_size)
    plt.xticks(fontsize=font_size)
    plt.ylabel('#Simulations', fontsize=font_size)
    plt.yticks(fontsize=font_size)

    n, bins, patches = plt.hist(data, bins=25, facecolor='#2ab0ff',edgecolor='#169acf', linewidth=1)
    plt.axvline(x = np.nanmean(data), color='r', linestyle='dashed', linewidth=1.3)

    if not timestr:
        timestr = time.strftime("%Y%m%d-%H%M%S")

    plt.savefig(f"{output}-{timestr}.pdf")
    plt.savefig(f"{output}-{timestr}.png", transparent=True)
    print(f"Outputs: {output}-{timestr}.png {output}-{timestr}.pdf")

def get_xtc_length_ns(xtc_file: str) -> float:
    try:
        with open(xtc_file, 'rb') as f:
            tag = f.read(8)                           # Tag: magic number and number of atoms
            n = 92 + unpack(">l",f.read(84)[-4:])[0]  # Size of frame in bytes

            f.seek(int(-5*n/4), 2)                    # This should contain a complete frame
            frame = f.read()                          # Read the remaining part in
            n = frame.index(tag)                      # Find the tag

            return 0.001*unpack(">f",frame[n+12:n+16])[0]
    except Exception as excpt:
        return -1

def get_step_id(xtc_file: str) -> int:
    return int(xtc_file.split(".xtc")[0].split("step")[-1])

def worker(entry: str) -> float:
    xtc_main = glob.glob(os.path.join(entry, "traj_comp.xtc"))
    xtc_steps = glob.glob(os.path.join(entry, "traj_comp.step*.xtc"))
    if xtc_main == []:
        return -2.0

    if len(xtc_main) != 1:
        logger.error(f"Found {len(xtc_main)} traj_comp.xtc in {entry}")
        return -2.0

    latest_xtc = ''
    if len(xtc_steps) > 0:
        latest_xtc = xtc_steps[0]
        for x in xtc_steps:
            if  get_step_id(x) > get_step_id(latest_xtc):
                latest_xtc = x

    main_length = get_xtc_length_ns(xtc_main[0])
    latest_xtc = get_xtc_length_ns(latest_xtc)

    # stopped = len(glob.glob(os.path.join(entry, "cg_stop"))) == 1

    # grep "Performance:      XXX" in mdrun.log

    return max(main_length, latest_xtc)

def print_failed(caught):
    traceback.print_exc(file=sys.stderr)

def get_all_length_sims(paths: List[str], pool_size: int) -> dict:
    start = time.time()
    list_dirs = []
    for x in paths:
        if not os.path.isdir(x):
            print(f"Skipping {x} (not a valid directory)")
            continue
        x = Path(x)
        curr = list(os.scandir(x))
        print(f"Found {len(curr)} in {x}")
        list_dirs += curr
    end = time.time() - start

    futures = {}
    res = {}

    print(f"Found {len(list_dirs)} in {end:.3f} seconds")

    start = time.time()
    with Pool(pool_size) as pool:
        for dir in list_dirs:
            if dir.is_dir():
                futures[dir.path] = pool.apply_async(worker, (dir.path,), error_callback = print_failed)
        
        for k,v in futures.items():
            res[k] = v.get()
    end = time.time() - start

    print(f"Processed {len(res)} in {end:.3f} seconds")
    return res

def main(args):
    results = get_all_length_sims(paths = args.dirs, pool_size = args.threads)
    data = []

    for k,v in results.items():
        if args.verbose:
            print("{:20s} {:10.2f} ns".format(k, v))
        data.append(v)

    data = np.array(data)
    output = "hist-mummi-c3"
    timestr = time.strftime("%Y%m%d-%H%M%S")
    np.savez(f"{output}-{timestr}.npz", data=data)

    print("")
    sims_notstarted = [k for k,v in results.items() if v < 0]
    sims_started = [k for k,v in results.items() if v > 0]
    print(f"# Simulations                    : {len(results):10d}")
    print(f"# CG simulations not started     : {len(sims_notstarted):10d}")
    print(f"# CG simulations running         : {len(sims_started):10d}")
    if args.max_sim_time > 0:
        sims_finished = [k for k,v in results.items() if v >= args.max_sim_time]
        print(f"# CG simulations > {args.max_sim_time:<5d} ns      : {len(sims_finished):10d}")

    if args.plot:
        plot(data[data > 0], output = output, max_ns = args.max_sim_time, timestr = timestr)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script that iterates in a directory containing GROMACS run directories")
    parser.add_argument("-d", "--dirs", help="Directories where to parse XTC", required=True, nargs='+')
    parser.add_argument("-t", "--threads", help="Number of threads to use", required=False, type=int, default=8)
    parser.add_argument("-m", "--max-sim-time", help="max simulation time to flag", required=False, type=int, default=600)
    parser.add_argument("-p", "--plot", help="Plot histogram", action='store_true')
    parser.add_argument("-v", "--verbose", help="Verbose mode", action='store_true')

    args = parser.parse_args()

    main(args)


