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
import os
from pathlib import Path
from struct import unpack
import traceback
import sys
import time
import glob
import warnings
import threading
from multiprocessing.pool import ThreadPool as Pool
from typing import List, Tuple

import numpy as np
import MDAnalysis as mda

# logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)

FILES_DELETE_SUCCESS = [
    "createsims.out",
	"cg_analysis.out",
	"analysis.tar.pylst",
	"analysis.tar.pytree",
	"martini_v2.x_new-rf-prod.mdp",
	"ppatch.npz",
	"state_prev.cpt",
	"'#*",
	"#traj_comp*",
	"restart",
	".traj_comp*",
	"traj_comp.step*.xtc",
    ".topol.tpr.*"
]

FILES_KEEP_FAILURE = [
	"struct_createsim.gro",
	"lipid_counts.npz",
	"lipids-water.gro",
	"createsims.out",
	"createsims.log",
	"system.top",
	"psim.npz",
	"macro_sim_dir/macro_createsims.out",
	"macro_sim_dir/snapshot.*"
]

def get_xtc_length_ns(xtc_file: str) -> float:
    try:
        with open(xtc_file, 'rb') as f:
            tag = f.read(8)                           # Tag: magic number and number of atoms
            n = 92 + unpack(">l",f.read(84)[-4:])[0]  # Size of frame in bytes

            f.seek(int(-5*n/4), 2)                    # This should contain a complete frame
            frame = f.read()                          # Read the remaining part in
            n = frame.index(tag)                      # Find the tag

            return 0.001*unpack(">f",frame[n+12:n+16])[0]
    except Exception as e:
        # try: mda
        res = mda_get_xtc_length_ns(xtc_file)
        if res == -1:
            print(f"Error: while reading XTC file {xtc_file}: {e}", file=sys.stderr)
        return res

# Somehow segfault after some iteration
def mda_get_xtc_length_ns(xtc_file: str) -> float:
    try:
        with warnings.catch_warnings(action="ignore"):
            u = mda.Universe(xtc_file)
            # we convert from ps to ns
            length = (u.trajectory.dt * u.trajectory.n_frames) / 1000
            return float(length)
    except Exception as e:
        print(f"Error: while reading {xtc_file} : {e}", file=sys.stderr)
        return -1

def get_step_id(xtc_file: str) -> int:
    return int(xtc_file.split(".xtc")[0].split("step")[-1])

def check_xtc_validity(entry: str, tolerance_xtc: float = 5) -> int:
    xtc_main = glob.glob(os.path.join(entry, "traj_comp.xtc"))
    xtc_steps = glob.glob(os.path.join(entry, "traj_comp.step*.xtc"))

    if xtc_main == []:
        # print(f"Warning: missing traj_comp.xtc in {entry}", file=sys.stderr)
        return 4

    if len(xtc_main) != 1:
        print(f"Warning: found {len(xtc_main)} traj_comp.xtc in {entry}", file=sys.stderr)
        return 3

    xtc_main = xtc_main[0]
    latest_xtc = ''
    if len(xtc_steps) > 0:
        latest_xtc = xtc_steps[0]
        for x in xtc_steps:
            if  get_step_id(x) > get_step_id(latest_xtc):
                latest_xtc = x
    else:
        # print(f"Warning: missing traj_comp.step*.xtc files", file=sys.stderr)
        return 2

    latest_xtc_ns = get_xtc_length_ns(latest_xtc)
    main_length = get_xtc_length_ns(xtc_main)

    if (latest_xtc_ns - tolerance_xtc <= main_length) and (main_length <= latest_xtc_ns + tolerance_xtc):
        return 0
    else:
        print(
            f"Warning: [{xtc_main}] => {main_length} ns  != {latest_xtc} => {latest_xtc_ns} ns (th = {tolerance_xtc})",
            file=sys.stderr
        )
        return 1

def successful_dir(path: str):
    log = os.path.isfile(os.path.join(path, "cg_analysis.log"))
    out = os.path.isfile(os.path.join(path, "cg_analysis.out"))
    return log or out

def expand_list(base_path: str, list_files: List[str]) -> List[str]:
    """
    Return a similar list to list_files but if list_files contains pattern like *
    """
    res = []
    for f in list_files:
        pattern_expanded = list(glob.glob(os.path.join(base_path, f)))
        res += pattern_expanded
    return res

def remove_files(dir: os.DirEntry, files_to_remove: List[str], bck_dir: str = None) -> dict:
    """
    Return a dict with the amount of bytes deleted and the number of files deleted
    {'total_file', 'files_deleted', 'size_deleted'}
    """
    abs_path = os.path.abspath(dir.path)
    patterns = expand_list(abs_path, files_to_remove)
    res = {"path": abs_path, "total_file": 0, "total_size": 0, "files_deleted": 0, "size_deleted": 0}

    if bck_dir is not None:
        new_bck_dir = os.path.join(bck_dir, os.path.relpath(dir.path))
        os.makedirs(new_bck_dir, exist_ok = True)

    for item in os.scandir(dir):
        if item.is_file():
            fname = os.path.join(abs_path, item.name)
            stat = item.stat()
            res["total_file"] += 1
            res["total_size"] += stat.st_size
            if fname in patterns:
                res["files_deleted"] += 1
                res["size_deleted"] += stat.st_size
                if bck_dir is None:
                    print(f"{fname},{stat.st_size}", file=sys.stderr)
                else:
                    try:
                        new_fname = os.path.join(new_bck_dir, item.name)
                        # Atomic if src and dst on same filesystem
                        os.rename(fname, new_fname)
                        print(f"{fname},{new_fname},{stat.st_size}", file=sys.stderr)
                    except Exception as e:
                        print(f"Error: rename {fname}: {e}", file=sys.stderr)
                    
    return res

def keep_files(dir: os.DirEntry, files_to_keep: List[str], bck_dir: str = None) -> dict:
    """
    Return a dict with the amount of bytes deleted and the number of files deleted
    {'total_file', 'files_deleted', 'size_deleted'}
    """
    res = {"path": "", "total_file": 0, "total_size": 0, "files_deleted": 0, "size_deleted": 0}
    if dir is None:
        return res
    abs_path = os.path.abspath(dir.path)
    patterns = expand_list(abs_path, files_to_keep)
    res["path"] = abs_path

    if bck_dir is not None:
        new_bck_dir = os.path.join(bck_dir, os.path.relpath(dir.path))
        os.makedirs(new_bck_dir, exist_ok = True)

    for item in os.scandir(dir):
        if item.is_file():
            fname = os.path.join(abs_path, item.name)
            stat = item.stat()
            res["total_file"] += 1
            res["total_size"] += stat.st_size

            if not fname in patterns:
                res["files_deleted"] += 1
                res["size_deleted"] += stat.st_size
                if bck_dir is None:
                    print(f"{fname},{stat.st_size}", file=sys.stderr)
                else:
                    try:
                        new_fname = os.path.join(new_bck_dir, item.name)
                        # Atomic if src and dst on same filesystem
                        os.rename(fname, new_fname)
                        print(f"{fname},{new_fname},{stat.st_size}", file=sys.stderr)
                    except Exception as e:
                        print(f"Error: rename {fname} => {new_fname}: {e}", file=sys.stderr)
        else:
            # we jump in the directory. Not recursive because max depth = 1
            abs_path_maccro = os.path.join(abs_path, item.name)
            patterns_maccro = expand_list(abs_path, files_to_keep)

            for macro_file in os.scandir(item):
                fname = os.path.join(abs_path_maccro, macro_file.name)
                stat = macro_file.stat()
                res["total_file"] += 1
                res["total_size"] += stat.st_size

                if fname in patterns_maccro:
                    res["files_deleted"] += 1
                    res["size_deleted"] += stat.st_size
                    if bck_dir is None:
                        print(f"{fname},{stat.st_size}", file=sys.stderr)
                    else:
                        try:
                            new_fname = os.path.join(new_bck_dir, macro_file.name)
                            # Atomic if src and dst on same filesystem
                            os.rename(fname, new_fname)
                            print(f"{fname},{new_fname},{stat.st_size}", file=sys.stderr)
                        except Exception as e:
                            print(f"Error: rename {fname}: {e}", file=sys.stderr)
    return res

def clean_entry(entry: os.DirEntry, bck_dir: str = None) -> dict:
    res = {}
    if successful_dir(entry.path):
        # Success
        res = remove_files(dir = entry, files_to_remove = FILES_DELETE_SUCCESS, bck_dir = bck_dir)
    else:
        # Failure
        failed_entry = None
        for x in os.scandir(entry):
            if x.is_dir():
                failed_entry = x

        if failed_entry is None:
            print(f"Warning: {entry.name} seems empty", file=sys.stderr)
        else:
            res = keep_files(dir = failed_entry, files_to_keep = FILES_KEEP_FAILURE, bck_dir = bck_dir)
    return res  

def worker(entry: os.DirEntry, tolerance_xtc: float, bck_dir: str = None) -> dict:
    success = successful_dir(entry.path)
    if success:
        valid_entry = check_xtc_validity(entry.path, tolerance_xtc)
    res = clean_entry(entry, bck_dir)
    if success:
        res["valid_xtc"] = valid_entry
    else:
        res["valid_xtc"] = 5 # 5 is code for failed simulation
    return res

def dir_not_empty(entry: os.DirEntry) -> bool:
    return list(os.scandir(entry)) != []

def dir_is_valid(entry: os.DirEntry) -> bool:
    """
    A valid directory is one with one of these files:
    - cg_analysis.out
    - createsims_failure
    """
    cg_analysis = os.path.isfile(os.path.join(entry.path, "cg_analysis.out")) or os.path.isfile(os.path.join(entry.path, "cg_analysis.log"))
    createsim_failure = os.path.isfile(os.path.join(entry.path, "createsims_failure"))
    return cg_analysis or createsim_failure

def get_all_length_sims(
        path: str,
        pool_size: int,
        tolerance: float,
        dry_run: bool = True,
        max_sims: int = None,
        bck_dir: str = None
        ) -> Tuple[dict, float]:
    if not os.path.isdir(path):
        print(f"Error: {path} is not a valid directory", file=sys.stderr)
        return ({}, 0)

    path = Path(path)
    start = time.time()
    list_dirs = list(os.scandir(path))

    if max_sims:
        print(f"Selected the top {max_sims} sims")
        list_dirs = list_dirs[:max_sims]

    end = time.time() - start
    futures = {}
    res = {}

    print(f"Found {len(list_dirs)} in {end:.3f} seconds")

    def print_failed(caught):
        traceback.print_exc(file=sys.stdout)

    if dry_run:
        print(f"!!! Dry run (files won't be deleted) !!!")

    start = time.time()
    if pool_size == 1:
        for dir in list_dirs:
            if dir_is_valid(dir):
                res[dir.path] = worker(dir, tolerance, bck_dir)
    else:
        with Pool(pool_size) as pool:
            for dir in list_dirs:
                if dir_is_valid(dir):
                    futures[dir.path] = pool.apply_async(worker, (dir, tolerance, bck_dir), error_callback = print_failed)

            for k,v in futures.items():
                res[k] = v.get()

    end = time.time() - start
    return res, end

def main(args):
    bck_dir = None
    if not args.dry_run:
        bck_dir = os.path.abspath(args.backup_dir)
        os.makedirs(bck_dir, exist_ok = True)
        print(f"Files will be moved to {bck_dir}")

    sims_dir = os.path.abspath(args.dir)

    if not os.path.isdir(sims_dir):
        print(f"Error: {sims_dir} is not a valid directory", file=sys.stderr)
        return

    results, time_taken = get_all_length_sims(
        path = sims_dir,
        pool_size = args.threads,
        tolerance = args.tolerance,
        dry_run = args.dry_run,
        max_sims = args.max_sims,
        bck_dir = bck_dir
    )

    print("")
    sims_compliant = [k for k,v in results.items() if v["valid_xtc"] == 0]
    sims_unvalid_xtc = [k for k,v in results.items() if v["valid_xtc"] == 1]
    sims_others = [k for k,v in results.items() if v["valid_xtc"] in [2,3,4]]
    sims_failed = [k for k,v in results.items() if v["valid_xtc"] == 5]
    th = len(results) / time_taken

    total_files = int(sum([v["total_file"] for k,v in results.items()]))
    total_size = sum([v["total_size"] for k,v in results.items()]) / (1024**3)
    deleted_files = int(sum([v["files_deleted"] for k,v in results.items()]))
    deleted_size = sum([v["size_deleted"] for k,v in results.items()]) / (1024**3)

    total_sims_xtc = len(sims_compliant)+len(sims_unvalid_xtc)

    print(f"Total simulations         : {len(results):>10d}")
    print(f"Failed simulations        : {len(sims_failed):>10d}")
    print(f"Processed in (s)          : {time_taken:>10.3f}")
    print(f"Throughput (sim/s)        : {th:>10.3f}")
    print(f"---------------------------")
    print(f"Total sims with XTC       : {total_sims_xtc:>10d}")
    print(f"Valid XTCs                : {len(sims_compliant):>10d}")
    print(f"Tolerance (ns)        +/- : {args.tolerance:>10.3f}")
    print(f"Unvalid XTCs              : {len(sims_unvalid_xtc):>10d}")
    print(f"Other sims status         : {len(sims_others):>10d}")
    print(f"---------------------------")
    print(f"# Total files             : {total_files:>10d}")
    print(f"# Deleted files           : {deleted_files:>10d}")
    if total_files > 0:
        print(f"  Ratio (%)               : {100*deleted_files/total_files:>10f}")
    print(f"---------------------------")
    print(f"# Total size (GB)         : {total_size:>10.3f}")
    print(f"# Deleted size (GB)       : {deleted_size:>10.3f}")
    if total_size > 0:
        print(f"  Ratio (%)               : {100*deleted_size/total_size:>10.3f}")


    timestr = time.strftime("%Y%m%d-%H%M%S")
    np.savez(f"data-cleanup-{timestr}.npz", **results)  # data is a dict here

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script that cleanup files from simulation directory")
    parser.add_argument("-d", "--dir", help="Directory where to clean up files", required=True, type=str, default="sims-cg")
    parser.add_argument("-t", "--threads", help="Number of threads to use", required=False, type=int, default=1)
    parser.add_argument("-e", "--tolerance", help="Tolerance threshold to check XTC validity", required=False, type=float, default=5)

    parser.add_argument(
        "-m",
        "--max-sims",
        help="Only runs for this number of simulations (useful for testing)",
        required=False,
        type=int,
        default=None
    )

    parser.add_argument("--dry-run", help="Perform a dry-run (do not delete files)", action="store_true")
    parser.add_argument("--backup-dir", '-b', help="Where to move all the files", required=False, type=str)

    args = parser.parse_args()

    if not args.dry_run:
        if args.backup_dir is None:
            parser.error('Without --dry-run, --backup-dir is expected')      

    main(args)
