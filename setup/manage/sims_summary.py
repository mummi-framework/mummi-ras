#!/usr/bin/env python3

from struct import unpack
import traceback
import sys
import numpy as np
import datetime
from typing import List
import re
import os

def save_db(filename: str, data: dict):
    data["created_at"] = datetime.datetime.now().strftime('%d.%m.%Y-%H:%M:%S')

    np.savez_compressed(filename, **data)

def collect_data(list_dirs: List[str]):
    sims_names = []
    frame_ns = []
    ns_per_day = []
    regexp = re.compile("Performance:")
    for dir in list_dirs:
        try:
            xtc = os.path.join(dir, "traj_comp.xtc")
            with open(xtc, 'rb') as f:
                tag = f.read(8)                           # Tag: magic number and number of atoms
                n = 92 + unpack(">l",f.read(84)[-4:])[0]  # Size of frame in bytes
                f.seek(int(-5*n/4), 2)                    # This should contain a complete frame
                frame = f.read()                          # Read the remaining part in
            n = frame.index(tag)                      # Find the tag
            frame_ns.append(float(0.001*unpack(">f",frame[n+12:n+16])[0]))
            perf = []
            with open(os.path.join(dir, "cg_analysis.out"), 'r') as f:
                for line in f:
                    if re.search(regexp, line):
                        new_perf = float([x for x in line.split(' ') if x != ''][1])
                        if perf != []:
                            print(f"Found two performance ns/day for {dir}: {perf} {new_perf}")
                        perf.append(new_perf)
            ns_per_day.append(np.mean(perf))
        except Exception as _:
            frame_ns.append(np.nan)
            ns_per_day.append(np.nan)
        finally:
            sims_names.append(xtc)
        print(sims_names[-1], frame_ns[-1], ns_per_day[-1])
    return {"structures" : sims_names, "last_frame_ns": frame_ns, "ns_per_day": ns_per_day}

if __name__ == "__main__":
    ts = datetime.datetime.now().strftime('%d%m%Y_%H%M%S')
    output = f"sims-summary-{ts}.npz"
    data = collect_data(list_dirs = sys.argv[1:])
    save_db(filename = output, data = data)

    perf = data["ns_per_day"]
    print(f"Mean: {np.nanmean(perf)} ns/day")
    print(f"Std : {np.nanstd(perf)} ns/day")
