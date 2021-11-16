# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
from typing import List, Dict
import time
import yaml

MODIFY_METHODS = ['add', 'remove']
JOB_TYPES = ['createsim', 'cg', 'backmapping', 'aa']
LIST_TYPES = ['queued', 'running']
CHKPT_FILE_PREFIX = 'wfmanager'

class Checkpoint():
    def __init__(self):
        self.checkpoint = {
            'flux': '',
            'iterCounterPC': 0,
            'iterCounterWF': 0,
            'patchCounter': 0,
            'ts': time.strftime('%Y-%m-%d %H:%M:%S'),
            'jobs_aa': {
                'type': 'aa',
                'jobCnt': 0,
                'nqueued': 0,
                'nrunning': 0,
                'queued': [],
                'running': {},
            },
            'jobs_backmapping': {
                'type': 'backmapping',
                'jobCnt': 0,
                'nqueued': 0,
                'nrunning': 0,
                'queued': [],
                'running': {},
            },
            'jobs_cg': {
                'type': 'cg',
                'jobCnt': 0,
                'nqueued': 0,
                'nrunning': 0,
                'queued': [],
                'running': {},
            },
            'jobs_createsim': {
                'type': 'createsim',
                'jobCnt': 0,
                'nqueued': 0,
                'nrunning': 0,
                'queued': [],
                'running': {},
            },
        }


    def load(self, chkpt_fpath: str):
        self.checkpoint = yaml.load(open(chkpt_fpath), Loader=yaml.UnsafeLoader)
        return self


    def dump(self, out_filepath: str):
        yaml.dump(self.checkpoint, open(out_filepath, 'w'), Dumper=yaml.Dumper, sort_keys=False)
        return self


    def modify(self, patches: List[str], method: str, job_type: str, list_type: str, **kwargs):
        assert method in MODIFY_METHODS
        assert job_type in JOB_TYPES
        assert list_type in LIST_TYPES

        modify_action = {
                ("add", "queued"):     _add_to_list_unique,
                ("add", "running"):    _add_to_running_dict,
                ("remove", "queued"):  _remove_from_list,
                ("remove", "running"): _remove_from_running_dict,
            }

        # modify job list and job count
        jobs_to_modify = self.checkpoint['jobs_' + job_type]
        modified_jobs = modify_action[(method, list_type)](jobs_to_modify[list_type], patches)
        jobs_to_modify[list_type] = modified_jobs
        jobs_to_modify['n' + list_type] = len(modified_jobs)
        return self


def _add_to_list_unique(original_list: List[str] , add_list: List[str]) -> List[str]:
    return list({**dict.fromkeys(original_list), **dict.fromkeys(add_list)})


def _remove_from_list(original_list: List[str], remove_list: List[str]) -> List[str]:
    return list(filter(lambda item: item not in remove_list, original_list ))


def _add_to_running_dict(original_dict: Dict[int, List[str]], add_list: List[str]) -> Dict[int, List[str]]:

    current_patches = [item[0] for item in original_dict.values() if item]

    i = 0
    for item in add_list:
        if item not in current_patches:
            while i in original_dict: i += 1
            original_dict[i] = [item]
            i += 1
    return original_dict


def _remove_from_running_dict(original_dict: Dict[int, List[str]], remove_list: List[str]) -> Dict[int, List[str]]:

    # create reverse dict (i.e {patch_id: job_id})
    reverse_dict = { patch_list[0]: job_id for job_id, patch_list in original_dict.items() if patch_list}

    # remove patches
    for patch_id in remove_list:
        job_id = reverse_dict.get(patch_id, None)
        original_dict.pop(job_id, None)
    
    return original_dict


if __name__ == "__main__":
    import sys 
    from os.path import dirname, abspath, join, isdir
    import argparse
    import yaml

    # get config from argparse
    parser = argparse.ArgumentParser(prog='update_checkpoint')
    subparsers = parser.add_subparsers(help='sub-command help', dest='method')

    create_parser = subparsers.add_parser('create')
    create_parser.add_argument('chkpt_dirpath'),

    modify_parser = subparsers.add_parser('modify')
    modify_parser.add_argument('-i', '--inplace', action="store_true"),
    modify_parser.add_argument('chkpt_path', help="filepath to .chk.yml file"),
    modify_parser.add_argument('method', choices=MODIFY_METHODS, help="method")
    modify_parser.add_argument('job_type', choices=JOB_TYPES, help="job type")
    modify_parser.add_argument('list_type', choices=LIST_TYPES, help="list type")
    modify_parser.add_argument('patch_filepath')
    config = vars(parser.parse_args(sys.argv[1:]))

    # out base filename
    out_chkpt_fname = f"{CHKPT_FILE_PREFIX}_{time.strftime('%Y%m%d-%H%M%S')}.chk.yml"
        
    # CREATE skeleton checkpoint file
    if config['method'] == 'create':
        dirpath = abspath(config['chkpt_dirpath'])

        if isdir(dirpath):
            Checkpoint().dump(join(dirpath, out_chkpt_fname))

        else:
            print(f'{dirpath} is not a directory', file=sys.stderr)
            exit(1)
    
    # MODIFY existing checkpoint file
    if config['method'] in ['add', 'remove']:

        # out filepath 
        in_chkpt_fpath = config['chkpt_path']
        chkpt_dir = abspath(dirname(in_chkpt_fpath))
        out_chkpt_fpath = in_chkpt_fpath if config['inplace'] else join(chkpt_dir, out_chkpt_fname )

        # read in checkpoint and patch files
        patches    = [patch.strip() for patch in open(config['patch_filepath']) if patch.strip()]

        Checkpoint()\
            .load(in_chkpt_fpath)\
            .modify(patches, **config)\
            .dump(out_chkpt_fpath)

