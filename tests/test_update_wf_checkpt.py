# Copyright (c) 2021, Lawrence Livermore National Security, LLC. All rights reserved. LLNL-CODE-827655.
# This work was produced at the Lawrence Livermore National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44) between the U.S. Department of Energy (DOE) and Lawrence Livermore National Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers, notice of U.S. Government Rights and license terms and conditions.
# ------------------------------------------------------------------------------
from mummi_ras.scripts import update_wf_chkpt


def test_add_to_list_unique():
    orig_list = ['a', 'b', 'c']
    add_list  = ['d', 'a', 'c', 'e']
    expected  = ['a', 'b', 'c', 'd', 'e']
    assert update_wf_chkpt._add_to_list_unique(orig_list, add_list) == expected 

def test_remove_from_list_unique():
    orig_list   = ['a', 'b', 'c']
    remove_list = ['d', 'a', 'c', 'e']
    expected    = ['b']
    assert update_wf_chkpt._remove_from_list(orig_list, remove_list) == expected

def test_remove_from_running_dict():
    orig_dict   = {1: ['a'], 2: ['b'], 3: ['c']}
    remove_list = ['a', 'c', 'd']
    expected    = {2: ['b']}
    assert update_wf_chkpt._remove_from_running_dict(orig_dict, remove_list) == expected 

def test_add_to_running_dict():
    orig_dict   = {1: ['a'], 2: ['b'], 3: ['c']}
    add_list    = ['a', 'c', 'd', 'e']
    expected    = {1: ['a'], 2: ['b'], 3: ['c'], 0: ['d'], 4: ['e']}
    assert update_wf_chkpt._add_to_running_dict(orig_dict, add_list) == expected 

