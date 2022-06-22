#!/bin/bash
# ------------------------------------------------------------------------------

spack load -r mummi

spack unload mummi
spack unload hwloc~cuda

spack load -r hwloc+cuda
spack load -r py-maestrowf@1.1.9
spack load -r flux-sched@0.16.0

spack load -r py-pytaridx@0.0.4
spack load -r py-h5py
spack load -r py-keras

# ------------------------------------------------------------------------------

