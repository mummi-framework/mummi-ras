#!/usr/bin/env bash
# ------------------------------------------------------------------------------

spack load py-pytest
spack load py-psutil
spack load py-filelock
spack load py-pyyaml
spack load py-numpy
spack load py-pytaridx
spack load py-redis
spack load py-maestrowf
spack load flux-sched
spack load py-cryptography@36.0.1
spack load py-scipy
spack load py-parmed
spack load py-matplotlib
spack load py-mdanalysis-mummi
spack load faiss
spack load py-keras
spack load py-h5py@2.8.0~mpi
spack load swig

spack load py-pip
# this is preventing proper functioning of ssh (on summit)
spack unload openssl

# ------------------------------------------------------------------------------
