## All Atomistic (AA) Module

### Overview

Once the backmapping calculation is successful, 
the converted AA configurations will be added into the MUMMI workflow queue
to be simulated using the AMBER MD simulation package. 
At first, the equilibrium step will be carried out, 
including the NVT and NPT simulations. 
The NVT simulation is performed for 500 ps with 
position restraints of heavy atoms of proteins and lipids. 
The NPT simulation is carried out for 1 ns without restraints. 
The production run lasts up to 50 ns. 
The MUMMI AA simulation module allocates one GPU for each simulation 
because the multi-GPU setup is inefficient due to slow communication 
across GPUs or across nodes. Thus, many single-GPU MD simulations are 
run in parallel to achieve better GPU utilization. 
The AA module starts/restarts and monitors MD simulation runs on GPU. 

When an MD simulation is running on GPU, 
a Python-based AA analysis module (`aa_analysis`) is executed on CPU that 
analyzes new AA trajectory snapshots as soon as they are generated.
The module can analyze both RAS and RAS+RAF simulations.
As a summary, the online AA analysis methods include:

1. Capture protein coordinates after a periodic boundary condition (PBC) transform
2. Calculate the RMSD for various protein selections
3. Calculate end-to-end distance for protein selections
4. Calculate tilt/rotation/z-depth of RAS+RAF relative to the membrane

The results are sent for feedback and success/failure is reported.

Below is an example, boilerplate run command for the AA analysis module:

```
python aa_analysis.py --simname pfpatch_test01 --trajname md.1.mdcrd
--toponame gro2amber.gro --fcount 0 --step 1 --locpath /aa_traj_files
--outpath /aa_output --noninteractive
```
