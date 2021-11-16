## Createsims Module 
The MuMMI Creatsims Module constructs a CG Martini simulation from a macro model patch. It supports constructions of CG simulations from a few different versions of macro model input data as well as simulation construction based on manual protein definition and random lipid placement. After initial simulation construction, Creatsims runs initial simulation equilibrium and converts final production ready input files from GROMACS (www.gromacs.org) for ddcMD ([github.com/LLNL/ddcMD](https://github.com/LLNL/ddcMD)) format.

At a high-level, Createsims follows these steps:
1.	Retrieve the selected macro model patch and parse to simulation frame of reference
2.	Convert lipid concentration individual lipids
3.	Selected and place proteins from input libraries of protein configurations
4.	Construct the system initial coordinates using a modified version of the insane.py builder
5.	Energy minimizes, pull proteins to membrane, and run initial equilibrium of system using a series of EM and MD simulations run with GROMACS
6.	Convert final GROMACS inputs to ddcMD using ddcMDconverter ([github.com/LLNL/ddcMDconverter](https://github.com/LLNL/ddcMDconverter)).

### Execution Requirements
The Creatsims Module is Python based and itself runs on 1 CPU, but for the GROMACS equilibrium runs managed by Creatsims can run on different compute resources, current default in MuMMI is 24 CPUs 4 threads each (96 threads total).  

### Example of standalone test launch 
```
  <source MuMMI environment>
  createsims
    --gromacs gmx --loglevel 1 --fstype taridx 
    --logpath ./ -c --mpi "gmx mdrun" 
    --patch <name of patch> --mdrunopt " -nt 96 -rdd 2.0 -ntomp 4 -dd 4 3 2" 
    --inpath '<macro model path>' --outpath ./ --simnum test_run_1
```