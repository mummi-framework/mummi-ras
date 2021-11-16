## Backmapping module

*Note*: This module is currently incomplete as a subset of files are stilll awaiting release permissions.

### Backmapping protocol

Due to the limitations of the CG model using the Martini force field, 
it is converted to the AA model using the CHARMM36 force field 
by the backmapping module. The backmapping procedure is implemented 
in the Python code in conjunction with Shell scripts, which have
been deposited into MUMMI_RESOURCE repo. This procedure retrieves 
an ML-selected  snapshot from the ddcMD trajectory by using 
our customized MDAnalysis Tool with the ddcMD trajectory support. 
The CG model is converted to the AA model using a modified version 
of the backward Python tool. Initial 3-step energy minimization is 
performed without nonbonded interactions between lipids and with 

* no distance restraints
* regular distance restraints
* extra strong distance restraints. 

It is followed by second energy minimization with all nonbonded 
interactions. Short position-restrained MD simulations are performed 
with increasing timestep from 0.2, 0.5, 1, to 2 fs. 
Both energy minimization and MD simulation are carried out 
on 18 designated CPU cores using GROMACS. The structure is checked, 
and stereochemistry is fixed. A final MD simulation with 
velocity Langevin dynamics without distance restraints 
is performed with the standard (unmodified) potential. 
If it is not successful, it will repeat the backmapping procedure 
up to 3 times. 

We use AMBER MD engine for the AA MD production run 
since the AMBER GPU code yields better performance. 
Thus, the input files are converted from GROMACS to AMBER 
using ParmEd. A time step of 4 fs is set for the AMBER MD simulations. 

### Python codes and Shell scripts

The backmapping modeule consists of two Python scirpts, 
backmapping.py and grotop.py, and a series of Shell scripts:

```buildoutcfg
|-- 1-EMa.mdp
|-- 1-EMb.mdp
|-- 1-EMc.mdp
|-- 2-EM.mdp
|-- 3-MD-0.0002.mdp
|-- 3-MD-0.0005.mdp
|-- 3-MD-0.001.mdp
|-- 3-MD-0.002.mdp
|-- CHL1.itp
|-- CHL1_CN_MARTINI_BACKMAPPING_ENFORCE.itp
|-- DIPE.itp
|-- DIPE_CN_MARTINI_BACKMAPPING_ENFORCE.itp
|-- Mapping
|   |-- CHL1.charmm36.map
|   |-- DIPE.charmm36.map
|   |-- PAPC.charmm36.map
|   |-- PAPS.charmm36.map
|   |-- POPC.charmm36.map
|   |-- POPE.charmm36.map
|   |-- POPS.charmm36.map
|   |-- SAPI.charmm36.map
|   |-- SSM.charmm36.map
|   |-- __init__.py
|   |-- __init__.pyc
|   |-- __pycache__
|   |   `-- __init__.cpython-36.pyc
|   |-- ace.charmm36.map
|   |-- ala.charmm36.map
|   |-- arg.charmm36.map
|   |-- asn.charmm36.map
|   |-- asp.charmm36.map
|   |-- cgw.charmm36.map
|   |-- cyf.charmm36.map
|   |-- cym.charmm36.map
|   |-- cyp.charmm36.map
|   |-- cys.charmm36.map
|   |-- gln.charmm36.map
|   |-- glu.charmm36.map
|   |-- gly.charmm36.map
|   |-- gtp.charmm36.map
|   |-- his.charmm36.map
|   |-- ile.charmm36.map
|   |-- leu.charmm36.map
|   |-- lys.charmm36.map
|   |-- met.charmm36.map
|   |-- mg.charmm36.map
|   |-- nma.charmm36.map
|   |-- phe.charmm36.map
|   |-- pro.charmm36.map
|   |-- ser.charmm36.map
|   |-- thr.charmm36.map
|   |-- trp.charmm36.map
|   |-- tyr.charmm36.map
|   |-- val.charmm36.map
|   `-- zn.charmm36.map
|-- PAPC.itp
|-- PAPC_CN_MARTINI_BACKMAPPING_ENFORCE.itp
|-- PAPS.itp
|-- PAPS_CN_MARTINI_BACKMAPPING_ENFORCE.itp
|-- POPC.itp
|-- POPC_CN_MARTINI_BACKMAPPING_ENFORCE.itp
|-- POPE.itp
|-- POPE_CN_MARTINI_BACKMAPPING_ENFORCE.itp
|-- POPS.itp
|-- POPS_CN_MARTINI_BACKMAPPING_ENFORCE.itp
|-- README
|-- SAPI.itp
|-- SAPI_CN_MARTINI_BACKMAPPING_ENFORCE.itp
|-- SSM.itp
|-- SSM_CN_MARTINI_BACKMAPPING_ENFORCE.itp
|-- ______TODO______
|-- backward.py
|-- buildNewDistanceMatrix_ras.sh
|-- buildNewDistanceMatrix_ras4A.sh
|-- buildNewDistanceMatrix_ras4Acraf.sh
|-- buildNewDistanceMatrix_rascraf.sh
|-- buildStaticDistanceMatrix_ras.sh
|-- buildStaticDistanceMatrix_ras4A.sh
|-- charmm36-jun2015_DOE_2016_Dec20_martinimod_March2020.ff
|   |-- atomtypes.atp
|   |-- cmap.itp
|   |-- ffbonded.itp
|   |-- ffnonbonded.itp
|   |-- forcefield.doc
|   |-- forcefield.itp
|   |-- gb.itp
|   |-- ions.itp
|   |-- merged.arn
|   |-- merged.c.tdb
|   |-- merged.hdb
|   |-- merged.n.tdb
|   |-- merged.r2b
|   |-- merged.rtp
|   |-- merged.vsd
|   |-- nbfix.itp
|   |-- spc.itp
|   |-- spce.itp
|   |-- tip3p.itp
|   |-- tip4p.itp
|   `-- watermodels.dat
|-- charmm36-jun2015_DOE_2016_Dec20_martinimod_March2020_addCYP_April2021.ff
|   |-- atomtypes.atp
|   |-- cmap.itp
|   |-- ffbonded.itp
|   |-- ffnonbonded.itp
|   |-- forcefield.doc
|   |-- forcefield.itp
|   |-- gb.itp
|   |-- ions.itp
|   |-- merged.arn
|   |-- merged.c.tdb
|   |-- merged.hdb
|   |-- merged.n.tdb
|   |-- merged.r2b
|   |-- merged.rtp
|   |-- merged.vsd
|   |-- nbfix.itp
|   |-- spc.itp
|   |-- spce.itp
|   |-- tip3p.itp
|   |-- tip4p.itp
|   `-- watermodels.dat
|-- dorun_campaign3_global.sh
|-- dorun_campaign3_local.sh
|-- dorun_campaign3_local.sh-orig
|-- extractSecondaryStructureFromRemark.sh
|-- fixup.mdp
|-- gro2amber.sh
|-- langevin.mdp
|-- master_make_new_distance_restraint_files.sh
|-- mdanalysis_chiral_cistrans.py
|-- mod_topandpdb_campaign3.sh
|-- neworder_protein_rascraf.itp
|-- origFranken_ras.pdb
|-- origFranken_ras4A.pdb
|-- origFranken_ras4Acraf.pdb
|-- origFranken_rascraf.pdb
|-- protein_ras.itp
|-- protein_ras4A.itp
|-- protein_ras4A_CN_MARTINI_BACKMAPPING_ENFORCE.itp
|-- protein_ras4A_CN_MARTINI_BACKMAPPING_ENFORCE_specialProteinOmega.itp
|-- protein_ras4Acraf.itp
|-- protein_ras4Acraf_CN_MARTINI_BACKMAPPING_ENFORCE.itp
|-- protein_ras4Acraf_CN_MARTINI_BACKMAPPING_ENFORCE_specialProteinOmega.itp
|-- protein_ras4Acraf_reorder.itp
|-- protein_ras_CN_MARTINI_BACKMAPPING_ENFORCE.itp
|-- protein_ras_CN_MARTINI_BACKMAPPING_ENFORCE_specialProteinOmega.itp
|-- protein_rascraf.itp
|-- protein_rascraf_CN_MARTINI_BACKMAPPING_ENFORCE.itp
|-- protein_rascraf_CN_MARTINI_BACKMAPPING_ENFORCE_specialProteinOmega.itp
|-- protein_rascraf_reorder.itp
|-- sinceCG.sh
|-- template.top
|-- toalignrascrafonly.ndx
`-- toalignrasonly.ndx
```

The backmapping.py is the main application and the main protocol is
implemented in the `sinceCG.sh` script for the backmapping calculations.

### Testing script

The following is a testing script for the backmapping module under 
the MuMMI farmework. It assumes that the input files (`topol.tpr`, 
`system.top`, `positions.tar`) for the backmapping calculation are in 
`$MUMMI_ROOT/sims-cg/${patch_id}`. The `backmapping.py` in the following
script takes patch id (pfpatch_000000000012) and frame id(1325000) 
by using options `-p` and `-f` respectively. User can specify 
where to perform the calculation by option `--outlocal`. If not
specified, the calculation will be performed in a subdirectory of
the `/tmp/mummi` by default.

```buildoutcfg
#!/bin/bash

#spack environment setup.
source ~/.mummi/config.mummi.sh
source $MUMMI_APP/setup/setup.env.sh

spack load -r py-mdanalysis-mummi@mda_1.0.1_ddcmd
spack load py-tqdm
spack load -r py-parmed

spack load -r gromacs@2019.6.mummifix +double
spack load -r gromacs@2019.6.mummifix ~double

python $MUMMI_APP/mummi_ras/transformations/backmapping/backmapping.py -p pfpatch_000000000012 -f 1325000 --outlocal /g/g92/zhang30/W/MuMMI/scratch/backmap_h2
```