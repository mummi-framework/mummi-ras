#!/bin/bash

echo "Run creatsims in current dir - WARMING this version has fixed input, macro model, and wrong object.data"
#echo "Run cretsims for patch $1 using variant $2"
pwd 


echo "creatsims - run for ML generated structure $1"
date
# Note, input structure needs to have the format structure_iter00_00xxxx... and be located in $MUMMI_ROOT/ML

source ~/.mummi/config.mummi.sh  # I setup new root 
source $MUMMI_APP/setup/setup.env.sh # WARNING using dev spach 
export MUMMI_SPACK_TARGET=$(spack arch --target)
# WARNING ddcmdconverter chagnes pip so need to toggle the venv 
deactivate
spack load ddcmdconverter target=$MUMMI_SPACK_TARGET
source $MUMMI_ROOT/$MUMMI_VENV/bin/activate
# Note needs to be after activate else paths get removed 
# spack load gridsim2d@v2023-03-09-precampaign3 target=$MUMMI_SPACK_TARGET
# spack load gridsim2d@v2023-03-09-precampaign3 arch=$(spack arch) # there are newer versions but likely this is the right one
#spack load gromacs@2022.6~cuda target=$MUMMI_SPACK_TARGET
# Change to CPU version 
spack load gromacs@2022.06.gpu target=$MUMMI_SPACK_TARGET
export OMP_NUM_THREADS=8

if [ -d "$1" ]; then
  echo "Dir $1 does exist assumign setup done and exiting"
  exit 1
fi
mkdir $1
cd $1


echo "creatsims - start"
date

# Copy in test protein and macro patch - not used now did before macro and ML was working 
#cp -r $MUMMI_RESOURCES/testcases/createsim/* .
#cp eA_s01_CG_merge.gro struct_createsim_init.gro # rename for creatsims 


# Run creatsims
python $MUMMI_APP/mummi_ras/transformations/createsim/createsims.py \
  --gromacs gmx --loglevel 1 --fstype simple \
  --logpath ./ -c --mpi "gmx mdrun" \
  --patch $1 --mdrunopt " -nb gpu -ntmpi 1 -ntomp $OMP_NUM_THREADS" --outpath ./ \
  --inpath $MUMMI_ROOT/ml 

# old CPU option --mdrunopt " -nt 96 -rdd 2.0 -ntomp 4 -dd 4 3 2"
## Createim - from specified GCMM patche $1 using RAS variants $2 is patch
#jsrun -n 1 --bind=none -c 24 -a 1 createsim \
#  --gromacs gmx --loglevel 1 --fstype mummi \
#  --logpath ./ -c --mpi "gmx mdrun" \
#  --patch $1 --mdrunopt " -nt 96 -rdd 2.0 -ntomp 4 -dd 4 3 2" --outpath ./ \
#  --inpath '/p/gpfs1/splash/campaign4/bootstrap-macro/good-sign/' --simnum $2 >> createsims.out  2>&1

echo "creatsims - end"
date

