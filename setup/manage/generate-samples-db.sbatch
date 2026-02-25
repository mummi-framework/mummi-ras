#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -A lrn005
#SBATCH -J mummi-gen-struct
#SBATCH -o %x-%j.out
#SBATCH -p extended
#SBATCH -t 20:00:00

## For round 2
# export MUMMI_CONF="/lustre/orion/lrn005/proj-shared/mummi-root-c3/config.round2.mummi.sh"

export MUMMI_CONF="/ccs/home/lpottier/mummi-framework-c3/config-fbucg.mummi.sh"

source $MUMMI_CONF
source $MUMMI_APP/setup/setup.env.sh
NSAMPLES=60000

# OUT_DIR=$MUMMI_ROOT
OUT_DIR="/lustre/orion/lrn005/proj-shared/mummi-root-c3/round-2"

python3 -m torch.distributed.launch --use-env $MUMMI_APP/tests/test-generators.py -n ${NSAMPLES} -o $OUT_DIR/pre-generated-c3-r2-${NSAMPLES}.npz

chmod g+rw $OUT_DIR/pre-generated-c3-r2-${NSAMPLES}.npz

# Saving a backup just in case
mkdir -p $OUT_DIR/backup-db/
chmod g+rw $OUT_DIR/backup-db/
cp $OUT_DIR/pre-generated-c3-r2-${NSAMPLES}.npz $OUT_DIR/backup-db/

echo "Generated file in $OUT_DIR/pre-generated-c3-r2-${NSAMPLES}.npz"