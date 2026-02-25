#!/usr/bin/env bash
#SBATCH -A lrn005
#SBATCH -J interpolator
#SBATCH -o %x-%j.out
#SBATCH -t 2:00:00
#SBATCH -p batch
#SBATCH -N 1

#SBATCH --mail-user=pottier1@llnl.gov
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

source ~/.mummi/config.mummi.sh
source $MUMMI_APP/setup/setup.env.sh

#   -f FRACTION, --fraction FRACTION
#                         Fraction of states to use
#   -lu LAMBDA_UPPER, --lambda-upper LAMBDA_UPPER
#                         Upper bound for lambda
#   -ll LAMBDA_LOWER, --lambda-lower LAMBDA_LOWER
#                         Lower bound for lambda
python3 $MUMMI_APP/tests/test-generate-interpolator.py -f "$1"
