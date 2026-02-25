#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Copyright (c) 2025, Lawrence Livermore National Security, LLC. All rights
# reserved. LLNL-CODE-827655. This work was produced at the Lawrence Livermore
# National Laboratory (LLNL) under contract no. DE-AC52-07NA27344 (Contract 44)
# between the U.S. Department of Energy (DOE) and Lawrence Livermore National
# Security, LLC (LLNS) for the operation of LLNL.  See license for disclaimers,
# notice of U.S. Government Rights and license terms and conditions.
# -----------------------------------------------------------------------------
#
# This script launches two feedback process:
#   - The createsims feedback that monitors createsims simulation status
#   - The CG2ML feedback process that converts CG frames into latent space coordinates

function version() {
  echo "$@" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }';
}

mummi_feedback_nnodes=$1
re='^[0-9]+$'
if ! [[ $mummi_feedback_nnodes =~ $re ]] ; then
  echo "(`hostname`: `date`) --> ERROR: $mummi_feedback_nnodes is not a number" >&2
  exit 1
fi

if [ $mummi_feedback_nnodes -eq 0 ]; then
  echo "(`hostname`: `date`) --> WARNING: Not Launching Feedback because mummi_feedback_nnodes = $mummi_feedback_nnodes" >&2
  exit 1
fi

# ------------------------------------------------------------------------------

# We place the feedback in the mlserver directory
mkdir -p $MUMMI_ROOT/mlserver

# ------------------------------------------------------------------------------
pushd $MUMMI_ROOT/mlserver

# CG2ML feedback

NUM_THREADS_FB_CG2ML=50 # BE CAREFUL TO NOT ALLOCATE MORE CORES THAN AVAILABLE

feedback_cg2ml_cmd="python3 $MUMMI_APP/mummi_ras/scripts/run_feedback_cg2ml.py --refresh 120 -t ${NUM_THREADS_FB_CG2ML}"

echo "(`hostname`: `date`) --> Launching Feedback CG2ML with $mummi_feedback_nnodes nodes ($NUM_THREADS_FB_CG2ML threads)"
echo "   Working directory:      `pwd`"
echo "   Feedback CG2ML command:   $feedback_cg2ml_cmd"

#if in python env, we make sure we activate it
if [[ ! -z "${VIRTUAL_ENV}" ]]; then
  echo "(`hostname`: `date`) --> Sourcing $VIRTUAL_ENV/bin/activate"
  source $VIRTUAL_ENV/bin/activate
fi

FLUX_CMD="flux"

FLUX_VERSION=$(version $(flux version | awk '/^commands/ {print $2}'))
MIN_VER_FLUX=0.55 # Below 0.55 we can use flux mini but >= that we cannot
MIN_VER_FLUX_LONG=$(version ${MIN_VER_FLUX})
# Remove leading 0 because they are interpreted as octal numbers in bash
if [[ "${FLUX_VERSION#00}" -lt "${MIN_VER_FLUX_LONG#00}" ]]; then
  FLUX_CMD="flux mini"
fi

jobid=$(${FLUX_CMD} submit -N "$mummi_feedback_nnodes" -n 1 -c $NUM_THREADS_FB_CG2ML \
    --job-name=fb-cg2ml. \
    --requires=rank:1 \
    --output=feedback-cg2ml-{{id}}.log \
    sh -c "$feedback_cg2ml_cmd")

echo "(`hostname`: `date`) --> Feedback CG2ML job (job ID: $jobid)"


NUM_THREADS_FB_CREATESIMS=6 # BE CAREFUL TO NOT ALLOCATE MORE CORES THAN AVAILABLE
feedback_createsims_cmd="python3 $MUMMI_APP/mummi_ras/scripts/run_feedback_createsims.py --refresh 60 -t ${NUM_THREADS_FB_CREATESIMS}"

echo "(`hostname`: `date`) --> Launching Feedback createsims with $mummi_feedback_nnodes nodes ($NUM_THREADS_FB_CREATESIMS threads)"
echo "   Feedback createsims command:   $feedback_createsims_cmd"

jobid=$(${FLUX_CMD} submit -N "$mummi_feedback_nnodes" -n 1 -c $NUM_THREADS_FB_CREATESIMS \
    --job-name=fb-createsims. \
    --requires=rank:1 \
    --output=feedback-createsims-{{id}}.log \
    sh -c "$feedback_createsims_cmd")

echo "(`hostname`: `date`) --> Feedback Createsims job (job ID: $jobid)"

popd > /dev/null 2>&1
# ------------------------------------------------------------------------------
