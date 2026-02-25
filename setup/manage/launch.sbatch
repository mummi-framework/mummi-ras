#!/usr/bin/env bash

if [[ "$#" -lt 2 || "$#" -gt 3 ]]; then
  echo "Usage: $0 NODES WALLTIME (for example, 12:00:00) JOBID (Slurm Job ID for dependent JOB)" >&2
  exit 1
fi

## For round 1
export MUMMI_CONF="/lustre/orion/lrn005/proj-shared/mummi-root-c3/config.round1.mummi.sh"

if [ -z "$MUMMI_CONF" ]; then
  mummi_config=$HOME/.mummi/config.mummi.sh
  echo "(`hostname`: `date`) --> \$MUMMI_CONF was not set, defaulting to config in \$HOME directory: $mummi_config" >&2
else
  mummi_config="$MUMMI_CONF"
  echo "(`hostname`: `date`) --> using \$MUMMI_CONF=$MUMMI_CONF for config " >&2
fi

if [ ! -f $mummi_config ]; then
  echo "(`hostname`: `date`) --> ERROR: Could not find ($mummi_config)" >&2
  exit 1
fi
source $mummi_config

NODES="$1"
TIMEOUT="$2" # for ex. 12:00:00
#QUEUE="extended" # max 64 nodes for 24 hours
QUEUE="batch"
DEP="$3"
#EMAIL="ingolfsson1@llnl.gov"

echo "(`hostname`: `date`) --> Job with $NODES nodes running for $TIMEOUT time will be scheduled on queue $QUEUE."
echo "(`hostname`: `date`) --> Email will be sent to $EMAIL."
echo "(`hostname`: `date`) --> MUMMI_ROOT is $MUMMI_ROOT."

if [ -z "${DEP}" ]; then
  sbatch -N $NODES -A lrn005 -J mummi-c3 -o %x-%j.out -t $TIMEOUT -p $QUEUE --mail-user=$EMAIL --mail-type=BEGIN --mail-type=END --mail-type=FAIL $MUMMI_APP/setup/launch_mummi.sh $NODES
else
  echo "(`hostname`: `date`) --> Job will be executed after job ${DEP} is done."
  sbatch -N $NODES -A lrn005 -J mummi-c3 -o %x-%j.out -t $TIMEOUT -p $QUEUE --mail-user=${EMAIL} --mail-type=BEGIN --mail-type=END --mail-type=FAIL --dependency=afterany:${DEP} $MUMMI_APP/setup/launch_mummi.sh $NODES
fi

echo "===================="
echo -e "To monitor your job\n\tsqueue -u \$USER --start"
echo -e "\tsacct -j <JOB_ID> --format=\"JobID,JobName,Partition, NNodes,AllocCPUS, State,FailedNode,Start,WorkDir,User\""
