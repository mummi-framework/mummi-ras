#!/usr/bin/env bash
# ------------------------------------------------------------------------------

jobid=$1    # $LSB_JOBID
myuser=$2   # or any user name you are running as.
pdshhosts=`lsfjobs -A HOSTLIST -j $jobid | grep " $myuser " | awk '{print $NF}'`

# Check to see that it looks reasonable
echo "  > Run ($3) on all nodes of jobs ($1) for user ($2), nodes: $pdshhosts"
pdsh -R ssh -w "$pdshhosts" $3

# ------------------------------------------------------------------------------