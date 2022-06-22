#!/usr/bin/env bash
# ------------------------------------------------------------------------------
source ~/.bash_profile
USRNAME=`whoami`

echo ""
echo "(`hostname`: `date`) --> Launching MuMMI for (`whoami`)"

# ------------------------------------------------------------------------------
# input parameters
# ------------------------------------------------------------------------------
mummi_nnodes=$1
mummi_macro_nnodes=0
mummi_redis_nnodes=0

re='^[0-9]+$'
if ! [[ $mummi_nnodes =~ $re ]] ; then
  echo "(`hostname`: `date`) --> ERROR: mummi_nnodes=$mummi_nnodes is not a number" >&2
  exit 1
fi

# ------------------------------------------------------------------------------
# load the mummi config from ($HOME/.mummi/config.mummi.sh)
# ------------------------------------------------------------------------------
mummi_config=$HOME/.mummi/config.mummi.sh
if [ ! -f $mummi_config ]; then
  echo "(`hostname`: `date`) --> ERROR: Could not find ($mummi_config)" >&2
  exit 1
fi
source $mummi_config

# ------------------------------------------------------------------------------
# apply an external lock on the production run
# so we can control no job runs before we want it to
# ------------------------------------------------------------------------------
mummi_lockfile=$MUMMI_ROOT/mummi.lock
if [ -f $mummi_lockfile ]; then
  echo "(`hostname`: `date`) --> ERROR: Lock file ($mummi_lockfile) found" >&2
  exit 1
fi

touch $mummi_lockfile

# ------------------------------------------------------------------------------
# load the mummi environment
# ------------------------------------------------------------------------------
source $MUMMI_APP/setup/setup.env.sh

# ------------------------------------------------------------------------------
# print the configuration
# ------------------------------------------------------------------------------
echo ''
echo '> LSF Job details'
echo '    Job ID             :' $LSB_JOBID
echo '    hostname           :' `hostname`
echo '    pwd                :' $PWD
echo '> MuMMI config'
echo '    MUMMI_VENV         :' $MUMMI_VENV
echo '    MUMMI_ROOT         :' $MUMMI_ROOT
echo '    MUMMI_APP          :' $MUMMI_APP
echo '    MUMMI_CORE         :' $MUMMI_CORE
echo '    MUMMI_RESOURCES    :' $MUMMI_RESOURCES
echo '> Run config'
echo '    mummi_nnodes       :' $mummi_nnodes
echo '    mummi_macro_nnodes :' $mummi_macro_nnodes
echo '    mummi_redis_nnodes :' $mummi_redis_nnodes

# ------------------------------------------------------------------------------
# Check all nodes and clean /var/tmp
# ------------------------------------------------------------------------------
if false; then
  echo ''
  echo '----> Checking the status of nodes'
  #echo '----> WARNING - removed as this one is always getting errors'
  # TODO: Likely will not work on Summit
  check_sierra_nodes
  # Worked for 600, Sun Nov11 tried 1200 x5 times - x3 errors x2 frezes

  echo ''
  echo '----> Cleaning up nodes before launch (all our relevant files on compute node)'
  $MUMMI_UTILS/run_with_pdsh.sh $LSB_JOBID $USRNAME $MUMMI_UTILS/clean_node.sh

  # Check nodes /tmp space - @TODO make this only print if problem
  echo ''
  echo '----> Checking nodes before launch (space on compute node)'
  $MUMMI_UTILS/run_with_pdsh.sh $LSB_JOBID $USRNAME $MUMMI_UTILS/check_node.sh
fi

# ------------------------------------------------------------------------------
# Prepare for the run!
# ------------------------------------------------------------------------------

# need access to mummi. but paths getting overridden? this is a temporary fix!
spack load py-pip
spack unload openssl

echo ''
python3 $MUMMI_APP/mummi_ras/scripts/create_organization.py

# -----------------------------------------------------------------------------
#PERM_READY_FILE=$MUMMI_ROOT/workspace/fix_perms_done.txt
#if [[ $MUMMI_RESTART != 'yes' ]] || [[ ! -f $PERM_READY_FILE ]]; then
#  $MUMMI_UTILS/fix_perms.sh $MUMMI_GROUP
#  touch $PERM_READY_FILE
#fi

# ------------------------------------------------------------------------------
# launch the various components of mummi
# ------------------------------------------------------------------------------

# step 1: launch flux
$MUMMI_CORE/setup/launch_flux.sh $mummi_nnodes
source $MUMMI_ROOT/flux/flux_client.env

# steps 2 and 3: launch redis and macro
$MUMMI_CORE/setup/launch_redis.sh $mummi_redis_nnodes
$MUMMI_APP/setup/launch_macro.sh $mummi_macro_nnodes

# step 4: launch workflow
if true; then
   # temporary fixes
   # todo: pass these node counts to wf through a config?
   export MUMMI_NNODES=$mummi_nnodes
   export MUMMI_REDIS_NNODES=$mummi_redis_nnodes
   export MUMMI_MACRO_NODES=$mummi_macro_nnodes
fi
$MUMMI_APP/setup/launch_workflow.sh $mummi_nnodes

# step 5: launch monitoring script
$MUMMI_CORE/setup/launch_monitor.sh $MUMMI_APP/specs/workflow/monitor.yaml


# ------------------------------------------------------------------------------
# After MuMMI - cleanup/stopping
# ------------------------------------------------------------------------------
if false; then
  # TODO: to be changed
  echo "----> Cleanup after MuMMI Run"
  flux exec -r $i pkill $FLUX_BOOTSTRAP # then we should not need bkill on launch_splash
  echo "> Flux should now stop by itself and bsub should return"

  sleep 2m

  echo ''
  echo '----> Cleaning up nodes after exiting (all our relevant files on compute node)'
  #mpirun -N 1 -n "$((SPLASH_NNODES + 1))" $MUMMI_CODE/resources/cleanup/cleanup
  $MUMMI_UTILS/run_with_pdsh.sh $LSB_JOBID $USRNAME $MUMMI_UTILS/clean_node.sh


  # Check nodes /tmp space - @TODO make this only print if problem
  echo ''
  echo '----> Checking nodes before launch (space on compute node)'
  $MUMMI_UTILS/run_with_pdsh.sh $LSB_JOBID $USRNAME $MUMMI_UTILS/check_node.sh


  # The following should not be needed!

  # Return partition
  # This did not work - for some reasion LSB_JPBID is not set here  ???
  echo "bkill bsub $LSB_JOBID"
  bkill $LSB_JOBID
fi

rm $mummi_lockfile

# ------------------------------------------------------------------------------
# end of the master script
# ------------------------------------------------------------------------------
