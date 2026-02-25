#!/usr/bin/env bash
# ------------------------------------------------------------------------------
if [ -f "~/.bash_profile" ]; then
  source ~/.bash_profile
fi
USRNAME=`whoami`

echo ""
echo "(`hostname`: `date`) --> Launching MuMMI for (`whoami`)"

# ------------------------------------------------------------------------------
# input parameters
# ------------------------------------------------------------------------------
mummi_nnodes=$1
mummi_feedback_nnodes=1
mummi_macro_nnodes=0
mummi_redis_nnodes=0
mummi_mlserver_nnodes=1
mummi_ucg_nnodes=0

if [[ $mummi_ucg_nnodes -gt 0 ]] && [[ $mummi_mlserver_nnodes -gt 0 ]]; then
  echo -n "(`hostname`: `date`) --> ERROR: For now we cannot run UCG and ML at the same time. "
  echo "mummi_ucg_nnodes = $mummi_ucg_nnodes and mummi_mlserver_nnodes = $mummi_mlserver_nnodes."
  exit 1
fi

re='^[0-9]+$'
if ! [[ $mummi_nnodes =~ $re ]] ; then
  echo "(`hostname`: `date`) --> ERROR: mummi_nnodes=$mummi_nnodes is not a number" >&2
  exit 1
fi

# ------------------------------------------------------------------------------
# load the mummi config from ($HOME/.mummi/config.mummi.sh)
# ------------------------------------------------------------------------------
if [ -z "$MUMMI_CONF" ]; then
  mummi_config=$HOME/.mummi/config.mummi.sh
  echo "(`hostname`: `date`) --> \$MUMMI_CONF was not set, defaulting to config in \$HOME directory: $mummi_config" >&2
else
  mummi_config="$MUMMI_CONF"
  echo "(`hostname`: `date`) --> using \$MUMMI_CONF=$MUMMI_CONF for config " >&2
fi


# mummi_config=$HOME/.mummi/config.mummi.sh
if [ ! -f $mummi_config ]; then
  echo "(`hostname`: `date`) --> ERROR: Could not find ($mummi_config)" >&2
  exit 1
fi
source $mummi_config
echo "(`hostname`: `date`) --> Loading MuMMI configuration ($mummi_config)"

# ------------------------------------------------------------------------------
# apply an external lock on the production run
# so we can control no job runs before we want it to
# ------------------------------------------------------------------------------
if false; then
  mummi_lockfile=$MUMMI_ROOT/mummi.lock
  if [ -f $mummi_lockfile ]; then
    echo "(`hostname`: `date`) --> ERROR: Lock file ($mummi_lockfile) found" >&2
    #exit 1
  fi
  touch $mummi_lockfile
fi

mkdir -p $MUMMI_ROOT
echo "(`hostname`: `date`) --> Copying helpers scripts for MuMMI" >&2
cp -f $MUMMI_APP/setup/manage/* $MUMMI_ROOT

# ------------------------------------------------------------------------------
# load the mummi environment
# ------------------------------------------------------------------------------
source $MUMMI_APP/setup/setup.env.sh

# Test if the packages are correctly installed in the virtualenv
if false; then
  echo ''
  python3 $MUMMI_APP/tests/test-packaging.py
fi

# ------------------------------------------------------------------------------
# print the configuration
# ------------------------------------------------------------------------------
echo ''
echo '> LSF Job details'
echo '    Job ID                :' $LSB_JOBID
echo '    hostname              :' $(hostname)
echo '    pwd                   :' $PWD
echo '> MuMMI config'
echo '    MUMMI_VENV            :' $MUMMI_VENV
echo '    MUMMI_ROOT            :' $MUMMI_ROOT
echo '    MUMMI_APP             :' $MUMMI_APP
echo '    MUMMI_CORE            :' $MUMMI_CORE
echo '    MUMMI_RESOURCES       :' $MUMMI_RESOURCES
echo '> Run config'
echo '    mummi_nnodes          :' $mummi_nnodes
echo '    mummi_feedback_nnodes :' $mummi_feedback_nnodes
echo '    mummi_macro_nnodes    :' $mummi_macro_nnodes
echo '    mummi_redis_nnodes    :' $mummi_redis_nnodes
echo '    mummi_mlserver_nnodes :' $mummi_mlserver_nnodes
echo '    mummi_ucg_nnodes      :' $mummi_ucg_nnodes

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

echo ''
python3 $MUMMI_APP/mummi_ras/scripts/create_organization.py

# ------------------------------------------------------------------------------
# launch the various components of mummi
# ------------------------------------------------------------------------------
# step 1: launch flux
$MUMMI_CORE/setup/launch_flux.sh $mummi_nnodes
source $MUMMI_ROOT/flux/flux_client.env

# steps 2 and 3: launch redis and macro
# Dropping Redis support for now.
# $MUMMI_CORE/setup/launch_redis.sh $mummi_redis_nnodes
$MUMMI_APP/setup/launch_macro.sh $mummi_macro_nnodes
$MUMMI_APP/setup/launch_feedback.sh $mummi_feedback_nnodes
$MUMMI_APP/setup/launch_mlserver.sh $mummi_mlserver_nnodes
$MUMMI_APP/setup/launch_ucg.sh $mummi_ucg_nnodes

# step 4: launch workflow
if true; then
  # temporary fixes
  # todo: pass these node counts to wf through a config?
  export MUMMI_NNODES=$mummi_nnodes
  export MUMMI_FB_NNODES=$mummi_feedback_nnodes
  export MUMMI_REDIS_NNODES=$mummi_redis_nnodes
  export MUMMI_MACRO_NNODES=$mummi_macro_nnodes
  export MUMMI_MLSERVER_NNODES=$mummi_mlserver_nnodes
  export MUMMI_UCGSERVER_NNODES=$mummi_ucg_nnodes
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
  # This did not work - for some reasons LSB_JPBID is not set here  ???
  echo "bkill bsub $LSB_JOBID"
  bkill $LSB_JOBID
fi

rm $mummi_lockfile 2>/dev/null

# ------------------------------------------------------------------------------
# end of the master script
# ------------------------------------------------------------------------------
