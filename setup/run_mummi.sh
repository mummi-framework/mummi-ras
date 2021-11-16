#!/usr/bin/env bash
source ~/.bash_profile


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
SETUP_CONFIG=$HOME/.mummi/config.mummi.sh
if [ ! -f "$SETUP_CONFIG" ]; then
  echo "> ERROR: Could not find $SETUP_CONFIG. This file should load necessary environment variables for MuMMI."
  exit 1
fi

source $SETUP_CONFIG

# ------------------------------------------------------------------------------
# check the input parameters
# ------------------------------------------------------------------------------
if ! [ -n "$MUMMI_ROOT" ]; then
  echo 'ERROR: please set MUMMI_ROOT.'
  exit 1
fi

re='^[0-9]+$'
if ! [[ $1 =~ $re ]] ; then
   echo "ERROR: nodes is not a number" >&2; exit 1
fi

# ------------------------------------------------------------------------------
# apply an external lock on the production run
# so we can control no job runs before we want it to
LOCK_FILE=$MUMMI_ROOT/mummi.lock
if [ -f $LOCK_FILE ]; then
  echo "> Lock file ($LOCK_FILE) found! Exiting MuMMI job."
  exit 1
fi

#LOCK_FILE=$MUMMI_ROOT/mummi.go
#$MUMMI_UTILS/wait_for_file.sh $LOCK_FILE

# ------------------------------------------------------------------------------
# set up the node counts for this run
# ------------------------------------------------------------------------------
#MUMMI_NNODES=`lsfjobs -A HOSTS -j $LSB_JOBID | grep $USRNAME | awk '{print $NF}'`
export MUMMI_NNODES=$1
export MUMMI_RESTART=yes           # yes   or   no

# ------------------------------------------------------------------------------
# launch MuMMI
# ------------------------------------------------------------------------------
USRNAME=`whoami`
echo ''
echo '----> Launching MuMMI for ('$USRNAME') on ('`hostname`')...'
echo `date`
echo ''
echo '    LSB_JOBID       :' $LSB_JOBID
echo '    MUMMI_NNODES    :' $MUMMI_NNODES
echo '    PWD             :' $PWD
echo ''
echo '    MUMMI_ROOT        :' $MUMMI_ROOT
echo '    MUMMI_REDIS_NNODES:' $MUMMI_REDIS_NNODES
echo '    MUMMI_MACRO_NNODES:' $MUMMI_MACRO_NNODES
echo ''
echo '    MUMMI_RESTART   :' $MUMMI_RESTART
echo ''
echo '    MUMMI_VENV      :' $MUMMI_VENV
echo '    MUMMI_APP       :' $MUMMI_APP
echo '    MUMMI_CORE      :' $MUMMI_CORE
echo '    MUMMI_RESOURCES :' $MUMMI_RESOURCES

# ------------------------------------------------------------------------------
# set up the environment (if needed)
# ------------------------------------------------------------------------------
if ! [ -n "$VIRTUAL_ENV" ] ; then
  source $MUMMI_APP/setup/setup.env.sh
else
  v="$(basename -- $VIRTUAL_ENV)"
  if ! [ "$v" = "$MUMMI_VENV" ]; then
    echo "ERROR: expected virtual environment ($MUMMI_VENV), but found ($v)."
    exit 1
  fi
fi

MUMMI_UTILS=$MUMMI_APP/setup/utils

# ------------------------------------------------------------------------------
# Check all nodes and clean /var/tmp
# ------------------------------------------------------------------------------
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


# ------------------------------------------------------------------------------
# Create the mummi root!
# ------------------------------------------------------------------------------
echo ''
echo '----> Creating MUMMI_ROOT('$MUMMI_ROOT')'
python $MUMMI_APP/mummi_ras/scripts/create_organization.py

PERM_READY_FILE=$MUMMI_ROOT/workspace/fix_perms_done.txt
if [[ $MUMMI_RESTART != 'yes' ]] || [[ ! -f $PERM_READY_FILE ]]; then
  $MUMMI_UTILS/fix_perms.sh $MUMMI_GROUP
  touch $PERM_READY_FILE
fi

# ------------------------------------------------------------------------------
# Launch flux
# ------------------------------------------------------------------------------
FLUX_ROOT=$MUMMI_ROOT/flux
FLUX_FILE=$FLUX_ROOT/flux.info
FLUX_BOOTSTRAP=$MUMMI_CORE/setup/flux/flux_bootstrap.sh

if [ -f $FLUX_FILE ] ; then
  $MUMMI_UTILS/create_backup_of_directory.sh $FLUX_ROOT
  rm -rf $FLUX_ROOT
fi

$MUMMI_CORE/setup/flux/flux_launch.sh $MUMMI_NNODES $FLUX_ROOT &

# now, wait for the flux info file
$MUMMI_UTILS/wait_for_file.sh $FLUX_FILE
if [ -f $FLUX_FILE ]; then
  echo "> Flux launched successfully."
else
  echo "> Failed to launch Flux."
  exit 1
fi

# Export
export FLUX_URI=`cat $FLUX_FILE`
export FLUX_NNODES=`flux getattr size`
echo "> FLUX_URI    :" $FLUX_URI
echo "> FLUX_SSH    :" $FLUX_SSH
echo "> FLUX_NNODES :" $FLUX_NNODES
echo " > Flux launched using jsrun"
echo " > flux modules: "
echo "`flux module list`"
echo " > version:"
echo "`flux version`"
echo "--------------------------------------------"
echo "INFO: Adjusting grace period to 2m..."
flux module reload job-exec kill-timeout=15.0m
echo "Sleeping for 10s..."
sleep 10s
unset PMI_LIBRARY

# ------------------------------------------------------------------------------
# Launch redis server
# ------------------------------------------------------------------------------
if [ $MUMMI_REDIS_NNODES -gt 0 ]; then
  
  echo '----> Launching Redis on' $MUMMI_REDIS_NNODES 'nodes'
 
  source $MUMMI_CORE/setup/redis/load_redis_env.sh
  source $MUMMI_CORE/setup/redis/start_all_redis_nodes.sh $MUMMI_REDIS_NNODES
  #sleep 10s
  #source $MUMMI_CORE/setup/redis/start_all_fakecg_nodes.sh $MUMMI_FAKECG_NNODES

  if [ $MUMMI_FAKECG_NNODES -gt 0 ]; then
    echo '-> Launching fake cg'
    for i in $(seq 1 $MUMMI_FAKECG_NNODES);
    do
      FAKECG_STDOUT="${MUMMI_REDIS_OUTPUTS}/fakecg${i}.stdout"
      echo $FAKECG_STDOUT
      flux mini submit --job-name fakecg$i --output=$FAKECG_STDOUT -n1 -N1 -c24 flux start python $MUMMI_CORE/tests/test-fake_cg.py
      # fdinatal -- commenting this out because it'll fail
    done
  fi

else
  echo '----> Disabled Redis'
fi



# sh $MUMMI_CODE/tests/jsrun-submit.sh

# ------------------------------------------------------------------------------
# Launch macro
# ------------------------------------------------------------------------------
if [ $MUMMI_MACRO_NNODES -gt 0 ]; then
   echo '----> Launching Macro on' $MUMMI_MACRO_NNODES 'nodes'
   cd $MUMMI_ROOT/macro
   $MUMMI_APP/setup/launch/launch_macro.sh
fi

# ------------------------------------------------------------------------------
# Launch MuMMI with Bash
# ------------------------------------------------------------------------------
echo '----> Launching MuMMI with shell'
cd $MUMMI_ROOT/workspace
$MUMMI_APP/setup/launch/launch_mummi.sh

# ------------------------------------------------------------------------------
# Launch MuMMI with Maestro
# ------------------------------------------------------------------------------
# echo '----> Launching MuMMI with Maestro'
#cp $MUMMI_CODE/resources/run_splash.yaml .
# cd $MUMMI_ROOT/workspace
# $MUMMI_CODE/setup/launch/launch_mummi.sh $MUMMI_CODE


# ------------------------------------------------------------------------------
# the monitor script will also check if exit is needed!
#if false; then
echo "---> Track MuMMI status"
python $MUMMI_CORE/mummi_core/scripts/monitor_mummi.py $MUMMI_APP/specs/workflow/monitor.yaml
#fi


sleep inf


if false; then

# ------------------------------------------------------------------------------
# After MuMMI - cleanup/stopping
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# end of the master script
# ------------------------------------------------------------------------------
