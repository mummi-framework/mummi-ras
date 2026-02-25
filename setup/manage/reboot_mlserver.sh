#!/usr/bin/env bash
# ------------------------------------------------------------------------------
if [ -f "~/.bash_profile" ]; then
  source ~/.bash_profile
fi
USRNAME=`whoami`

echo ""
echo "(`hostname`: `date`) --> Launching MuMMI ML server only for (`whoami`)"

# ------------------------------------------------------------------------------
# input parameters
# ------------------------------------------------------------------------------
#mummi_nnodes=$1
mummi_mlserver_nnodes=1
re='^[0-9]+$'
if ! [[ $mummi_mlserver_nnodes =~ $re ]] ; then
  echo "(`hostname`: `date`) --> ERROR: mummi_mlserver_nnodes=$mummi_nnodes is not a number" >&2
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


if [ ! -f $mummi_config ]; then
  echo "(`hostname`: `date`) --> ERROR: Could not find ($mummi_config)" >&2
  exit 1
fi
source $mummi_config

# ------------------------------------------------------------------------------
# load the mummi environment
# ------------------------------------------------------------------------------
source $MUMMI_APP/setup/setup.env.sh
source $MUMMI_ROOT/flux/flux_client.env

# Just a safety measure
wfm=$(flux jobs --count=80000 | grep mlserver)
if [[ "$?" -eq 0 ]]; then
  echo "(`hostname`: `date`) --> A ML server already runs. Abort"
  echo "(`hostname`: `date`) $wfm"
  exit 1
fi

$MUMMI_APP/setup/launch_mlserver.sh $mummi_mlserver_nnodes