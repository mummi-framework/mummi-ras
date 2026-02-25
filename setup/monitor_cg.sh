l#!/usr/bin/env bash

if [ -f "~/.bash_profile" ]; then
  source ~/.bash_profile
fi
USRNAME=`whoami`

echo ""
echo "(`hostname`: `date`) --> Launching CG monitoring for (`whoami`)"
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

while true; do
  echo "(`hostname`: `date`) --> waiting for 45m..."
  sleep 25m
  echo "(`hostname`: `date`) --> done waiting. ready for next check."
  $MUMMI_APP/setup/manage/check_sims_length.sh
done

