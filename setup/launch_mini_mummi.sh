#!/usr/bin/env bash
# ------------------------------------------------------------------------------
if [ -f "~/.bash_profile" ]; then
  source ~/.bash_profile
fi
USRNAME=`whoami`

export MUMMI_MINI=1

echo ""
echo "(`hostname`: `date`) --> Launching Mini MuMMI for (`whoami`)"

$MUMMI_APP/setup/launch_mummi.sh "$1"