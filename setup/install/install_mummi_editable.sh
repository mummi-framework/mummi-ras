#!/bin/bash

# disabling sourcing of this script!
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
  echo "> ERROR: Please do not source this script!"
  return
fi

if [ -n "$VIRTUAL_ENV" ] ; then
  :
else
  echo '> ERROR: Please do not attempt to install editable MuMMI outside your local virtual environment.
        This could cause issues for everyone by accidently polluting spack environment.'
  exit
fi

# ------------------------------------------------------------------------------
echo ''
echo '>' `date`
echo "> This script (re)installs editable MuMMI for (`whoami`) in virtual environment ($VIRTUAL_ENV)!"
echo ''

# spack find py-pip
# if [ $? -eq 0 ]; then
#   echo "Unloading Spack pip"
#   spack unload py-pip
# fi

TIMEOUT=600
# ------------------------------------------------------------------------------
echo "> Using pip: `command -v pip`"
# Problematic on Frontier where compute nodes don't have access to internet.
# echo "> Upgrading pip"
# pip install --require-virtualenv --default-timeout=$TIMEOUT --upgrade pip
echo '> Removing existing build of MuMMI'
pip uninstall --yes mummi_core mummi_ras

# ------------------------------------------------------------------------------
#spack load swig

export CC=`command -v mpicc`
export CXX=`command -v mpic++`

# ------------------------------------------------------------------------------
echo '> Building MuMMI_core'
#cd $MUMMI_CORE
rm -rf $MUMMI_CORE/build 2>/dev/null
rm -rf $MUMMI_CORE/mummi_core.egg-info 2>/dev/null
echo "Current pip: `command -v pip`"
pip install --require-virtualenv --default-timeout=$TIMEOUT -e $MUMMI_CORE --prefix $VIRTUAL_ENV

echo '> Building MuMMI_ras'
#cd $MUMMI_APP
rm -rf $MUMMI_APP/build 2>/dev/null
rm -rf $MUMMI_APP/mummi_ras.egg-info 2>/dev/null
pip install --require-virtualenv --default-timeout=$TIMEOUT -e $MUMMI_APP --prefix $VIRTUAL_ENV

# ------------------------------------------------------------------------------
#spack unload swig

# ------------------------------------------------------------------------------
HOST=`hostname --long`
if [[ $HOST == *summit* ]]; then

    ENV_PYTHON="#\!/bin/env python"
    HEADER_FILES=( $( grep -I "opt/spack" $VIRTUAL_ENV/bin/* | awk -F":" ' { print $1 } ' ))
    HEADER_FILES+=( $( grep -I "$VIRTUAL_ENV" $VIRTUAL_ENV/bin/* | awk -F":" ' { print $1 } ' ))
    for f in "${HEADER_FILES[@]}"
    do
        echo "Replacing python header in (${f})"
        sed -i.bak --expression "1s@^.*@$ENV_PYTHON@" $f
    done
fi

echo ''
echo '> Installed MuMMI in ('$VIRTUAL_ENV')'
# ------------------------------------------------------------------------------

