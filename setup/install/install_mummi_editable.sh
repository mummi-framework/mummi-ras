#!/bin/bash
# ------------------------------------------------------------------------------


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

# ------------------------------------------------------------------------------
echo "> Using pip: `which pip`"
echo '> Removing existing build of MuMMI'

pip uninstall mummi_core mummi_ras

# ------------------------------------------------------------------------------
export CC=`which mpicc`
export CXX=`which mpic++`
spack load swig@4.0.2/lk5op3c

# ------------------------------------------------------------------------------
pushd $MUMMI_CORE >/dev/null

echo "> Building MuMMI_core (`pwd`)"

rm -rf build 2>/dev/null
rm -rf mummi_core.egg-info 2>/dev/null
pip install -v -e .

popd >/dev/null


pushd $MUMMI_APP >/dev/null

echo "> Building MuMMI_ras (`pwd`)"
rm -rf build 2>/dev/null
rm -rf mummi_ras.egg-info 2>/dev/null
pip install -v -e .

popd >/dev/null

# ------------------------------------------------------------------------------
spack unload swig


# ------------------------------------------------------------------------------
if [[ $MUMMI_HOST == summit ]]; then
    sh $MUMMI_APP/setup/utils/replace_python_in_venv_bins.sh 
fi

echo ''
echo '> Installed MuMMI in ('$VIRTUAL_ENV')'
# ------------------------------------------------------------------------------

