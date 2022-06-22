#!/usr/bin/env bash

# ------------------------------------------------------------------------------
function mummi_add2path() {

  if [ -n $PATH ] ; then
    if ( echo $PATH | grep -q $1 ) ; then
      :
    else
      export PATH=$1:$PATH
    fi
  else
    export PATH=$1
  fi
}

function mummi_add2pythonpath() {

  if [ -n $PYTHONPATH ] ; then
    if ( echo $PYTHONPATH | grep -q $1 ) ; then
      :
    else
      export PYTHONPATH=$1:$PYTHONPATH
    fi
  else
    export PYTHONPATH=$1
  fi
}

function mummi_add2ldlibpath() {

  if [ -n $LD_LIBRARY_PATH ] ; then
    if ( echo $LD_LIBRARY_PATH | grep -q $1 ) ; then
      :
    else
      export LD_LIBRARY_PATH=$1:$LD_LIBRARY_PATH
    fi
  else
    export LD_LIBRARY_PATH=$1
  fi
}

# ------------------------------------------------------------------------------
