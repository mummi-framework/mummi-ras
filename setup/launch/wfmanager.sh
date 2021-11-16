#!/bin/bash

pushd $MUMMI_ROOT/workspace

# Launching wfmanager
autobind-24 python $MUMMI_APP/mummi_ras/scripts/run.wfmanager.py >> wfmanager.log 2>&1
popd
