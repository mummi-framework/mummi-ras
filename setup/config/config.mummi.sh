#!/bin/bash

# ------------------------------------------------------------------------------
# likely to be set for each run
# ------------------------------------------------------------------------------
# path where you want to run mummi (workspace)
export MUMMI_ROOT=/usr/workspace/mummiusr/mummi_demoroot_20220724

# ------------------------------------------------------------------------------
# likely to be set only once
# ------------------------------------------------------------------------------
# virtual environment name of your choice
# no need to change (unless you want to do a complete reinstall)
export MUMMI_VENV=mummi-ras-demo

# path where you clone mummi repositories
# no need to change (unless you move your repositories)
export MUMMI_RESOURCES=/usr/workspace/mummiusr/mummi_resources
export MUMMI_CORE=/usr/workspace/mummiusr/mummi_framework/mummi-core
export MUMMI_APP=/usr/workspace/mummiusr/mummi_framework/mummi-ras

# ------------------------------------------------------------------------------
