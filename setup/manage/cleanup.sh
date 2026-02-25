#!/usr/bin/env bash

RM="rm -v -r -f"
CP="cp -p"

$RM *.out
$RM config.yaml
$RM mlserver/*.log
$RM ucgserver/*.log
$RM monitoring/*
$RM workspace/create_sim*
$RM workspace/createsim-*
$RM workspace/cganalysis-*
$RM workspace/wfmanager.chk.bak.*
$RM workspace/wfmanager*.log
$RM workspace/wfmanager*.err
$RM workspace/*.csv
$RM workspace/run_backup*
$RM flux/*
$RM backup/*
$RM createsim/*
$RM cg/*
