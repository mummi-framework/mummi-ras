#!/usr/bin/env bash

job_createsim=$(grep -A 3 jobs_createsim workspace/wfmanager.chk)
job_cg=$(grep -A 3 jobs_cg workspace/wfmanager.chk)
job_aa=$(grep -A 3 jobs_aa workspace/wfmanager.chk)
job_bm=$(grep -A 3 jobs_backmapping workspace/wfmanager.chk)

echo -e "> MuMMI run in $MUMMI_ROOT"
echo -e " > Createsim    running=$(echo $job_createsim | cut -d' ' -f7) queued=$(echo $job_createsim | cut -d' ' -f5)"
echo -e " > CG           running=$(echo $job_cg | cut -d' ' -f7) queued=$(echo $job_cg | cut -d' ' -f5)"
echo -e " > AA           running=$(echo $job_aa | cut -d' ' -f7) queued=$(echo $job_aa | cut -d' ' -f5)"
echo -e " > Backmapping  running=$(echo $job_bm | cut -d' ' -f7) queued=$(echo $job_bm | cut -d' ' -f5)"
