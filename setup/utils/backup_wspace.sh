#!/usr/bin/env bash
# ------------------------------------------------------------------------------

pushd $MUMMI_ROOT/workspace 
timestamp=`date +%Y-%m-%d_%H-%M-%S`

d=run_backup_$timestamp
mkdir $d
echo "backup in ${d}"

mv *flux.sh $d/   2>/dev/null
mv *log $d/       2>/dev/null
mv *.err $d/      2>/dev/null
mv *.out $d/      2>/dev/null

cp jobtracker.history.csv $d/  2>/dev/null
cp wfmanager.chk $d/           2>/dev/null
mv wfmanager.chk.20210* $d/    2>/dev/null

if true; then
    # also backup the ml chkpoints
    mkdir -p $d/ml/macro
    mkdir -p $d/ml/cg
    cp $MUMMI_ROOT/ml/macro/Sampler*  $d/ml/macro/  2>/dev/null
    cp $MUMMI_ROOT/ml/cg/Sampler* $d/ml/cg/         2>/dev/null
fi

echo "backup finished"
popd

# ------------------------------------------------------------------------------
