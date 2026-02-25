#!/usr/bin/env bash
# ------------------------------------------------------------------------------

pushd $MUMMI_ROOT/workspace > /dev/null 2>&1
timestamp=`date +%Y-%m-%d_%H-%M-%S`

d=run_backup_$timestamp
mkdir $d
echo "(`hostname`: `date`) --> backup in ${d}"

mv *flux.sh $d/   2>/dev/null
mv *log $d/       2>/dev/null
mv *.err $d/      2>/dev/null
mv *.out $d/      2>/dev/null

mv createsim $d/      2>/dev/null
mv cg $d/             2>/dev/null
mv aa $d/             2>/dev/null
mv backmapping $d/    2>/dev/null

cp jobtracker.history.csv $d/  2>/dev/null
cp wfmanager.chk $d/           2>/dev/null
mv wfmanager.chk.bak.* $d/    2>/dev/null

if true; then
    # also backup the ml chkpoints
    mkdir -p $d/ml/macro
    mkdir -p $d/ml/cg
    cp $MUMMI_ROOT/ml/macro/Sampler*  $d/ml/macro/  2>/dev/null
    cp $MUMMI_ROOT/ml/cg/Sampler* $d/ml/cg/         2>/dev/null
fi

echo "(`hostname`: `date`) --> backup finished"
popd > /dev/null 2>&1

# ------------------------------------------------------------------------------
