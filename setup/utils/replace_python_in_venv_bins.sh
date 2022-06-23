#!/bin/bash
# ------------------------------------------------------------------------------

# this script looks at all entry points and replaces the path to python
# on Summit, we have seen that python path points to spack and gets truncated

env_python="#\!/bin/env python"

pushd $VIRTUAL_ENV/bin

headers=( $( grep -I "opt/spack" * | awk -F":" ' { print $1 } ' ))
headers+=( $( grep -I "$VIRTUAL_ENV" * | awk -F":" ' { print $1 } ' ))
for f in "${headers[@]}"
do
    echo "Replacing python header in (${f})"
    sed -i.bak --expression "1s@^.*@$env_python@" $f
done

popd 
# ------------------------------------------------------------------------------

