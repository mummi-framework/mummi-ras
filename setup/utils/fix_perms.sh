#!/usr/bin/env bash
# ------------------------------------------------------------------------------

group=$1
mdepth=30

umask 007
echo '> changing permissions on all files owned by ('$USER') to be accessible by group ('$group')'
cd $MUMMI_ROOT
for entry in `ls .`; do
  echo '  > fixing perms in' $entry
  if [[ -d $entry ]]; then
    find $entry -maxdepth $mdepth -user $USER ! -group $group | xargs -r -d '\n' chgrp $group
    find $entry -maxdepth $mdepth -user $USER | xargs -r -d '\n' chmod g=u
  else
    chgrp $group $entry
    chmod g=u $entry
  fi
done

# ------------------------------------------------------------------------------
