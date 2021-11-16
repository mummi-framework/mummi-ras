#!/usr/bin/env bash
# ------------------------------------------------------------------------------

file=$1

exit_counter=0
while [ ! -f $file ]
do
  echo "($file) does not exist yet."
  sleep 10s

  exit_counter=$((exit_counter + 1))
  if [ "$exit_counter" -eq "100" ]; then
    echo "Timeout: Failed to find file ($file)."
    exit 1
  fi
done
echo "Found file ($file)."

# ------------------------------------------------------------------------------
