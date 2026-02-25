#!/usr/bin/env bash

RM="rm -v -r -f"
# Remove all checkpoints exctpt the two most recents
find workspace/wfmanager.chk.bak.* -type f -printf '%T@\t%p\n' | sort -t $'\t' -g |  head -n -2 | cut -d $'\t' -f 2- | xargs -r $RM

