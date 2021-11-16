#!/usr/bin/env bash
# ------------------------------------------------------------------------------

dir=$1
postfix=0
backup_dir=$(dirname $(readlink -f "${dir}"))/backup
backup_file=${backup_dir}/$(basename ${dir})

if [[ -d "${dir}" ]] ; then
    while [[ -d "${backup_file}_$postfix" ]]
    do
      postfix=$(( $postfix + 1 ))
    done

    dir2=${backup_file}_$postfix

    if [[ ! -d "${backup_dir}" ]] ; then
        mkdir "${backup_dir}"
    fi

    cp -rf $dir $dir2
    echo "  > Moving ("$dir") to ("$dir2")"
fi

# ------------------------------------------------------------------------------
