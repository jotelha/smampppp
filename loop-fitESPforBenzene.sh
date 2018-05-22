#!/bin/bash

# just a wrapper to load python 2 module necessary for Horton

args=$(getopt -n "$0" -l "cost-infile-prefix:,charge-outfile-prefix:" -o "i:o:" -- "$@")

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

eval set -- "$args"

COST_INFILE_PREFIX="cost.h5"

EXTRA_OPTIONS=""
echo "Got $# arguments: " $@
while true; do
  case "$1" in
    -i | --cost-infile-prefix ) COST_INFILE_PREFIX="$2"; shift; shift ;;
    -o | --charge-outfile-prefix )  CHARGE_OUTFILE_PREFIX="$2"; shift; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

source "$HOME/.bash_profile"

FILES=${COST_INFILE_PREFIX}.*.h5
echo "Looping over all ${FILES}..."
for f in $FILES; do
    suffix=${f#${COST_INFILE_PREFIX}.}
    parameter=${suffix%.h5}
    logfile="${CHARGE_OUTFILE_PREFIX}.${parameter}.log"
    cmd="./fitESPforBenzene.py '$f' 2>&1 | tee ${logfile}"
    echo "Exectuting '$cmd'..."
    eval "$cmd"
done
