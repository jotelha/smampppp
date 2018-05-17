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

module purge
module load horton/2.1.0b3

FILES=${COST_INFILE_PREFIX}.*
echo "Looping over all ${FILES}..."
for f in $FILES; do
    suffix=${f#${COST_INFILE_PREFIX}.}
    logfile="${CHARGE_OUTFILE_PREFIX}.${suffix}"
    cmd="./fitESPforBenzene.py '$f' 2>&1 | tee ${logfile}"
    echo "Exectuting '$cmd'..."
    eval "$cmd"
done
