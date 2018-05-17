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

module purge
module load horton/2.1.0b3

if [ -n "$WEIGHTS_OUTFILE_CUBE" ] ; then
  EXTRA_OPTIONS="$EXTRA_OPTIONS --wsave '$WEIGHTS_OUTFILE_CUBE'"
fi

LNRHOREF=-18
until [  $LNRHOREF -gt -2 ]; do
    cmd="horton-esp-cost.py '${ESP_INFILE_CUBE}' '${COST_OUTFILE_HDF5}.${LNRHOREF}' --pbc 000 --overwrite --wdens ${DENS_INFILE_CUBE}:${LNRHOREF}:0.8 ${EXTRA_OPTIONS}"
    echo "Exectuting '$cmd'..."
    let LNRHOREF+=1
    eval "$cmd"
done
#cmd="horton-esp-cost.py '$ESP_INFILE_CUBE' '$COST_OUTFILE_HDF5' --pbc 000 --overwrite --wdens :${LNRHOREF}:0.8"
# echo "Exectuting '$cmd'..."

# horton-esp-cost.py "$ESP_INFILE_CUBE" "$COST_OUTFILE_HDF5" --pbc 000 --overwrite $EXTRA_OPTIONS
