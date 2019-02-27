#!/bin/bash

# just a wrapper to load python 2 module necessary for Horton

args=$(getopt -n "$0" -l "sign,esp-infile-cube:,cost-outfile-hdf5:,dens-infile-cube:,weights-outfile-cube:" -o "se:o:d:w:" -- "$@")

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

eval set -- "$args"

ESP_INFILE_CUBE="esp.cube"
COST_OUTFILE_HDF5="cost.h5"

EXTRA_OPTIONS=""
echo "Got $# arguments: " $@
while true; do
  case "$1" in
    -e | --esp-infile-cube )   ESP_INFILE_CUBE="$2"; shift; shift ;;
    -o | --cost-outfile-hdf5 ) COST_OUTFILE_HDF5="$2"; shift; shift ;;
    -d | --dens-infile-cube )  DENS_INFILE_CUBE="$2"; shift; shift ;;
    -w | --weights-outfile-cube ) WEIGHTS_OUTFILE_CUBE="$2"; shift; shift;;
    -s | --sign )              EXTRA_OPTIONS="$EXTRA_OPTIONS --sign"; shift;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

source "$HOME/.bash_profile"

module purge
module load horton/2.1.0b3

if [ -n "$DENS_INFILE_CUBE" ] ; then
  EXTRA_OPTIONS="$EXTRA_OPTIONS --wdens '$DENS_INFILE_CUBE'"
fi
if [ -n "$WEIGHTS_OUTFILE_CUBE" ] ; then
  EXTRA_OPTIONS="$EXTRA_OPTIONS --wsave '$WEIGHTS_OUTFILE_CUBE'"
fi

cmd="horton-esp-cost.py '$ESP_INFILE_CUBE' '$COST_OUTFILE_HDF5' --pbc 000 --overwrite $EXTRA_OPTIONS"
echo "Exectuting '$cmd'..."

eval "$cmd"
# horton-esp-cost.py "$ESP_INFILE_CUBE" "$COST_OUTFILE_HDF5" --pbc 000 --overwrite $EXTRA_OPTIONS
