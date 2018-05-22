#!/bin/bash
#set -e
#set -o pipefail

# sample call:
# ./loop-esp-cost.sh --esp-infile-cube benzene/benzene.esp.cube \
#  --dens-infile-cube benzene/benzene.rho.cube \
#  --cost-outfile-hdf5 benzene/wdens/benzene.cost.lnrhoref \
#  --weights-outfile-cube benzene/wdens/benzene.weights.lnrhoref \
#  --sign 2>&1 | tee benzene/wdens/benzene.loop-esp-cost.log

LONG_ARGS="sign,esp-infile-cube:,cost-outfile-hdf5:,dens-infile-cube:"
LONG_ARGS="${LONG_ARGS},weights-outfile-cube:,lnrhoref-min:,lnrhoref-max:"
args=$(getopt -n "$0" -l "${LONG_ARGS}" -o "se:o:d:w:" -- "$@")

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

eval set -- "$args"

ESP_INFILE_CUBE="esp.cube"
COST_OUTFILE_HDF5="cost.h5"
DENS_INFILE_CUBE=
LNRHOREF_MIN=-18
LNRHOREF_MAX=3
EXTRA_OPTIONS=""
echo "Got $# arguments: " $@
while true; do
  case "$1" in
    -e | --esp-infile-cube )   ESP_INFILE_CUBE="$2"; shift; shift ;;
    -o | --cost-outfile-hdf5 ) COST_OUTFILE_HDF5="$2"; shift; shift ;;
    -d | --dens-infile-cube )  DENS_INFILE_CUBE="$2"; shift; shift ;;
    -w | --weights-outfile-cube ) WEIGHTS_OUTFILE_CUBE="$2"; shift; shift;;
    --lnrhoref-min )           LNRHOREF_MIN="$2"; shift; shift;;
    --lnrhoref-max )           LNRHOREF_MAX="$2"; shift; shift;;
    -s | --sign )              EXTRA_OPTIONS="$EXTRA_OPTIONS --sign"; shift;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

source "$HOME/.bash_profile"

module purge
module load horton/2.1.0b3

LNRHOREF=$LNRHOREF_MIN
echo "Looping whole integer values from ln(rho_ref) = ${LNRHOREF_MIN}" \
     " to ${LNRHOREF_MAX}... "
until [  $LNRHOREF -gt $LNRHOREF_MAX ]; do
    if [ -n "$WEIGHTS_OUTFILE_CUBE" ] ; then
      CUR_EXTRA_OPTIONS="$EXTRA_OPTIONS --wsave '${WEIGHTS_OUTFILE_CUBE}.${LNRHOREF}.cube'"
      # .cube extension is necessary in the end for Horton to recognize file format
    fi

    cmd="horton-esp-cost.py '${ESP_INFILE_CUBE}' '${COST_OUTFILE_HDF5}.${LNRHOREF}.h5'"
    cmd="${cmd} --pbc 000 --overwrite --wdens '${DENS_INFILE_CUBE}:${LNRHOREF}:0.8'"
    cmd="${cmd} ${CUR_EXTRA_OPTIONS}"
    echo "Exectuting '$cmd'..."
    let LNRHOREF+=1
    eval "$cmd" &
done
wait
