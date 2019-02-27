#!/bin/bash

# just a wrapper to load python 2 module necessary for Horton and loop over 
# a set of cost functions

args=$(getopt -n "$0" -l "infile-prefix:,outfile-prefix:,total-charge:" -o "i:o:q:" -- "$@")

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

eval set -- "$args"

# adapt file name parameters here!
INFILE_PREFIX="parameterized.cost.file.prefix"
OUTFILE_PREFIX="parameterized.outfile.prefix"

PDB_INFILE="system100.pdb"
TOP_INFILE="system100.top"
ATOMS_IN_CHARGE_GROUPS="atoms_in_charge_group.csv"
CHARGE_GROUPS_TOTAL_CHARGE="charge_group_total_charge.csv"
ATOMS_OF_SAME_CHARGE="atoms_of_same_charge.csv"
# FITTED_POINT_CHARGES_TXT="fitted_point_charges.txt"
# FITTED_POINT_CHARGES_TOP="fitted_point_charges.top"
QTOT=0

EXTRA_OPTIONS=""
echo "Got $# arguments: " $@
while true; do
  case "$1" in
    -i | --infile-prefix ) INFILE_PREFIX="$2"; shift; shift ;;
    -o | --outfile-prefix )  OUTFILE_PREFIX="$2"; shift; shift ;;
    -q | --total-charge )  QTOT="$2"; shift; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

source "$HOME/.bash_profile"

module purge
module load gpaw
# expects a parameter to be in the filename, i.e.
# "cost.file.lnrhoref.-9.h5"

FILES=${INFILE_PREFIX}.*.h5
echo "Looping over all ${FILES}..."
for f in $FILES; do
    suffix=${f#${INFILE_PREFIX}.}
    parameter=${suffix%.h5}
    logfile="${OUTFILE_PREFIX}.${parameter}.log"
    fitted_point_charges_csv="${OUTFILE_PREFIX}.${parameter}.csv"
    fitted_point_charges_top="${OUTFILE_PREFIX}.${parameter}.top"

    # A COMMAND OF THE FOLLOWING PATTERN IS CONSTRUCTED
    # fitESPconstrained.py system100.cost.h5  system100.pdb \
    # system100.top atoms_in_charge_group.csv \
    # charge_group_total_charge.csv atoms_of_same_charge.csv \
    # fitted_point_charges.txt fitted_point_charges.top \
    # --qtot 0 --verbose 2>&1 | tee esp-fit-constrained.log
    cmd="./fitESPconstrained.py '$f' '${PDB_INFILE}' '${TOP_INFILE}' \
        '${ATOMS_IN_CHARGE_GROUPS}' '${CHARGE_GROUPS_TOTAL_CHARGE}' \
        '${ATOMS_OF_SAME_CHARGE}' '${fitted_point_charges_top}' \
        '${fitted_point_charges_csv}' --qtot ${QTOT} --verbose \
        2>&1 | tee ${logfile}" 
    echo "Exectuting '$cmd'..."
    eval "$cmd" &
done
wait
