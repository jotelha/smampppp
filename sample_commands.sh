#!/bin/bash -x
# sample work flow on NEMO for extracting ESP and electron density 
# from a GPAW  simulation

# run (parallel mpirun) gpaw and (non-parallel, but resource-intense)
# horton-esp-cost.py steps as jobs in interactive session started by

msub -I -l walltime=2:00:00 -l nodes=1:ppn=20 # start interactive 2h session

module purge
module load gpaw/1.3.0
# extract .gpw file from SMAMP molecule of charge +6
mpirun -n 20 gpaw-python gpw_from_traj.py -c 6 molecule.traj smamp.gpw

# write ESP and electron densities to standard files
mpirun -n 20 gpaw-python esp_from_gpw.py smamp.gpw
# usage: esp_from_gpw.py [-h]
#                        infile [outfile.cube] [outfile.csv] [outfile_rho.cube]
#                        [outfile_rho_pseudo.cube]
# 
# Extracts the electrostatic potential (ESP) from a GPAW .gpw restart file given
# as command line argument.
# 
# positional arguments:
#   infile
#   outfile.cube          Electrostatic potential in GAUSSIAN-native .cube
#                         format, default 'esp.cube'
#   outfile.csv           Electrostatic potential and x,y,z coordinates as four-
#                         valued lines of .8 digits precision mantissa notation,
#                         default 'esp.csv'
#   outfile_rho.cube      All-electron density in GAUSSIAN-native .cube format,
#                         default 'rho.cube'
#   outfile_rho_pseudo.cube
#                         All-electron density in GAUSSIAN-native .cube format,
#                         default'rho_pseudo.cube'
# 
# optional arguments:
#   -h, --help            show this help message and exit


# with all-electron density for weighting:
# construct united-atom cube files by truncating implicit atoms:
./aa2ua_cube.py system100.pdb system100.top esp.cube esp_ua.cube
./aa2ua_cube.py system100.pdb system100.top rho.cube rho_ua.cube
# construct cost function with weighting from electron density:

module purge
module load horton/2.1.0b
horton-esp-cost.py smamp_esp_hartree_ua.cube smamp_esp_cost_wdens_ua.h5 --pbc 000 --wdens rho_ua.cube --overwrite
horton-esp-fit.py -q 6 smamp_esp_cost_wdens_ua.h5 smamp_esp_charges_wdens_ua.h5 --overwrite
# horton-esp-fit.py possesses the possibility to impose symmetry constraints via
#   --ridge RIDGE         The thikonov regularization strength used when solving
#                         the charges. [default=0.0]
#   --symmetry SYMMETRY SYMMETRY
#                         Perform a symmetry analysis on the charges. This
#                         option requires two arguments argument in the
#                         following order: the cube file used to construct the
#                         ESP cost function and a CIF file with the generators
#                         of the symmetry of this system and a primitive unit
#                         cell.
# TODO: we have to figure out how to impose the symmetry constraints here

# if --overwrite not set, HORTON will not update possibly existin output files

# map fitted charges back onto the original united-atoms in GROMACS topology file

module purge
module load gpaw/1.3.0
./mapCharges.py system100.pdb system100.top smamp_esp_charges_wdens_ua.h5 smamp_esp_charges_wdens_ua.top
# usage: mapCharges.py [-h] [-i INSERTION_RULES] [-aa]
#                      [infile.pdb] [infile.top] [infile.h5] [outfile.top]
# 
# Maps point charges obtained by GPAW and HORTON on the original GROMACS
# topology initially modified by insertHbyList.py
# 
# positional arguments:
#   infile.pdb            Original .pdb file, before insertion of implicit
#                         hydrogen.
#   infile.top            Original GROMACS .top file
#   infile.h5             Point charges in hdf5 format by horton-esp-fit.py
#   outfile.top           GROMACS .top output filewith updated charges according
#                         to given .hdf5
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   -i INSERTION_RULES, --insertion-rules INSERTION_RULES
#                         A string representation of a python dictionary,
#                         describing how many implicit hydrogens have been
#                         inserted at which atom. Example:
#                         {'CD4':1,'CD3':1,'CA2':2,'CA3':2,'CB2':2,'CB3':2}
#   -aa, --all-atoms      Determines whether the charges are for all-atoms (or
#                         united-atoms, default) representation
# 


# ./aa2ua_cube.py system100.pdb system100.top rho_pseudo.cube rho_pseudo_ua.cub

# similar procedure for an all-atom fit, just using the unmodified cube files:
module purge
module load horton/2.1.0b3 
horton-esp-cost.py smamp_esp_hartree.cube smamp_esp_cost_wdens.h5 --pbc 000 --wdens rho.cube --overwrite
horton-esp-fit.py -q 6 smamp_esp_cost_wdens.h5 smamp_esp_charges_wdens.h5 --overwrite

module purge
module load gpaw/1.3.0
./mapCharges.py -aa system100.pdb system100.top smamp_esp_charges_wdens.h5 smamp_esp_charges_wdens.top

# attention: since horton is based on python 2.7 while all other software 
# depends on python 3.6, the module purgin and reloading is necessary
# if steps are carried out sequentially within one session as shown here
