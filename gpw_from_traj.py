#!/usr/bin/env python
""" Reconstructs quantum system from GPAW .traj file given as command line argument and stores it in .gpw format. """

from ase.io import read
from ase.io import write
from gpaw import Mixer
from gpaw import GPAW, FermiDirac

from gpaw import GPAW
from gpaw import restart

from gpaw import GPAW
from ase.optimize.bfgslinesearch import BFGSLineSearch #Quasi Newton
from ase.units import Bohr

import os.path
import argparse

parser = argparse.ArgumentParser(description='Reconstructs quantum system from'
    ' GPAW .traj file given as command line argument and stores it in'
    ' .gpw format.')
parser.add_argument('-c', '--charge',metavar='CHARGE', type=float, 
    nargs='?', const=1.0, default=0.0, help="The system's total charge,"
    " 0.0 if not specified, 1.0 if specified without value.")
parser.add_argument('-b', '--box', metavar=('X','Y','Z'), type=float,
    nargs=3, default=[30.0,30.0,30.0], help="The bounding box' measures in"
    " Angstrom, default [30.0, 30.0, 30.0]")
parser.add_argument('infile_traj')
parser.add_argument('outfile_gpw')
args = parser.parse_args()

traj_file = args.infile_traj
gpw_file  = args.outfile_gpw
charge    = args.charge
box       = args.box
print("Reading structure of charge {} (e) in box {} (Angstrom) from input"
    " file '{}'".format(charge,box, traj_file))
struc = read(traj_file) 
# to read last frame, use "filename@-1" 
# or read(traj_file, index=-1)
# However, that's implicitly done by default, as stated on
# https://wiki.fysik.dtu.dk/ase/ase/io/io.html:
# "The last configuration will be returned by default"

calc  = GPAW(xc='PBE', h=0.2, charge=+1,
             spinpol=True, convergence={'energy': 0.001},
             mixer=Mixer(beta=0.25, nmaxold=10, weight=1.0), occupations=FermiDirac(width=0.1))
struc.set_cell(box)
struc.set_pbc([0,0,0])
struc.center()
struc.set_calculator(calc)

Epot  = struc.get_potential_energy()

if gpw_file == '':
    gpw_file = os.path.splitext(traj_file)[0] + '.gpw'

print("Writing output file '{}'...".format(gpw_file))
calc.write(gpw_file, mode='all') 
# The file .gpw is a binary file containing wave functions, densities, positions and everything else (also the 
# parameters characterizing the PAW calculator used for the calculation).
# source: https://wiki.fysik.dtu.dk/gpaw/documentation/manual.html#restarting-a-calculation
