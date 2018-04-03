#!/usr/bin/env python
""" Reconstructs quantum system from GPAW .traj file given as command line argument and stores it in .gpw format. """

from ase.io import read
from ase.io import write

from gpaw import GPAW
from gpaw import restart

from gpaw import GPAW
from ase.optimize.bfgslinesearch import BFGSLineSearch #Quasi Newton
from ase.units import Bohr

import os.path
import argparse

parser = argparse.ArgumentParser(description='Reconstructs quantum system from GPAW .traj file given as command line argument and stores it in .gpw format.')
parser.add_argument('-c', '--charge',metavar='INTEGER_CHARGE',type=int,nargs='?', const=1, default=0)
parser.add_argument('infile')
parser.add_argument('outfile', nargs='?',default='')
args = parser.parse_args()

traj_file = args.infile
gpw_file = args.outfile
charge = args.charge
print("Reading structure of charge {} from input file '{}'".format(charge, traj_file))
struc = read(traj_file)

calc  = GPAW(xc='PBE', h=0.2, charge=charge,
             spinpol=True, convergence={'energy': 0.001})
struc.set_cell([25,25,25])
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
