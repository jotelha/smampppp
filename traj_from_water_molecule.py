#!/usr/bin/env python
""" Constructs quantum system for ASE-internal water molecule descript in and extracts ESP"""
from ase.io import read
from ase.io import write

from gpaw import GPAW
from gpaw import restart

from ase.build import molecule
from ase.optimize.bfgslinesearch import BFGSLineSearch #Quasi Newton
from ase.units import Bohr
from ase.units import Hartree
import numpy as np

import os.path
import sys

import argparse

parser = argparse.ArgumentParser(description='Constructs quantum system' 
        ' for ASE-internal water molecule desciption and extracts electrostatic'
        ' potential (ESP)')
#parser.add_argument('infile')
parser.add_argument('-c', '--charge',metavar='INTEGER_CHARGE',type=int,nargs='?', const=1, default=0)
parser.add_argument('outfile_cube', nargs='?', metavar='outfile.cube',
        default='esp.cube', help="Electrostatic potential in GAUSSIAN-native"
        " .cube format, default 'esp.cube'")
parser.add_argument('outfile_csv', nargs='?', metavar='outfile.csv',
        default='esp.csv', help="Electrostatic potential and x,y,z coordinates"
                " as four-valued lines of .8 digits precision mantissa"
                " notation, default 'esp.csv'")
parser.add_argument('outfile_rho_cube', nargs='?', metavar='outfile_rho.cube',
        default='rho.cube', help="All-electron density in GAUSSIAN-native .cube"
        " format, default 'rho.cube'")
parser.add_argument('outfile_rho_pseudo_cube', nargs='?',
        metavar='outfile_rho_pseudo.cube', default='rho_pseudo.cube',
        help="All-electron density in GAUSSIAN-native .cube format, default"
        "'rho_pseudo.cube'")

args = parser.parse_args()

charge = args.charge

struc = molecule('H2O')
struc.set_pbc([0,0,0])
struc.set_cell([10,10,10])
struc.center()

calc  = GPAW(xc='PBE', h=0.2, charge=charge,
             spinpol=True, convergence={'energy': 0.001})
#struc.set_cell([25,25,25])
#struc.center()
struc.set_calculator(calc)

dyn = BFGSLineSearch(struc, trajectory='h2o.traj',
                     restart='bfgs_ls.pckl', logfile='BFGSLinSearch.log')
dyn.run(fmax=0.05)
