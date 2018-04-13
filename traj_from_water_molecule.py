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
        ' for ASE-internal water molecule desciption within non-periodic'
        ' 5x5x5 Angstrom box and optimizes structure')
#parser.add_argument('infile')
# parser.add_argument('-c', '--charge',metavar='INTEGER_CHARGE',type=int,nargs='?', const=1, default=0)
parser.add_argument('outfile_traj', nargs='?', metavar='outfile.traj',
        default='h2o.traj', help="DFT-optimized structure, default 'h2o.traj'")
parser.add_argument('outfile_pckl', nargs='?', metavar='bfgs_ls.pckl',
        default='bfgs_ls.pckl', help="Restart file, default 'bfgs_ls.pckl'")
parser.add_argument('logfile', nargs='?', metavar='BFGSLinSearch.log',
        default='BFGSLinSearch.log', help="Optimization's log, default 'BFGSLinSearch.log'")

args = parser.parse_args()

struc = molecule('H2O')
struc.set_pbc([0,0,0])
struc.set_cell([5,5,5])
struc.center()

calc  = GPAW(xc='PBE', h=0.2, charge=0,
             spinpol=True, convergence={'energy': 0.001})
struc.set_calculator(calc)

dyn = BFGSLineSearch(struc, trajectory=args.outfile_traj,
                     restart=args.outfile_pckl, logfile=args.logfile)
dyn.run(fmax=0.05)
