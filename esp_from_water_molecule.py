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

parser = argparse.ArgumentParser( \
        description='Constructs quantum system for ASE-internal water molecule description and extracts ESP')
parser.add_argument('-c', '--charge',metavar='INTEGER_CHARGE',type=int,nargs='?', const=1, default=0)
parser.add_argument('outfile_cube', nargs='?', metavar='outfile.cube', default='phi_grid.cube', \
        help="Electrostatic potential in GAUSSIAN-native .cube format, default 'phi_grid.cube'")
parser.add_argument('outfile_csv', nargs='?', metavar='outfile.csv', default='phi_grid.csv', \
        help="Electrostatic potential and x,y,z coordinates as four-valued lines of .8 digits precision " \
                + " mantissa notation, default 'phi_grid.csv'")
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

Epot  = struc.get_potential_energy()

# https://wiki.fysik.dtu.dk/gpaw/devel/electrostatic_potential.html tells us, the
# get_electrostatic_corrections() method will return an array of integrated corrections with the unit
#     [corrections] = eV Angstrom^3
# However, 
# https://wiki.fysik.dtu.dk/gpaw/tutorials/ps2ae/ps2ae.html?highlight=get_electrostatic_potential#gpaw.utilities.ps2ae.PS2AE.get_electrostatic_potential
# states, the interpolated ESP PS2AE.get_electrostatic_potential(ae=True, rcgauss=0.02) is given in
#     [U_interpolated] = eV
# No unit information has been found on the gpaw.calculator.GPAW.get_electrostatic_potential() called here,
# however we assume homogeneous units throughout GPAW for now
phi = calc.get_electrostatic_potential()
# potential query comes from gpaw/hamiltonian.py
#     def get_electrostatic_potential(self, dens):
#        self.update(dens)
#
#        v_g = self.finegd.collect(self.vHt_g, broadcast=True)
#        v_g = self.finegd.zero_pad(v_g)
#        if hasattr(self.poisson, 'correction'):
#            assert self.poisson.c == 2
#            v_g[:, :, 0] = self.poisson.correction
#        return v_g
#
# A comment from the same file ...
#    The XC-potential and the Hartree potential are evaluated on the fine grid, and the sum is then restricted 
#    to the coarse grid.
# ... and a note from https://wiki.fysik.dtu.dk/gpaw/algorithms.html?highlight=fine%20grid ...
#    Finite-difference (FD):
#    Uniform real-space orthorhombic grids. Two kinds of grids are involved in the calculations: 
#    A coarse grid used for the wave functions and a fine grid (23=8 times higher grid point density) used for 
#    densities and potentials. The pseudo electron density is first calculated on the coarse grid from the wave 
#    functions, and then interpolated to the fine grid, where compensation charges are added for achieving 
#    normalization. The effective potential is evaluated on the fine grid (solve the Poisson equation and calculate 
#    the exchange-correlation potential) and then restricted to the coarse grid where it needs to act on the wave 
#    functions (also on the coarse grid).
# ... tell us: potential has twice as many grid points in each spatial dimension as the actual number of coarse grid
# points queried by "calc.get_number_of_grid_points()"
nX = phi.shape # = 2 * calc.get_number_of_grid_points()
X = struc.cell.diagonal()
x_grid = np.linspace(0,X[0],nX[0])
y_grid = np.linspace(0,X[1],nX[1])
z_grid = np.linspace(0,X[2],nX[2])

x_grid3,y_grid3,z_grid3=np.meshgrid(x_grid,y_grid,z_grid)

# https://theochem.github.io/horton/2.1.0b3/lib/mod_horton_units.html?highlight=units#module-horton.units
# apparently, Horton internally uses atomic units.
# If this applies strictly, we have electron mass m_e = 1, electron charge e = 1,
# reduced Planck's constant h_bar = 1 and Coulomb force constant k_e = 1 / (4 Pi eps_0 ) = 1 per definition
# Furthermore, it should expect
# length in Bohr (a_0) , defined as 4 Pi eps_0 h_bar^2 / (m_e e^2) = 1
# energy in Hartree (E_h), defined as m_e e^4 / (4 Pi eps_0 h_bar)^2 = 1
# electric potential, defined as E_h / e = 1

# thus, GPAW potential in units of "eV" 
# are to be converted to units of "E_h / e = m_e e^3 / (4 Pi eps_0 h_bar)^2"
#     U_hor = U_gpw * E_h / (e*eV)
# we use 
#    ase.units.Hartree = 27.211386024367243 (eV)
phi_hartree = phi / Hartree

# put potential in grid points and xyz-coordinates in csv-file format (four %.8e values, seperated by whitespace)
#as expected by resp FORTRAN code 2.1 (October 1994 Jim Caldwell)
dat = np.vstack( ( phi_hartree.flatten(), x_grid3.flatten()/Bohr, y_grid3.flatten()/Bohr, z_grid3.flatten()/Bohr ) ).T
# spatial units are converted to Bohr. What unit is the potential?
# Division (not multiplication)  is necessary here, as ase.units.Bohr is defined as
#     u['Bohr'] = (4e10 * pi * u['_eps0'] * u['_hbar']**2 / u['_me'] / u['_e']**2)  # Bohr radius
# with unit [ Bohr ] = Angstrom / Bohr
#     ase.units.Bohr = 0.5291772105638411 (Ang)

write(args.outfile_cube, struc, data=phi_hartree) # apparently the native GAUSSIAN format for ESP, readible by Horton
np.savetxt(args.outfile_csv,dat,fmt='%.8e',delimiter=' ')

#grid = calc.hamiltonian.gd.get_grid_point_coordinates()

#dat2 = np.vstack( ( phi.flatten(), np.concatenate(([0],grid[0,:,:,:].flatten().T),axis=0).T, \
#        np.concatenate(([0],grid[1,:,:,:].flatten().T),axis=0).T, \
#        np.concatenate(([0],grid[2,:,:,:].flatten().T),axis=0).T ) ).T

#np.savetxt(args.outfile_csv +'.compare',dat2,fmt='%.8e',delimiter=' ')
