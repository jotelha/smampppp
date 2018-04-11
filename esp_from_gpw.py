#!/usr/bin/env python
""" Extracts the electrostatic potential (ESP) from a GPAW .gpw restart file given as command line argument. """
#setup the gpaw calculation

from ase.io import read
from ase.io import write

from gpaw import GPAW
from gpaw import restart

from ase.units import Bohr
from ase.units import Hartree
import numpy as np

import argparse
import io

# from mailing list
# https://listserv.fysik.dtu.dk/pipermail/gpaw-users/2012-February/001277.html
# In PAW, the effective potential is vHartree + vxc + vbar, where vbar are 
# fixed contributions from the atomic setups.

# https://listserv.fysik.dtu.dk/pipermail/gpaw-users/2012-February/001283.html
# For that case [periodic], the average of the Coulomb potential is zero, where the 
 # Coulomb potential is the potential we get from solving the Poisson 
# equation "nabla v_coulomb=4 pi rho", where rho is the pseudo charge 
# density.  The total effective potential is the v_coulomb+v_XC.

# vt_sg is KS potential on fine grid (?), sG on coarse grid
def get_vtsg_potential(calc): # no conversion!
        dens = calc.density
        
        # the following lines have been adapted from GPAW
        # code without fully understanding their role
        calc.initialize_positions()
        dens.interpolate_pseudo_density()
        dens.calculate_pseudo_charge()
        calc.hamiltonian.update(dens)

        v_g = calc.hamiltonian.finegd.collect(calc.hamiltonian.vt_sg, broadcast=True)
        v_g = calc.hamiltonian.finegd.zero_pad(v_g)
        v_g = v_g[0] + v_g[1] # sum for spin up & down
        return v_g

parser = argparse.ArgumentParser(description='Extracts the electrostatic'
        ' potential (ESP) from a GPAW .gpw restart file given as command line'
        ' argument.')
parser.add_argument('infile')
parser.add_argument('outfile_vHtg_cube', nargs='?', metavar='vHtg.cube', 
       default='vHtg.cube', help="Electrostatic potential (vHt_g) in GAUSSIAN-native"
       " .cube format, default 'vHtg.cube'")
parser.add_argument('outfile_vtsg_cube', nargs='?', metavar='vtsg.cube',
        default='vtsg.cube', help="Electrostatic potential (vt_sg) in GAUSSIAN-native"
        " .cube format, default 'vtsg.cube'")
parser.add_argument('outfile_vHtg_gesp', nargs='?', metavar='vHtg.gesp', 
       default='vHtg.gesp', help="Electrostatic potential and x,y,z coordinates"
               " as four-valued lines of 16.7E formatted numbers"
               " default 'vHtg.gesp'")
parser.add_argument('outfile_rho_cube', nargs='?', metavar='rho.cube',
        default='rho.cube', help="All-electron density in GAUSSIAN-native .cube"
        " format, default 'rho.cube'")
parser.add_argument('outfile_rho_pseudo_cube', nargs='?', 
        metavar='rho_pseudo.cube', default='rho_pseudo.cube',
        help="Pseudo electron density in GAUSSIAN-native .cube format, default"
        "'rho_pseudo.cube'")

args = parser.parse_args()

gpw_file = args.infile
print("Reading input file '{}'".format(gpw_file))

# The file .gpw is a binary file containing wave functions, densities, 
# positions and everything else (also the parameters characterizing 
# the PAW calculator used for the calculation).
# source: https://wiki.fysik.dtu.dk/gpaw/documentation/manual.html#restarting-a-calculation
struc, calc = restart(gpw_file)

vHtg = calc.get_electrostatic_potential()
# potential query comes from gpaw/paw.py
##       def get_electrostatic_potential(self):
#        """Return the electrostatic potential.
#
#        This is the potential from the pseudo electron density and the
#        PAW-compensation charges.  So, the electrostatic potential will
#        only be correct outside the PAW augmentation spheres.
#        """
#
#        ham = self.hamiltonian
#        dens = self.density
#        self.initialize_positions()
#        dens.interpolate_pseudo_density()
#        dens.calculate_pseudo_charge()
#        return ham.get_electrostatic_potential(dens) * Ha
# ...which in turn calls get_electrostatic_potential(dens) from hamiltonian.py...
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
nX = vHtg.shape # = 2 * calc.get_number_of_grid_points()
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
vHtg_hartree = vHtg / Hartree                   # some (what?) potential in Hartree
vtsg_hartree = get_vtsg_potential(calc)   # Hartree potential in Hartree
# 20180410, by Richard Leute:
#   Richard:
#      Zu meiner Frage, das Elektrostatische Feld (E) hätte ich aus dem Hartree
#      Potential (H) abgeleitet (E = -grad(H)). Jetzt sollte ich mir nur sicher
#      sein was das Hartree Potential in GPAW ist. Bekomme ich das mit, 
#      'calc.hamiltonian.vt_sG'? Das habe ich so aus der mailing list genommen 
#      https://listserv.fysik.dtu.dk/pipermail/gpaw-users/2012-February/001275.html
#      Die Einheiten sind dann Hartree und Bohr?
#   Michael:
#      ja, sollte genauso sein.


# put potential in grid points and xyz-coordinates in csv-file format (four %.8e values, seperated by whitespace)
#as expected by resp FORTRAN code 2.1 (October 1994 Jim Caldwell)
pos = struc.get_positions()
dat = np.vstack( ( vHtg_hartree.flatten(), x_grid3.flatten()/Bohr, y_grid3.flatten()/Bohr, z_grid3.flatten()/Bohr ) ).T
# spatial units are converted to Bohr. What unit is the potential?
# Division (not multiplication)  is necessary here, as ase.units.Bohr is defined as
#     u['Bohr'] = (4e10 * pi * u['_eps0'] * u['_hbar']**2 / u['_me'] / u['_e']**2)  # Bohr radius
# with unit [ Bohr ] = Angstrom / Bohr
#     ase.units.Bohr = 0.5291772105638411 (Ang)

write(args.outfile_vHtg_cube, struc, data=vHtg_hartree) # apparently the native GAUSSIAN format for ESP, readible by Horton
write(args.outfile_vtsg_cube, struc, data=vtsg_hartree)


n_atoms = len(struc)
n_dat = dat.shape[0]
l1 = "{:5d}{:5d}{:5d}".format(n_atoms,n_dat,0)
ghd = io.BytesIO()
np.savetxt(ghd, pos, fmt='%16.7E', delimiter='', header=l1, comments='')
ghd = ('\n' + ' '*16).join(ghd.getvalue().decode().splitlines())
np.savetxt(args.outfile_vHtg_gesp, dat, fmt='%16.7E', delimiter='' ,header=ghd, comments='')
# horton/io/cube.py, line 65:
#     all coordinates in a cube file are in atomic units
# in line with
# http://theochem.github.io/horton/2.1.0/lib/mod_horton_units.html?highlight=units#module-horton.units
# referencing
#    B. J. Mohr and B. N. Taylor, CODATA recommended values of the fundamental 
#    physical constants: 1998, Rev. Mod. Phys. 72(2), 351 (2000)
# which recommends 0.529e-10 m as the length measure in a.u. (atomic units)
# thus the only conclusion can be that Horton expects Bohr as length unit
# and Hartree as energy unit. Conversions for outputing in this file are correct.


# the code below could be used to interpolate the potential onto a coarser or finer grid
# however, I got errors (most of the time)

# Transformer:
# t = PS2AE(calc, h=0.4)

# phi_grid_ps = t.get_electrostatic_potential(ae=False) # pseude wave-functions
# phi_grid_ae = t.get_electrostatic_potential() # all-electron wave-functions

# x_grid3,y_grid3,z_grid3=np.meshgrid(t.gd.coords(0)*Bohr,t.gd.coords(1)*Bohr,t.gd.coords(2)*Bohr)

# dat_ae = np.array( [ phi_grid_ae.flatten(), x_grid3.flatten(), y_grid3.flatten(), z_grid3.flatten() ] )
# np.savetxt('phi_grid_ae.csv',dat_ae.T,fmt='%.8e',delimiter=' ')

# dat_ps = np.array( [ phi_grid_ps.flatten(), x_grid3.flatten(), y_grid3.flatten(), z_grid3.flatten() ] )
# np.savetxt('phi_grid_ps.csv',dat_ps.T,fmt='%.8e',delimiter=' ')



# https://wiki.fysik.dtu.dk/gpaw/tutorials/bader/bader.html#bader-analysis
rho_pseudo      = calc.get_pseudo_density()
rho             = calc.get_all_electron_density()
# https://wiki.fysik.dtu.dk/gpaw/tutorials/all-electron/all_electron_density.html:
# As the all-electron density has more structure than the pseudo-density, it is
# necessary to refine the density grid used to represent the pseudo-density.
# This can be done using the gridrefinement keyword of the
# get_all_electron_density method:
#
# >>> n = calc.get_all_electron_density(gridrefinement=2)
#
# Current only the values 1, 2, and 4 are supported (2 is default).

# https://wiki.fysik.dtu.dk/gpaw/tutorials/bader/bader.html
# gives an example on how to convert and extract the electron densities:
rho_pseudo_per_bohr_cube = rho_pseudo * Bohr**3
rho_per_bohr_cube = rho * Bohr**3
write(args.outfile_rho_cube, struc, data=rho_per_bohr_cube) 
write(args.outfile_rho_pseudo_cube, struc, data=rho_pseudo_per_bohr_cube) 

# looking at the horton parsing script 
#    horton/scripts/espfit.py
# one finds...

# def parse_wdens(arg):
#     '''Parse the argument to the --wdens option of horton-espfit.py'''
#     if arg is None:
#         return
#     words = arg.split(':')
#     lnrho0 = -9
#     sigma = 0.8
#     if len(words) == 0:
#         fn_cube = None
#     elif len(words) == 1:
#         fn_cube = words[0]
#     elif len(words) == 2:
#         fn_cube = words[0]
#         lnrho0 = float(words[1])
#     elif len(words) == 3:
#         fn_cube = words[0]
#         lnrho0 = float(words[1])
#         sigma = float(words[2])
#     else:
#         raise ValueError('The argument to --wdens may at most contain three fields separated by a colon.')
#     if len(fn_cube) == 0:
#         fn_cube = None
#     return fn_cube, lnrho0, sigma

# ...

# def load_rho(coordinates, numbers, fn_cube, ref_ugrid, stride, chop):
#     '''Load densities from a file, reduce by stride, chop and check ugrid
# 
#        **Arguments:**
# 
#        coordinates
#             An array with shape (N, 3) containing atomic coordinates.
# 
#        numbers
#             A vector with shape (N,) containing atomic numbers.
# 
#        fn_cube
#             The cube file with the electron density.
# 
#        ref_ugrid
#             A reference ugrid that must match the one from the density cube
#             file (after reduction).
# 
#        stride
#             The reduction factor.
# 
#        chop
#             The number of slices to chop of the grid in each direction.
#     '''
#     if fn_cube is None:
#         # Load the built-in database of proatoms
#         natom = len(numbers)
#         numbers = np.unique(numbers)
#         proatomdb = ProAtomDB.from_refatoms(numbers, max_cation=0, max_anion=0, agspec='fine')
#         # Construct the pro-density
#         rho = np.zeros(ref_ugrid.shape)
#         for i in xrange(natom):
#             spline = proatomdb.get_spline(numbers[i])
#             ref_ugrid.eval_spline(spline, coordinates[i], rho)
#     else:
#         # Load cube
#         mol_rho = IOData.from_file(fn_cube)
#         rho = mol_rho.cube_data
#         ugrid = mol_rho.grid
#         # Reduce grid size
#         if stride > 1:
#             rho, ugrid = reduce_data(rho, ugrid, stride, chop)
#         # Compare with ref_ugrid (only shape)
#         if (ugrid.shape != ref_ugrid.shape).any():
#             raise ValueError('The densities file does not contain the same amount if information as the potential file$
#     return rho

# ... not performing any conversions on the electron density read from a .cube file,
# thus the density unit should be expected as Bohr^-3
