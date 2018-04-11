#!/usr/bin/env python
# cube2gesp, N atoms, M grid points
from ase.io.cube import read_cube_data
from ase.units import Bohr
import numpy as np
import io

import argparse

parser = argparse.ArgumentParser(description='Convert .cube format to'
        ' Fortran-RESP-readible .gesp format.')
parser.add_argument('infile_cube', metavar='infile.cube',
        help="Input in .cube format (Bohr and Hartree)")
parser.add_argument('outfile_gesp', nargs='?', metavar='outfile.gesp',
       default='outfile.gesp', help="Output in .gesp format")
parser.add_argument('-n', nargs='?',
        metavar='N', default=-1, const=99999, type=int,
        help="Maximum number of grid points")
args = parser.parse_args()

# cube file measures in Bohr, but ASE converts input to Angstrom
# however, field is not converted while distances are
cube_data, cube_atoms = read_cube_data(args.infile_cube)

nX = cube_data.shape
X = cube_atoms.cell.diagonal() # measures of cell in spatial dimensions
x_grid = np.linspace(0,X[0],nX[0])
y_grid = np.linspace(0,X[1],nX[1])
z_grid = np.linspace(0,X[2],nX[2])

x_grid3,y_grid3,z_grid3=np.meshgrid(x_grid,y_grid,z_grid)

# atom positions in N x 3 array (rows x cols)
pos = cube_atoms.get_positions()
# potential on grid in M x 4 array (rows x cols)
dat = np.vstack( ( cube_data.flatten(), x_grid3.flatten()/Bohr, y_grid3.flatten()/Bohr, z_grid3.flatten()/Bohr ) ).T

n_atoms = len(cube_atoms)
n_dat = dat.shape[0]
l1 = "{:5d}{:6d}{:5d}".format(n_atoms,n_dat,0) # first line
ghd = io.BytesIO() # .gesp header buffer
np.savetxt(ghd, pos, fmt='%16.7E', delimiter='', header=l1, comments='')
ghd = ('\n' + ' '*16).join(ghd.getvalue().decode().splitlines())
np.savetxt(args.outfile_gesp, dat, fmt='%16.7E', delimiter='' ,header=ghd, comments='')

