#!/usr/bin/env python
# cube2gesp, N atoms, M grid points
from ase.io.cube import read_cube_data
from ase.units import Bohr
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import LinearNDInterpolator

import numpy as np
import io

import argparse

def main():
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
    
    cube2gesp(args.infile_cube, args.outfile_gesp, args.n)

def cube2gesp(infile_cube, outfile_gesp, N = -1):
    # cube file measures in Bohr, but ASE converts input to Angstrom
    # however, field is not converted while distances are
    cube_data, cube_atoms = read_cube_data(infile_cube)

    # nX = cube_data.shape
    # X = cube_atoms.cell.diagonal() # measures of cell in spatial dimensions
    # x_grid = np.linspace(0,X[0],nX[0])
    # y_grid = np.linspace(0,X[1],nX[1])
    # z_grid = np.linspace(0,X[2],nX[2])
    
    # x_grid3,y_grid3,z_grid3=np.meshgrid(x_grid,y_grid,z_grid)
    
    # general approach for creating a grid of coordinates
    # in 3 dimensions and N points in each spatial dimension,
    # X.shape == (3, N, N, N)
    dim = cube_data.ndim # usually 3
    print("Read .cube file of shape {} in box {}.".format(cube_data.shape,cube_atoms.cell.diagonal()))
    X_lin = []
    X = np.empty((dim,) + cube_data.shape)
    for i in range(0,dim):
        X_lin.append( np.linspace(0, cube_atoms.cell[i,i], cube_data.shape[i]) )
    X_list=np.meshgrid(*X_lin,indexing='ij')
    X = np.asarray(X_list)
    Unit_Cell = X[:,1,1,1].prod() # for uniform grid
    Integral = cube_data.flatten().sum()*Unit_Cell
    print("Reconstructed uniform grid of shape {}.".format(X.shape))
    
    if N > 0:
        print("Interpolate .cube data of shape {} and {} points in total onto grid"
            " of less than {} points.".format(cube_data.shape, 
            np.prod(cube_data.shape), N) )
        gridInterpolator = RegularGridInterpolator( tuple(X_lin), cube_data, 
            method="linear", bounds_error=True )
        
        NX = np.asarray(cube_data.shape)
        r = (N / np.prod(cube_data.shape))**(1/dim)
        Nx = np.floor(r * NX)
        x_lin = []
        for i in range(0,dim):
            x_lin.append( np.linspace(0, cube_atoms.cell[i,i], Nx[i]) )
        x_list=np.meshgrid(*x_lin,indexing='ij')
        # X1d = n
        X = np.asarray(x_list)
        
        print(" Original grid: {:32,d} points, {} by spatial dimensions".format(int(NX.prod()), NX.astype(int)))
    
        cube_data = gridInterpolator( X.T )    
        print("Coarsened grid: {:32,d} points, {} by spatial dimensions".format(int(Nx.prod()), Nx.astype(int)))
    
        print("Coarsened grid by factor {:.3f}.".format(r))
        
        unit_cell = X[:,1,1,1].prod() # for uniform grid
        integral = cube_data.flatten().sum()*unit_cell
    
        print("on  original grid: data integral {:.3e} with {:.3e} unit cell".format(Integral, Unit_Cell))
        print("on coarsened grid: data integral {:.3e} with {:.3e} unit cell".format(integral, unit_cell))
    
    # atom positions in N x 3 array (rows x cols)
    pos = cube_atoms.get_positions()
    # potential on grid in M x 4 array (rows x cols)
    # dat = np.vstack( ( cube_data.flatten(), x_grid3.flatten()/Bohr, y_grid3.flatten()/Bohr, z_grid3.flatten()/Bohr ) ).T
    dat = np.vstack( ( cube_data.flatten(), X[0].flatten()/Bohr, X[1].flatten()/Bohr, X[2].flatten()/Bohr ) ).T
    
    n_atoms = len(cube_atoms)
    n_dat = dat.shape[0]
    l1 = "{:5d}{:6d}{:5d}".format(n_atoms,n_dat,0) # first line
    ghd = io.BytesIO() # .gesp header buffer
    np.savetxt(ghd, pos, fmt='%16.7E', delimiter='', header=l1, comments='')
    ghd = ('\n' + ' '*16).join(ghd.getvalue().decode().splitlines())
    np.savetxt(outfile_gesp, dat, fmt='%16.7E', delimiter='' ,header=ghd, comments='')


# source: https://swcarpentry.github.io/python-novice-inflammation/10-cmdline
# When you import a Python file, __name__ is the name of that file (e.g., when 
# importing readings.py, __name__ is 'readings'). However, when running a 
# script in bash, __name__ is always set to '__main__' in that script so that 
# you can determine if the file is being imported or run as a script.

if __name__ == '__main__':
    main()
