#!/usr/bin/env python
# gesp2cube, N atoms, M grid points
# from ase.io.cube import read_cube_data
from ase.io import write
from ase.units import Bohr
from ase import Atoms
#from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import NearestNDInterpolator
import numpy as np
# import io

import argparse

def main():
    parser = argparse.ArgumentParser(description='Convert'
            ' Fortran-RESP-readible .gesp format. to .cube format')
    parser.add_argument('infile_gesp', metavar='infile.gesp',
            help="Input in .gesp format")
    parser.add_argument('outfile_cube', nargs='?', metavar='outfile.cube',
            default='outfile.cube', help="Output in .cube format"
            " (A.U., i.e. Bohr and Hartree)")
    parser.add_argument('-n', nargs='?',
            metavar='N', default=-1, const=99999, type=int,
            help="Maximum number of grid points")
    args = parser.parse_args()

    gesp2cube(args.infile_gesp, args.outfile_cube, args.n)
    
def gesp2cube(infile_gesp,outfile_cube,N=-1,method='nearest'):
    # dat = np.vstack( ( cube_data.flatten(), X[0].flatten()/Bohr, X[1].flatten()/Bohr, X[2].flatten()/Bohr ) ).T
    f = open(infile_gesp, "r")
    l1 = f.readline() # read header
    l1words = l1.split() # header fields separated by white space
    n_atoms = int(l1words[0]) # header format: 1st field number of atoms
    n_dat = int(l1words[1]) # 2nd field number of data points
    
    l_atoms = [ f.readline() for i in range(n_atoms) ] # list of atoms
    # s_atoms = ''.join(l_atoms)
    p_atoms = np.loadtxt(l_atoms) # position of atoms
    
    l_dat = f.readlines()
    dat = np.loadtxt(l_dat) # all remaining lines 
    
    if len(dat) != n_dat:
        raise Exception("Number of expected data points {:d} does not agree with"
            " actually read amount {:d}!".format(n_dat,len(dat)))
    
    gesp_data = dat[:,0] # first valua of each line is actual (potential) value
    X = dat[:,1:] * Bohr # reamaining values are grid points, from Bohr to Ang
    X_min = X.min(axis=0) # minimum coordinates in each spatial direction
    X_max = X.max(axis=0) # dito for maximum
    box = X_max - X_min
    vol = box.prod()
    
    # just in case: correct coordinates by offset 
    X = X - X_min
    p_atoms = p_atoms - X_min
    dim = X.shape[1] # spatial dimensions, usually three
    
    # cube file measures in Bohr, but ASE converts input to Angstrom
    # however, field is not converted while distances are
    # gesp_data, cube_atoms = read_gesp_data(args.infile_cube)
    
    # nX = gesp_data.shape
    # X = cube_atoms.cell.diagonal() # measures of cell in spatial dimensions
    # x_grid = np.linspace(0,X[0],nX[0])
    # y_grid = np.linspace(0,X[1],nX[1])
    # z_grid = np.linspace(0,X[2],nX[2])
    
    # x_grid3,y_grid3,z_grid3=np.meshgrid(x_grid,y_grid,z_grid)
    
    # general approach for creating a grid of coordinates
    # in 3 dimensions and N points in each spatial dimension,
    # X.shape == (3, N, N, N)
    
    print("Read .gesp file of {} atoms and {} data points in box {} Ang of volume"
        " {:.3e} Ang^3.".format(n_atoms,n_dat,box,vol) )
    
    
    if N > 0:
        n_max = N
        print("Specified maximum number of grid points as {} by option -n"
            .format(n_max))
    else:
        n_max = n_dat
        print("Number of grid points will be less or equal to number of .gesp data"
            " points {}".format(n_max))
    #r = (vol / (n_dat**(1/dim)-1)**dim))**(1/dim) # determine uniform grid spacing r
    r = (vol / n_max)**(1/dim) # determine uniform grid spacing r
    XN = np.floor((box/r)) # ratio gives #volume elements, add 1 for #grid points
    
    X_lin = []
    for i in range(0,dim):
        X_lin.append( np.linspace(0, box[i], XN[i] ) )
    
    X_list=np.meshgrid(*X_lin,indexing='ij')
    X_grid = np.asarray(X_list)
    Unit_Cell = X_grid[:,1,1,1].prod() # for uniform grid, or use just X better ?
    Integral = gesp_data.flatten().sum()*Unit_Cell
    
    print("Constructed uniform grid of shape {} with spacing {:.3e}."
            .format(X_grid.shape,r))
    #print("On grid: data integral {} with {} Ang^3 unit cell"
    #        .format(Integral, Unit_Cell))
    
    print("Interpolate .gesp data at {:d} points in total onto grid"
            " of {:d} <= {:d} points.".format(n_dat, np.prod(X_grid[0,:].shape), n_max))
    
    if method == 'nearest':
        print("Nearest value interpolation.")
        nearestInterpolator = NearestNDInterpolator( X, gesp_data )    
        print("...")
        cube_data = nearestInterpolator( X_grid.T )
    else:
        print("Linear interpolation.")
        linearInterpolator = LinearNDInterpolator( X, gesp_data, fill_value=0.0 )    
        print("...")
        cube_data = linearInterpolator( X_grid.T )

    
    unit_cell = X_grid[:,1,1,1].prod()
    integral = cube_data.flatten().sum()*unit_cell
    
    print("On original .gesp space (assuming uniform volume discretization):"
            "data integral {:.3e} with {:.3e} Ang^3 unit cell".format(Integral, Unit_Cell))
    print("On regular .cube grid: data integral {:.3e} with {:.3e} unit cell".format(integral, unit_cell))
    
    
    atoms = Atoms(positions=p_atoms,cell=box)
    write(outfile_cube, atoms, data=cube_data, format='cube') 

if __name__ == '__main__':
    main()


