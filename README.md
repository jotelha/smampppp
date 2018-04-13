# Overview

## Units & File formats

Native units in
  * ASE
    * [length] = Angstrom
    * [energy] = eV
    * [electron density] = Angstrom^-3
  * GPAW (internally): "atomic units", i.e. in gpaw.hamiltonian
    * [length] = Bohr
    * [energy] = Hartree
    * [electron density] = Bohr^-3
  
Both systems measure charge in e = 1 (unit charge), thus [electostatic 
potential] = [energy]/[charge] same as unit of energy)

The GAUSSION's .cube file format uses same atomic unit system as GPAW 
internally. Thus, ASE data has to be converted before outputting to 
.cube-file. See "esp_from_gpw.py" for examples.

Horton and RESP internally uses the same units, i.e. Bohr and Hartree.
Horton reads .cube files as they are, without performing any conversion.

## Potentials

GPAW offers one function `get_electrostatic_potential()` within gpaw/paw.py
evaluating the variable `vHt_g` of a calculator's hamiltonian by calling
`calculator.hamiltonian.get_electrostatic_potential(...)`

# ESP fitting

## by Horton

TODO: fill in

## RESP with Antechamber

on the example of water. RESP fitting based on

From a PDB (or other structure format description)
such as `h2o.pdb`
    
    ATOM          1  OH  OSP3    1       4.013   0.831  -9.083  1.00  0.00
    ATOM      2 1HH  OSP3    1       4.941   0.844  -8.837  1.00  0.00
    ATOM      3 2HH  OSP3    1       3.750  -0.068  -9.293  1.00  0.00
    TER

 
generate an antechamber input file `h2o.ac` 

    CHARGE      0.00 ( 0 )
    Formula: H2 O1
    ATOM      1  OH  OSP     1       4.013   0.831  -9.083  0.000000        oh
    ATOM      2  HH1 OSP     1       4.941   0.844  -8.837  0.000000        ho
    ATOM      3  HH2 OSP     1       3.750  -0.068  -9.293  0.000000        ho
    BOND    1    1    2    1     OH  HH1
    BOND    2    1    3    1     OH  HH2

via

    antechamber -i h2o.pdb -fi pdb -o h2o.ac -fo ac
    
Second, generate an RESP input file 

    Resp charges for organic molecule

     &cntrl

     nmol = 1,
     ihfree = 1,
     ioutopt = 1,
     qwt = 0.00050,

     &end
        1.0
    Resp charges for organic molecule
        0    3
        8    0
        1    0
        1    2

with

    respgen -i  h2o.ac -o h2o.resp1.in -e 2 -f resp1

which also offers the additional options

    Usage: respgen -i input file name(ac)
               -o output file name
               -l maximum path length (default is -1, only recommand to use
                  when the program takes long time to finish or causes core dump.)
                  If applied, a value of 8 to 10 should good)
               -f output file format (resp1 or resp2) 
                  resp0 - evaluation the current charges 
                  resp1 - first stage resp fitting 
                  resp2 - second stage resp fitting
                  iresp0 - evaluation the current charges for polarizable model
                  iresp1- first stage of i_resp fitting 
                  iresp2- second stage of i_resp fitting
                  resp3 - one-stage resp fitting
                  resp4 - calculating ESP from point charges
                  resp5 - no-equalization
               -e equalizing atomic charge, default is 1
                  0 not use 
                  1 by atomic paths
                  2 by atomic paths and structural information, i.e. E/Z confirgurations
               -a additional input data (predefined charges, atom groups etc).)
               -n number of conformations (default is 1)
               -w weight of charge constraint, in default, 0.0005 for resp1 and 0.001 fore resp2

With such an input file, run

    resp -O -i h2o.resp1.in -e h2o.gesp -o h2o.resp1.out -p h2o.resp1.pch -t h2o.resp1.chg

which performs the actual ESP fitting. `h2o.resp1.out` is the log file,
`h2o.resp.pch` contains formatted information on point charges before and
after fitting, and `h2o.resp1.chg` just contains the fitted point charges.

`h2o.gesp` is the file containing information about the ESP on a grid.
The format requirements are strict, as specified in the source `resp.f`

    C Unit 10 input of ESP's  (mandatory)
    C
    CC      natoms,nesp (2i5)
    C      natoms,nesp (i5,i6)
    C              X , Y , Z  .   FORMAT (17X,3E16.7)
    C      QUPOT , X , Y , Z  .   FORMAT (1X,4E16.7)
    C
    C          QUPOT = THE QUANTUM MECHANICAL ELECTROSTATIC
    C                  POTENTIAL ( A.U )
    C
    C          X,Y,Z = THE COORDINATE AT WHICH THE POTENTIAL
    C                  IS CALCULATED ( A.U )
    C
    C      NOTE : THE PROGRAM G80UCSF WRITES IN THIS FORMAT BUT THE
    C             OUTPUT OF G90 MUST BE TRANSLATED (PROGRAM BOHR).

This header excerpt originates from `q4md/resp-2.4` by 
http://q4md-forcefieldtools.org, 2013. Compaired with
the 1994 RESP v2.1 as contained within the official AmberTools17 
`amber16/AmberTools/src/etc/resp.F` by Caldwell, J. and Bayly, C.,
it allows for ESP data on M_max = 999,999 grid points (i6) instead of 
previously 99,999 points (i5) and introduces some other modifications.

In the following an example from 
amber16/AmberTools/examples/resp_charge_fit/water/esp_wat.dat:

       3  295    0
                      -0.9982123E-32   0.0000000E+00   0.2313846E+00
                      -0.2610123E-32   0.1494187E+01  -0.9255383E+00
                      -0.1829851E-15  -0.1494187E+01  -0.9255383E+00
      -0.4207115E-01  -0.9982123E-32   0.0000000E+00   0.3935248E+01
      -0.4040255E-01   0.1851931E+01   0.0000000E+00   0.3439024E+01
      -0.3659220E-01   0.9259657E+00   0.1603820E+01   0.3439024E+01
    ...

The first line contains the number N of atoms (or point charges) and the number
M of entries for describing the potential on a grid. The zero is neglected.
Next follow the spatial coordinates of all N atoms. Eventually M lines
each give a potential value (1st column) at some point in space (2nd - 4th).
All non-integer numbers are of format %16.7E. The following python snippet 

    # cube2gesp, N atoms, M grid points
    from ase.io.cube import read_cube_data
    from ase.units import Bohr
    import numpy as np
    import io

    # cube file measures in Bohr, but ASE converts input to Angstrom 
    # however, field is not converted while distances are
    cube_data, cube_atoms = read_cube_data('infile.cube')
    
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
    
    n_atoms = len(struc)
    n_dat = dat.shape[0]
    l1 = "{:5d}{:6d}{:5d}".format(n_atoms,n_dat,0) # first line
    ghd = io.BytesIO() # .gesp header buffer 
    np.savetxt(ghd, pos, fmt='%16.7E', delimiter='', header=l1, comments='')
    ghd = ('\n' + ' '*16).join(ghd.getvalue().decode().splitlines())
    np.savetxt('outfile.gesp', dat, fmt='%16.7E', delimiter='' ,header=ghd, comments='')

illustrates how to convert a .cube file to the desired .gesp format by means of
ASE and numpy. Note that the number of points must not exceed M_max.
This conversion is automized by `cube2gesp.py`, 

    usage: cube2gesp.py [-h] [-n [N]] infile.cube [outfile.gesp]
    
    Convert .cube format to Fortran-RESP-readible .gesp format.
    
    positional arguments:
      infile.cube   Input in .cube format (Bohr and Hartree)
      outfile.gesp  Output in .gesp format
    
    optional arguments:
      -h, --help    show this help message and exit
      -n [N]        Maximum number of grid points

which accounts for the hard limit on point numbers in RESP.f by interpolating
.cube file data onto a grid of at most N points if requested.

AmberTools (with antechamber) and RESP are loadable as modules on our NEMO 
group environment. See `examples/sample_script.sh` for RESP fitting workflow on
potentials generated by GPAW. 

Attention: As of 2018/04/13, a scaling factor artifact appears within the 
fitted point charges, apparently depending on the ESP data's grid spacing. 
Not yet fixed.



