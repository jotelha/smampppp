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

INSERT SHORT PARAGRAPH ON WHICH INTERNAL GPAW VARIABLE REPRESENTS WHAT KIND OF POTENTIAL

# ESP fitting

## by GPAW and Horton

The following bash script snippet illustrates how to extract ESP and
electron density from a GPAW DFT run and feed it to Horton:

```
#!/bin/bash -x
# sample work flow on NEMO for extracting ESP and electron density 
# from a GPAW  simulation

# run (parallel mpirun) gpaw and (non-parallel, but resource-intense)
# horton-esp-cost.py steps as jobs in interactive session started by

msub -I -l walltime=2:00:00 -l nodes=1:ppn=20 # start interactive 2h session

module purge
module load gpaw/1.3.0
# extract .gpw file from SMAMP molecule of charge +6
mpirun -n 20 gpaw-python gpw_from_traj.py -c 6 molecule.traj smamp.gpw
# usage: gpw_from_traj.py [-h] [-c [CHARGE]] [-b X Y Z] infile_traj outfile_gpw
#
# Reconstructs quantum system from GPAW .traj file given as command line
# argument and stores it in .gpw format.
# 
# positional arguments:
#   infile_traj
#   outfile_gpw
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   -c [CHARGE], --charge [CHARGE]
#                         The system's total charge, 0.0 if not specified, 1.0
#                         if specified without value.
# -b X Y Z, --box X Y Z
#                         The bounding box' measures in Angstrom, default [25.0,
#                         25.0, 25.0]

# write ESP and electron densities to standard files
mpirun -n 20 gpaw-python esp_from_gpw.py smamp.gpw
# usage: esp_from_gpw.py [-h] infile [vHtg.cube] [rho.cube] [rho_pseudo.cube]
#
# Extracts the electrostatic potential (ESP) from a GPAW .gpw restart file given
# as command line argument.
#
# positional arguments:
#   infile
#   vHtg.cube        Electrostatic potential (vHt_g) in GAUSSIAN-native .cube
#                    format, default 'vHtg.cube'
#   rho.cube         All-electron density in GAUSSIAN-native .cube format,
#                    default 'rho.cube'
#   rho_pseudo.cube  Pseudo electron density in GAUSSIAN-native .cube format,
#                    default'rho_pseudo.cube'
#
# optional arguments:
#   -h, --help       show this help message and exit

# with all-electron density for weighting:
# construct united-atom cube files by truncating implicit atoms:
./aa2ua_cube.py system100.pdb system100.top vHtg.cube vHtg_ua.cube
./aa2ua_cube.py system100.pdb system100.top rho.cube rho_ua.cube
# the original gromacs .pdb and .top files are required to count 
# the atmons to truncate for UA representation

# construct cost function with weighting from electron density:
module purge
module load horton/2.1.0b
horton-esp-cost.py vHtg_ua.cube smamp_esp_cost_wdens_ua.h5 --pbc 000 --wdens rho_ua.cube --overwrite
horton-esp-fit.py -q 6 smamp_esp_cost_wdens_ua.h5 smamp_esp_charges_wdens_ua.h5 --overwrite
# horton-esp-fit.py possesses the possibility to impose symmetry constraints via
#   --ridge RIDGE         The thikonov regularization strength used when solving
#                         the charges. [default=0.0]
#   --symmetry SYMMETRY SYMMETRY
#                         Perform a symmetry analysis on the charges. This
#                         option requires two arguments argument in the
#                         following order: the cube file used to construct the
#                         ESP cost function and a CIF file with the generators
#                         of the symmetry of this system and a primitive unit
#                         cell.

# TODO: we have to figure out how to impose the symmetry constraints here
# if --overwrite not set, HORTON will not update possibly existin output files

# map fitted charges back onto the original united-atoms in GROMACS topology file
module purge
module load gpaw/1.3.0
./mapCharges.py system100.pdb system100.top smamp_esp_charges_wdens_ua.h5 smamp_esp_charges_wdens_ua.top
# usage: mapCharges.py [-h] [-i INSERTION_RULES] [-aa]
#                      [infile.pdb] [infile.top] [infile.h5] [outfile.top]
# 
# Maps point charges obtained by GPAW and HORTON on the original GROMACS
# topology initially modified by insertHbyList.py
# 
# positional arguments:
#   infile.pdb            Original .pdb file, before insertion of implicit
#                         hydrogen.
#   infile.top            Original GROMACS .top file
#   infile.h5             Point charges in hdf5 format by horton-esp-fit.py
#   outfile.top           GROMACS .top output filewith updated charges according
#                         to given .hdf5
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   -i INSERTION_RULES, --insertion-rules INSERTION_RULES
#                         A string representation of a python dictionary,
#                         describing how many implicit hydrogens have been
#                         inserted at which atom. Example:
#                         {'CD4':1,'CD3':1,'CA2':2,'CA3':2,'CB2':2,'CB3':2}
#   -aa, --all-atoms      Determines whether the charges are for all-atoms (or
#                         united-atoms, default) representation
# 

# similar procedure for an all-atom fit, just using the unmodified cube files:
module purge
module load horton/2.1.0b3 
horton-esp-cost.py vHtg.cube smamp_esp_cost_wdens.h5 --pbc 000 --wdens rho.cube --overwrite
horton-esp-fit.py -q 6 smamp_esp_cost_wdens.h5 smamp_esp_charges_wdens.h5 --overwrite

module purge
module load gpaw/1.3.0
./mapCharges.py -aa system100.pdb system100.top smamp_esp_charges_wdens.h5 smamp_esp_charges_wdens.top

# attention: since horton is based on python 2.7 while all other software 
# depends on python 3.6, the module purgin and reloading is necessary
# if steps are carried out sequentially within one session as shown here
```

### Looping over parameters

For Benzene, first, use

```
./loop-esp-cost.sh --esp-infile-cube benzene/dft/benzene_vHtg.cube \
  --cost-outfile-hdf5 benzene/nowdens_neg/benzene.cost.lnrhoref \
  --weights-outfile-cube benzene/nowdens_neg/benzene.weights.lnrhoref \
  --sign 2>&1 | tee benzene/nowdens_neg/benzene.loop-esp-cost.log
```

to loop over the parameter ln(rho_ref) as defined within the script, then use

```
./loop-fitESPforBenzene.sh -i benzene/nowdens_neg/benzene.cost.lnrhoref \
  -o benzene/nowdens_neg/benzene.charges.lnrhoref \
  2>&1 | tee benzene/nowdens_neg/loop-fitESPforBenzene.log 
```

to fit charges for all cost function filess.

In general, use

```
./loop-esp-cost.sh --esp-infile-cube sandbox/system100.vHtg_ua.cube \
  --dens-infile-cube sandbox/system100.rho_ua.cube \
  --cost-outfile-hdf5 sandbox/wdens/system100.cost.lnrhoref \
  --weights-outfile-cube sandbox/wdens/system100.weights.lnrhoref \
  --sign > sandbox/wdens/system100.loop-esp-cost.log 2>&1  &
```

for looping over the reference density logarithm, constructing weights based upon the 
true all-electron density stored in `sandbox/system100.rho_ua.cube`. Otherwise,
leave out this parameter to construct HORTON's "pro-density" by superposing
single-atom electron densities taken from internal precomputed tables:

```
./loop-esp-cost.sh --esp-infile-cube sandbox/system100.vHtg_ua.cube \
  --cost-outfile-hdf5 sandbox/nowdens/system100.cost.lnrhoref \
  --weights-outfile-cube sandbox/nowdens/system100.weights.lnrhoref \
  --sign > sandbox/nowdens/system100.loop-esp-cost.log 2>&1  &
```

Afterwards, use a command like

```
./loop-fitESPconstrained.sh -i sandbox/nowdens/system100.cost.lnrhoref \
  -o sandbox/nowdens/system100.fit -q 6 \
  > sandbox/nowdens/system100.loop-fitESPconstrained.log 2>&1 &
```

to fit point charges for all cost function files under toal charge -q constraint.
      
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
group environment. See [examples/sample_script.sh](./examples/sample_script.sh) 
for RESP fitting workflow on potentials generated by GPAW. 

Attention: As of 2018/04/13, a scaling factor artifact appears within the 
fitted point charges, apparently depending on the ESP data's grid spacing. 
Not yet fixed.
