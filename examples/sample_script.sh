# to be run on interactive NEMO job,
# invoked by the like of
# $ msub -I -l walltime=2:00:00 -l nodes=1:ppn=20 

module purge
module load gpaw/1.3.0

mpirun -n 20 gpaw-python ../traj_from_water_molecule.py h2o.555.traj
mpirun -n 20 gpaw-python ../gpw_from_traj.py -c 0 -b 5 5 5 h2o.555.traj h2o.555.gpw

# write ESP and electron densities to .cube files
mpirun -n 20 gpaw-python ../esp_from_gpw.py h2o.555.gpw h2o.vHtg.cube h2o.rho.cube h2o.rho_pseudo.cube

# RESP 2.1 ca n only process a maximum of 99,999 points
../cube2gesp.py h2o.vHtg.cube h2o.vHtg.gesp -n 99999
# output: 
#    Read .cube file of shape (48, 48, 48) in box [ 4.99998803  4.99998803  4.99998803].
#    Reconstructed uniform grid of shape (3, 48, 48, 48).
#    Interpolate .cube data of shape (48, 48, 48) and 110592 points in total onto grid of less than 99999 points.
#     Original grid:                          110,592 points, [48 48 48] by spatial dimensions
#    Coarsened grid:                           97,336 points, [ 46.  46.  46.] by spatial dimensions
#    Coarsened grid by factor 0.967.
#    on  original grid: data integral -3.9730541047951107 with 0.0012039634944905316 unit cell
#    on coarsened grid: data integral -4.00725739186612 with 0.0013717322566638183 unit cell

# Amber's RESP now requires quite a bit of 'format juggling':
../traj2pdb.py h2o.555.traj h2o.pdb

module load amber/16
antechamber -i h2o.pdb -fi pdb -o h2o.ac -fo ac
respgen -i  h2o.ac -o h2o.resp1.in -e 2 -f resp1

resp -O -i h2o.resp1.in -e h2o.vHtg.gesp -o h2o.resp1.out \
    -p h2o.resp1.pch -t h2o.resp1.chg -s h2o.resp1.gesp

# TODO: evaluate results and compare to reference case in examples/res,
# from AmberTools/examples/resp_charge_fit/water


# for comparison, on coarser grid
mkdir 999
cd 999
../../cube2gesp.py ../h2o.vHtg.cube h2o.vHtg.gesp.999 -n 999  
#    Read .cube file of shape (48, 48, 48) in box [ 4.99998803  4.99998803  4.99998803].
#    Reconstructed uniform grid of shape (3, 48, 48, 48).
#    Interpolate .cube data of shape (48, 48, 48) and 110592 points in total onto grid of less than 999 points.
#     Original grid:                          110,592 points, [48 48 48] by spatial dimensions
#    Coarsened grid:                              729 points, [ 9.  9.  9.] by spatial dimensions
#    Coarsened grid by factor 0.208.
#    on  original grid: data integral -3.9730541047951107 with 0.0012039634944905316 unit cell
#    on coarsened grid: data integral -4.041950511296964 with 0.24413887087595793 unit cell
resp -O -i ../h2o.resp1.in -e h2o.vHtg.gesp.999 -o h2o.resp1.out \
    -p h2o.resp1.pch -t h2o.resp1.chg -s h2o.resp1.gesp

cd ..
mkdir 9999
cd 9999
../../cube2gesp.py ../h2o.vHtg.cube h2o.vHtg.gesp.9999 -n 9999
#    Read .cube file of shape (48, 48, 48) in box [ 4.99998803  4.99998803  4.99998803].
#    Reconstructed uniform grid of shape (3, 48, 48, 48).
#    Interpolate .cube data of shape (48, 48, 48) and 110592 points in total onto grid of less than 9999 points.
#     Original grid:                          110,592 points, [48 48 48] by spatial dimensions
#    Coarsened grid:                            9,261 points, [ 21.  21.  21.] by spatial dimensions
#    Coarsened grid by factor 0.449.
#    on  original grid: data integral -3.9730541047951107 with 0.0012039634944905316 unit cell
#    on coarsened grid: data integral -3.9580503300627727 with 0.015624887736061303 unit cell

# Fitting results for grid of 9999 points

#    Resp charges for organic molecule
#    
#    iqopt   irstrnt  ihfree     qwt
#      0      1      1        0.000500
#    
#     rel.rms   dipole mom       Qxx      Qyy      Qzz
#        0.87731   0.08685   0.18590   0.31757  -0.50347
#    
#              Point charges before & after optimization
#        NO   At.No.    q0           q(opt)   IVARY  d(rstr)/dq
#         1     8   0.000000      -0.056815    0    0.004347
#         2     1   0.000000       0.028407    0    0.000000
#         3     1   0.000000       0.028407    2    0.000000
#    
#            Statistics of the fitting:
#      The initial sum of squares (ssvpot)                      0.005
#      The residual sum of squares (chipot)                     0.004
#      The std err of estimate (sqrt(chipot/N))               0.00210
#      ESP relative RMS (SQRT(chipot/ssvpot))                 0.87731
#    
#     Dipole Moment (Debye)=   0.08685

# and for 99999 points

#    Resp charges for organic molecule
#    
#    iqopt   irstrnt  ihfree     qwt
#      0      1      1        0.000500
#    
#     rel.rms   dipole mom       Qxx      Qyy      Qzz
#        0.81024   0.13846   0.29637   0.50630  -0.80268
#    
#              Point charges before & after optimization
#        NO   At.No.    q0           q(opt)   IVARY  d(rstr)/dq
#         1     8   0.000000      -0.090578    0    0.003706
#         2     1   0.000000       0.045289    0    0.000000
#         3     1   0.000000       0.045289    2    0.000000
#    
#            Statistics of the fitting:
#      The initial sum of squares (ssvpot)                      0.131
#      The residual sum of squares (chipot)                     0.086
#      The std err of estimate (sqrt(chipot/N))               0.00297
#      ESP relative RMS (SQRT(chipot/ssvpot))                 0.81024
#    
#     Dipole Moment (Debye)=   0.13846

# apparently agree UP TO A SCALING FACTOR.

# Conclusion, as of 2018/04/13:
#   TODO1:  Double-check grid interpolation 
#     (although integrals agree after interpolation, as shown in "cube2gesp.py"
#     outputs above).
#   TODO2:  Double-check required input units by RESP.f, although code comments,
#     explicitly mention A.U., which I would interpret as Bohr for distance
#     and Hartree (per elementary charge) for potential.


# $ module load resp/2.4 # overrides AmbertTool-native RESP 2.1 by newer verson
# should be able to process a maximum of 999,999 points.
# $ ../cube2gesp.py vHtg.cube vHtg.gesp.999999 -n 999999
# However, apparently requires modification to input file h2o.resp1.in.
#   TODO3: Get input right for 2.4 to work an run on
# $ resp -O -i h2o.resp1.in -e vHtg.gesp.999990 -o h2o.resp1.out.999999 \
# >    -p h2o.resp1.pch.999999 -t h2o.resp1.chg.999999
# Refer to README via
# $ less "$(dirname $(which resp))/README-2.4.txt"
