# run Benzene DFT and compare
# HuLuYangs' cost functions 
# weighted by all-electron density
# and superposed single-atom electron densities

cd benzene/dft # from repository root

# interactive session
msub -I -l walltime=2:00:00 -l nodes=1:ppn=20
# job started ...
cd $MOAB_SUBMITDIR
module purge
module load gpaw
mpirun -n 20 gpaw-python benzene_dft_by_gpaw.py \
  2>&1 | tee benzene_dft_by_gpaw.log

# after successfull dft and creation of exp and density .cube files
logout # to end interactive session

cd ../.. # back to repository root
# if necessary
mkdir benzene/wdens_neg
mkdir benzene/nowdens_neg

# with weights from all electron density and inverted sign convention
./loop-esp-cost.sh --esp-infile-cube benzene/dft/benzene_vHtg.cube \
  --dens-infile-cube benzene/dft/benzene_rho.cube \
  --cost-outfile-hdf5 benzene/wdens_neg/benzene.cost.lnrhoref \
  --weights-outfile-cube benzene/wdens_neg/benzene.weights.lnrhoref \
  --sign 2>&1 | tee benzene/wdens_neg/benzene.loop-esp-cost.log

module purge
module load gpaw
./loop-fitESPforBenzene.sh -i benzene/nowdens/benzene.cost.lnrhoref \
  -o benzene/nowdens/benzene.charges.lnrhoref \
  2>&1 | tee benzene/nowdens/loop-fitESPforBenzene.log 

# with weights from superposed single-atom electron densities 
# and inverted sign convention
./loop-esp-cost.sh --esp-infile-cube benzene/dft/benzene_vHtg.cube \
  --cost-outfile-hdf5 benzene/nowdens_neg/benzene.cost.lnrhoref \
  --weights-outfile-cube benzene/nowdens_neg/benzene.weights.lnrhoref \
  --sign 2>&1 | tee benzene/nowdens_neg/benzene.loop-esp-cost.log

module purge
module load gpaw
./loop-fitESPforBenzene.sh -i benzene/nowdens_neg/benzene.cost.lnrhoref \
  -o benzene/nowdens_neg/benzene.charges.lnrhoref \
  2>&1 | tee benzene/nowdens_neg/loop-fitESPforBenzene.log

# evaluation of data in
# BenzeneTest.ipynb
