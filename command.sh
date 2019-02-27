#!/bin/bash -x 
./loop-esp-cost.sh --esp-infile-cube system100.vHtg_ua.cube \
  --dens-infile-cube system100.rho_ua.cube \
  --cost-outfile-hdf5 system100.cost.lnrhoref \
  --sign \
  --weights-outfile-cube system100.weights.lnrhoref \
  2>&1 | tee system100.loop-esp-cost.log
