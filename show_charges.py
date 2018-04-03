#!/usr/bin/env python
""" Prints the point charges from a HORTON .h5 file """

import numpy as np

import h5py
import argparse

parser = argparse.ArgumentParser( \
        description='Prints the point charges from a HORTON .h5 file')
parser.add_argument('infile_h5', metavar='infile.h5', \
        help="Input file in HDF 5 format, usually the output 'charges.h5' of horton-esp-fit.py")
args = parser.parse_args()

hdf5_file = args.infile_h5

hdf5 = h5py.File(hdf5_file,'r')

print("Keys: %s" % hdf5.keys())
key0 = list(hdf5.keys())[0]

dat = list(hdf5[key0])

print("In '{}' we have...".format(key0))
print(dat)

