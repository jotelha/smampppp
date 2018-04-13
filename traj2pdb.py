#!/usr/bin/env python
import ase.io
import argparse
# import parmed as pmd

parser = argparse.ArgumentParser(
        description='Convert .traj format to .pdb (or other ASE-output'
        ' supported) format.')
parser.add_argument('infile_traj', metavar='infile.traj')
parser.add_argument('outfile_pdb', metavar='outfile.pdb', 
        help="Extension determines output format")
args = parser.parse_args()

struc = ase.io.read(args.infile_traj)
struc.write(args.outfile_pdb)

#pmd_struct = pmd.load_file(args.infile_traj)
#pmd_struct.write_pdb(args.outfile_pdb)
