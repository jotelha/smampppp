#!/usr/bin/env python
""" Maps point charges obtained by GPAW and HORTON on the original'
        ' GROMACS topology initially modified by insertHbyList.py """
## jlh 2018/04/02

import ast
import h5py
import ase.io
from ase.io.cube import read_cube_data
import parmed as pmd
from parmed import gromacs
from insertHbyList import insertHbyList

import argparse

parser = argparse.ArgumentParser(\
    description='Maps point charges obtained by GPAW and HORTON on the original'
        ' GROMACS topology initially modified by insertHbyList.py')
#parser.add_argument('-c', '--charge',metavar='INTEGER_CHARGE',
#        type=int,nargs='?', const=1, default=0)
#parser.add_argument('infile', nargs='?')
parser.add_argument('infile_pdb', nargs='?', metavar='infile.pdb',
        default='system.pdb',
        help="Original .pdb file, before insertion of implicit hydrogen.")
parser.add_argument('infile_top', nargs='?', metavar='infile.top',
        default='system.top', help="Original GROMACS .top file")
parser.add_argument('infile_cube', nargs='?', metavar='infile.cube',
        default='esp.cube',
        help="ESP descrition (or other scalar field) in all-atom cube file.")
parser.add_argument('outfile_cube', nargs='?', metavar='outfile.cube', 
        default='esp_fitted_system.top', help="Output truncated by atoms only"
        "present in all-atoms description")
parser.add_argument('-i','--insertion-rules',
        default="{'CD4':1,'CD3':1,'CA2':2,'CA3':2,'CB2':2,'CB3':2}",
        help="A string representation of a python dictionary, describing how "
        "many implicit hydrogens have been inserted at which atom. Example: "
        "{'CD4':1,'CD3':1,'CA2':2,'CA3':2,'CB2':2,'CB3':2}")
args = parser.parse_args()

#implicitHbondingPartners={'CD4':1,'CD3':1,'CA2':2,'CA3':2,'CB2':2,'CB3':2}
print('Using replacement rules "{}"...'.format(args.insertion_rules))
implicitHbondingPartners = ast.literal_eval(args.insertion_rules)

# pdb_file = 'system100.pdb'
# top_file = 'system100.top'
# hdf5_file = 'smamp_charges.h5'

# top_outfile = 'smamp_esp_charges.top'

infile_pdb = args.infile_pdb
infile_top = args.infile_top
infile_cube  = args.infile_cube

outfile_cube = args.outfile_cube

ase_struct=ase.io.read(infile_pdb)
pmd_struct = pmd.load_file(infile_pdb)
pmd_top = gromacs.GromacsTopologyFile(infile_top,parametrize=False)
# throws some warnings on angle types, does not matter for bonding info
pmd_top.strip(':SOL,CL') # strip water and electrolyte from system
pmd_top.box = pmd_struct.box # Needed because .prmtop contains box info
pmd_top.positions = pmd_struct.positions

new_ase_struct, new_pmd_struct, names = insertHbyList(ase_struct,pmd_top,
        implicitHbondingPartners,1.0)


surplus_atoms = len(new_ase_struct) - len(ase_struct)
print("{} atoms are going to be truncated from file {}...".format(surplus_atoms,infile_cube))
# hdf5 = h5py.File(infile_h5,'r')
cube_data, cube_atoms = read_cube_data(infile_cube)
ase.io.write(outfile_cube, cube_atoms[1:len(ase_struct)], data=cube_data)
# ATTENTION: this script just truncates atoms based on total count difference
# in UA and AA representations
