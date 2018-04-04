## jlh 2018/02/13
# takes dictionary of united atoms and implicit hydrogen number,
# adding a certain number of explicit hydrogens around these.
# Needs and returns corresponding ASE and ParmEd representations of the system

# TODO: think about the case when one has to add more than two H-atoms


import numpy as np
import matscipy as msp
from matscipy.neighbours import neighbour_list
from ase.data import atomic_numbers
import ase.io
from ase.neighborlist import NeighborList
from ase.visualize import view
import parmed as pmd
from parmed import gromacs
import nglview as nv
import sys

from insetHbyList import insertHbyList

# define a dictionary, marking how many H atoms are missing at which
# bonding partner explicitly:
implicitHbondingPartners={'CD4':1,'CD3':1,'CA2':2,'CA3':2,'CB2':2,'CB3':2}
# every occurence of these atoms (in sample case 3 times for each residue)
# will be processed independently

#loop over all MD steps
for i in range(100,101):
    print('test i is {}, or '.format(i), i)
    ase_struct=ase.io.read('system{}.pdb'.format(i))
    pmd_struct = pmd.load_file('system{}.pdb'.format(i))
    pmd_top = gromacs.GromacsTopologyFile('system{}.top'.format(i),parametrize=False)
    # throws some warnings on angle types, does not matter for bonding info
    pmd_top.strip(':SOL,CL') # strip water and electrolyte from system

    pmd_top.box = pmd_struct.box # Needed because .prmtop contains box info
    pmd_top.positions = pmd_struct.positions

    # placing choice not ideal yet, other tools such as VMD or Avogadro do
    # not necessarily infer correct bonding for new H-atoms...
    new_ase_struct, new_pmd_top=insertHbyList(ase_struct,pmd_top,implicitHbondingPartners,1.0)
    new_ase_struct.write('ase_pdbH_{}.pdb'.format(i))
    new_ase_struct.write('ase_pdbH_{}.traj'.format(i))

    #visualize the atom structure in ase-gui
    #view(new_ase_struct)

    # ... hence we can use an explicit topology file to visualize connectivity as "understood" by this tool
    new_pmd_top.write_pdb('pmd_pdbH_{}.pdb'.format(i))
    test_pmd = pmd.load_file('pmd_pdbH_{}.pdb'.format(i))

    new_pmd_top.write_psf('pmd_pdbH_{}.psf'.format(i)) # some topology format, un functionality similar to GROMACS' .top, but readable by VMD
    # in VMD, first load .psf, then .pdb to suppress (wrong/unwanted) connectivity inference from distances and display
    # bonds as explicitly defined.
