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
#gromacs.GROMACS_TOPDIR = "/home/jotelha/gromacs/2016.4/share/gromacs/top"

def insertHbyList(ase_struct,pmd_top,implicitHbondingPartners,bond_length=1.0):
    # make copies of passed structures as not to alter originals:
    new_pmd_top = pmd_top.copy(pmd.Structure)
    new_ase_struct = ase_struct.copy()
    # make copied atoms accessible by unchangable indices (standard list)
    originalAtoms = [a for a in new_pmd_top.atoms]

    implicitHbondingPartnersIdxHnoTuples = [(a.idx,implicitHbondingPartners[k])
                                            for a in pmd_top.atoms
                                            for k in implicitHbondingPartners.keys()
                                            if a.name == k]
    implicitHbondingPartnersIdxHnoDict = dict(implicitHbondingPartnersIdxHnoTuples)

    # build numbered neighbour list "manually"
    #i: list of atoms e.g. [0,0,0,1,1,2,2,...]
    i = np.array([b.atom1.idx for b in new_pmd_top.bonds])
    #j: list of bonding partners corresponding to i [1,2,3,...] ==> 0 has bondpartners 1,2 and 3; etc.
    j = np.array([b.atom2.idx for b in new_pmd_top.bonds])
    r = new_ase_struct.positions

    for k,Hno in implicitHbondingPartnersIdxHnoDict.items(): # for all atoms to append hydrogen to
        print('Adding {} H-atoms to {} (#{})...'.format(Hno,originalAtoms[k].name,k))
        for h in range(0,Hno):
            r = new_ase_struct.positions
            bondingPartners = j[i==k]
            print('bondingPartners', bondingPartners)
            partnerStr=''
            for p in bondingPartners:
                if partnerStr == '':
                    partnerStr = originalAtoms[p].name
                else:
                    partnerStr += ', ' + originalAtoms[p].name
            print('Atom {} already has bonding partners {}'.format(originalAtoms[k].name,partnerStr))
            dr = (r[j[i==k]] - r[k]).mean(axis=0)
            # my understanding: dr is vector
            # from atom k's position towards the geometrical center of mass
            # it forms with its defined neighbours
            # r0 is a vector offset into the opposit direction:
            dr = dr/np.linalg.norm(dr) #normalized vector in direction dr

            #calculate an orthogonal vector 'dr_ortho' on dr
            #and push the H atoms in dr+dr_ortho and dr-dr_ortho
            #if one has to add more than two H atoms introduce dr_ortho_2 = dr x dr_ortho
            dr_ortho = np.cross(dr,np.array([1,0,0]))
            if np.linalg.norm(dr_ortho) < 0.1 :   #if dr and (1,0,0) have almost the same direction
                dr_ortho = np.cross(dr,np.array([0,1,0]))
            # (1-2*h) = {1,-1} for h={0,1}
            h_pos_vec = (dr+(1-2*h)*dr_ortho)/np.linalg.norm(dr+(1-2*h)*dr_ortho)
            r0 = r[k] - bond_length*h_pos_vec

            new_ase_struct += ase.Atom('H', r0) # add atom in ase structure
            n_atoms = len(new_ase_struct)

            #introduce a corrector step for a added atom which is too close to others
            #do as many corrector stepps until all atoms are more than 1\AA appart
            c_step = 0
            while True:
                nl = NeighborList(cutoffs=[.5]*len(new_ase_struct),
                                  skin=0.09, self_interaction=False,
                                  bothways=True)
                nl.update(new_ase_struct)
                indices, offsets = nl.get_neighbors(-1)
                indices = np.delete(indices, np.where(indices==k))
                if len(indices) == 0:
                    break

                elif c_step > 15:
                    print('programm needs more than 15 corrector steps for H atom {} at atom {}'.format(n_atoms, k))
                    sys.exit(15)
                    break

                print('too close atoms', indices)
                c_step += 1
                print('correcter step {} for H {} at atom {}'.format(c_step, n_atoms-1, k))
                # if indices not empty -> the atom(-1)=r_H is to close together
                # with atom a_close=indices[0], it is a H-atom belonging to atom 'k'=r_k .
                #correctorstep: corr_step = (r_H-a_close)/|(r_H-a_close)|
                #corrected_pos: corr_pos = ((r_H-r_k) + corr_step)/|((r_H-r_k) + corr_step)|
                #new H position: new_r_H = r_k + corr_pos
                r_H, r_k, a_close = np.take(new_ase_struct.get_positions(),[-1,k,indices[0]], axis=0)
                #print('r_H, r_k, a_close', r_H, r_k, a_close)
                corr_step = (r_H-a_close)/np.linalg.norm((r_H-a_close)) #maybe introduce here a skaling Faktor s=0.3 or somthing like that to make tiny corrections and don't overshoot.
                corr_pos = ((r_H-r_k) + corr_step)/np.linalg.norm((r_H-r_k) + corr_step)
                new_r_H = r_k + bond_length * corr_pos
                #correct the H position to new_r_H in new_ase_struct
                trans = np.zeros([n_atoms,3])
                trans[-1] = new_r_H - r_H
                new_ase_struct.translate(trans)

            #view(new_ase_struct)
            #sys.exit()

            i = np.append(i,k) # manually update numbered neighbour lists
            j = np.append(j,len(new_ase_struct)-1)

            # update pmd topology
            bondingPartner = originalAtoms[k] # here we need the original numbering,
            # as ParmEd alters indices when adding atoms to the structure
            nameH = '{}{}'.format(h+1,bondingPartner.name) # atom needs a unique name
            print('Adding H-atom {} at position [ {}, {}, {} ]'.format(nameH,r0[0], r0[1], r0[2]))
            new_H = pmd.Atom(name=nameH, type='H', atomic_number=1)
            new_H.xx = r0[0] # ParmEd documentation not very helpful, did not find any more compact assignment
            new_H.xy = r0[1]
            new_H.xz = r0[2]
            # do not understand ParmEd that well, apparently we need the Bond object in order to update topology
            new_Bond = pmd.Bond(bondingPartner,new_H)
            new_H.bond_to(bondingPartner) # not sure, whether this is necessary
            new_pmd_top.bonds.append(new_Bond)
            new_pmd_top.add_atom_to_residue(new_H,bondingPartner.residue)
            originalAtoms.append(new_H) # add atom to the bottom of "index-stiff" list
    return new_ase_struct, new_pmd_top


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
