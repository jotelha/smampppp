## jlh 2018/02/13
# takes dictionary of united atoms and implicit hydrogen number,
# adding a certain number of explicit hydrogens around these.
# Needs and returns corresponding ASE and ParmEd representations of the system
# Ugly and annoying: 
# ASE identifies atoms only by number, no functionality to store atom name and residue 
#  New atoms are added at the end of the whole atom list
# Parmed inserts new atoms not at the end of the list, nut within their respective residues.
# Initially, ordering in ASE and Parmed are equal, but when structure is modified
# we have to track every change in order to be able to map Parmed structure to ase structure

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

import logging
#gromacs.GROMACS_TOPDIR = "/home/jotelha/gromacs/2016.4/share/gromacs/top"

# returns new_ase_struct and new_mpd_struct with hydrogens added
# additionally, two lists, "names" and "residues" are returned in order to make ASE atoms identifiable
def insertHbyList(ase_struct,pmd_top,implicitHbondingPartners,bond_length=1.0,debug=False):
    # make copies of passed structures as not to alter originals:
    new_pmd_top = pmd_top.copy(pmd.Structure)
    new_ase_struct = ase_struct.copy()

    # names stores the String IDs of all atoms as to facilitate later ASE ID -> Atom name mapping
    names = [ a.name for a in pmd_top.atoms ]
    residues = [ a.residue.name for a in pmd_top.atoms ]

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
        logging.info('Adding {} H-atoms to {} (#{})...'.format(Hno,originalAtoms[k].name,k))
        for h in range(0,Hno):
            r = new_ase_struct.positions
            bondingPartners = j[i==k]
            logging.info('bondingPartners {}'.format(bondingPartners))
            partnerStr=''
            for p in bondingPartners:
                if partnerStr == '':
                    partnerStr = originalAtoms[p].name
                else:
                    partnerStr += ', ' + originalAtoms[p].name
            logging.info('Atom {} already has bonding partners {}'.format(originalAtoms[k].name,partnerStr))
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
                    logging.info('programm needs more than 15 corrector steps for H atom {} at atom {}'.format(n_atoms, k))
                    sys.exit(15)
                    break

                logging.info('too close atoms {}'.format(indices))
                c_step += 1
                logging.info('correcter step {} for H {} at atom {}'.format(c_step, n_atoms-1, k))
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
            logging.info('Adding H-atom {} at position [ {}, {}, {} ]'.format(nameH,r0[0], r0[1], r0[2]))
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
            
            names.append(nameH) # append name of H-atom
            residues.append(bondingPartner.residue.name) # H is in same residue as bonding partner
    return new_ase_struct, new_pmd_top, names, residues