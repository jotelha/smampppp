#!/usr/bin/env python
""" Fits (united-atom) point charges onto (all-atom) ESP obtained by 
    GPAW and HORTON under certain charge group and symmetry constraints 
    as required by GROMOS force fields """
# script which should read in hortons old cost function and add new constraints

# PLAN:
# Write functions which read input files for constraints and use them to
# construct the new cost function. DONE
# After reading in and construction of the
# cost function, solve it and return the fitted charges. DONE

# ENHANCEMENT:
# Check if the constraints are applied properly or if there occured an error. DONE
# Check if the constraints are meaningfull or if they exclude each other...

# TODO:
# maybe write a function which can visualize the symmetry constraints
# to proof the input. Highlight atoms which are constrained or sth.
# like this. DONE IN TEXT FORM

import warnings
import logging

import numpy as np
import pandas as pd

import ase.io
import parmed as pmd
from parmed import gromacs
from insertHbyList import insertHbyList 



### Function definitions ###

def unconstrainedMinimize(A_matrix, b_vector, C_scalar, debug = False):
    """
    Minimization of the cost function.
    Parameters
    ----------
    A_matrix, b_vector, C_scalar:
        from hortons esp cost function. x A x^T - 2 B x + C
    debug: bool
        True:  print debug information
        False: dont print debug information
        default=False

    Return
    ------
    x:  M+N array
        Vector of M optimized charges and N lagrange multipliers
    A: numpy.ndarray
        matrix of horton minimization and lagrange condition
    B: numpy.array
        array of horton and lagrange conditions
    """

    N = b_vector.shape[0]
    #M = q_vector.shape[0]

    npv = float(np.version.version[2:])

    if debug:
        logging.info("A {}: \n {}".format(A_matrix.shape,A_matrix))
        logging.info("B {}: \n {}".format(b_vector.shape,b_vector))
        logging.info("C {}: \n {}".format(C_scalar.shape,C_scalar))
 
    A = A_matrix
    B = b_vector
    C = C_scalar

    x = np.linalg.solve(A, B)

    return x, A, B

def constrainedMinimize(A_matrix, b_vector, C_scalar, D_matrix = None,
                        q_vector = np.array([0]), debug = False):
    """
    Minimization of the cost function.
    Parameters
    ----------
    A_matrix, b_vector, C_scalar:
        from hortons esp cost function. x A x^T - 2 B x + C
    D_matrix:
        matrix of additional constraints
    q_vector:
        vector of constants which should be fulfilled by the additional
        constraints.
    debug: bool
        True:  print debug information
        False: dont print debug information
        default=False

    Return
    ------
    x:  M+N array
        Vector of M optimized charges and N lagrange multipliers
    A: numpy.ndarray
        matrix of horton minimization and lagrange condition
    B: numpy.array
        array of horton and lagrange conditions

    TODO
    ----
        --> check the Lagrange multipliers for consistency, e.g. too large
            Lagrange multipliers can be a hint for not fulfilled constraints?!
        --> for old numpy versions, version < 1.13.0, there is no np.block
            Instead use np.bmat in the same way...
    """
    
    logging.debug("")
    logging.debug("### constrainedMinimize ###")
    logging.debug("")

    N = b_vector.shape[0]
    M = q_vector.shape[0]

    npv = float(np.version.version[2:])

    if not isinstance(D_matrix,np.ndarray):
        D_matrix = np.atleast_2d( np.ones(N) )
        
    #if debug:
    logging.debug("{:d} unknowns, {:d} equality constraints".format(N,M))
    logging.debug("A {}: \n {}".format(A_matrix.shape,A_matrix))
    logging.debug("B {}: \n {}".format(b_vector.shape,b_vector))
    logging.debug("C {}: \n {}".format(C_scalar.shape,C_scalar))
    logging.debug("D {}: \n {}".format(D_matrix.shape,D_matrix))
    logging.debug("q {}: \n {}".format(q_vector.shape,q_vector))
    if npv < 13:
        logging.info('\nWARNING:\n Your numpy version {} is old,\n I am falling from '
            'np.block() back to np.bmat()\n\n'.format(np.version.version) )
        A = np.bmat([[ 2*np.atleast_2d(A_matrix), np.atleast_2d(D_matrix).T ],
                     [ np.atleast_2d(D_matrix), np.atleast_2d(np.zeros((M,M)))]])

    else:
        #maybe combine these steps to one np.block operation like for bmat
        A_upper = np.block(
            [ 2*np.atleast_2d(A_matrix), np.atleast_2d(D_matrix).T ])
        A_lower = np.block(
            [ np.atleast_2d(D_matrix), np.atleast_2d(np.zeros((M,M)))])
        A = np.block( [ [ A_upper ], [ A_lower ] ] )

    logging.debug("block A ({}): \n {}".format(A.shape,A))

    if npv < 13:
        B = np.bmat( [2*np.atleast_1d(b_vector), np.atleast_1d(q_vector)] ).T

    else:
        B = np.block( [2*np.atleast_1d(b_vector), np.atleast_1d(q_vector)] )

    logging.debug("block B ({}): \n {}".format(B.shape,B))

    C = C_scalar
    
    rank_A = np.linalg.matrix_rank(A)
    rank_AB = np.linalg.matrix_rank(np.hstack((A,np.atleast_2d(B).T))) 
    if rank_A == rank_AB:
        logging.info("A rank == [A,B] rank == {:d}, solvable".format(rank_A))
        if rank_A == A.shape[0]:
            logging.info("A rank == A dim, unique solution exists.")
        else:
            logging.info("A rank < A dim {:d}, no unique solution.".format(A.shape[0]))
    else:
        logging.info("A rank {:d} != [A,B] rank {:d}, unsolvable".format(
            rank_A, rank_AB))

    # What does numpy do for ambiguous systems?
    x = np.linalg.solve(A, B)

    return x, A, B


def constructPairwiseSymmetryConstraints(charges, N, symmetry=1.0, debug=False):
    """
    Function to construct D_matrix and q_vector for a pairwise symmetry

    Parameters:
    -----------
    charges: (1D or 2D numpy.ndarray int) 
            or (list of int) or (list of list of int)
        List of ASE charge indices which should be equal. 
        For more than two charges
        there will be (N-1) constraints for pairwise equal charges. If you give
        a list of lists with equal charges each sublist will be forced to have
        equal charges.
    N: int
        Total number of atoms in your system. Needed to fix the size of the
        D_matrix and q_vector to the system size (number of atoms)
    symmetry: int
         1: the pairwise charges are equal (symmetric), q_1 = q_2
        -1: the pairwise charges are the negative of each other (antisymmetric),
            q_1 = -q_2
        default=1
    debug: bool
        True:  print debug information
        False: dont print debug information
        default=False

    Return
    ------
    D: np.ndarray, dim=2
        D_matrix which carries the constraints of pairwise symmetric charges
    q: np.array
        q_vector carying the total charge of q_1+q_2 (always zero)

    TODO
    ----
        implement 2D charge array such that one can input at once all pairwise
        symmetries with an symmetry array. DONE
    """

   
    logging.debug("")
    logging.debug("### constructPairwiseSymmetryConstraints ###")
    logging.debug("")
        
    #charges = np.atleast_2d(charges)
    #D = np.ones((1, N))
    #q = np.ones((1))
    D = []
    q = []
    
    # loop over sets of equally charged atoms
    for charge_list in charges:
        M = len(charge_list)-1 # number of pairwise constraints for the current list

        #symmetry = symmetry*np.ones(M)

        logging.debug("charge list ({}): {}".format(charge_list.shape,charge_list))
        logging.debug("{:d} unknowns, {:d} pairwise equality constraints".format(N,M))
        logging.debug("symmetry: {}".format(symmetry))

        D_single = np.atleast_2d(np.zeros((M,N)))
        q_single = np.atleast_1d(np.zeros(M))
        D_single[:,charge_list[0]] = 1 # 1st atom in list used in every constraint

        for j in range(M):
            D_single[j,charge_list[j+1]] = -1.0*symmetry

        logging.debug("D_single ({}):\n{}".format(D_single.shape,D_single))
        logging.debug("q_single ({}): {}".format(q_single.shape,q_single))

        #add up D_single and q_single to D and q
        D.append(D_single)
        q.append(q_single)
        
    D,q = concatenated_constraints(D,q)
    return D, q

### old constructPairwiseSymmetryConstraints ###
#def constructPairwiseSymmetryConstraints(charges, N, symmetry=1.0, debug=False):
#    """
#    Function to construct D_matrix and q_vector for a pairwise symmetry
#
#    Parameters:
#    -----------
#    charges: list of ints
#        List of charge indices which should be equal. For more than two charges
#        there will be (N-1) constraints for pairwise equal charges.
#    N: int
#        Total number of atoms in your system. Needed to fix the size of the
#        D_matrix and q_vector to the system size (number of atoms)
#    symmetry: int, -1 or 1
#         1: the pairwise charges are equal (symmetric), q_1 = q_2
#        -1: the pairwise charges are the negative of each other (antisymmetric),
#            q_1 = -q_2
#        default=1
#    debug: bool
#        True:  print debug information
#        False: dont print debug information
#        default=False
#
#    Return
#    ------
#    D: np.ndarray, dim=2
#        D_matrix which carries the constraints of pairwise symmetric charges
#    q: np.array
#        q_vector carying the total charge of q_1+q_2 (always zero)
#
#    TODO
#    ----
#        implement 2D charge array such that one can input at once all pairwise
#        symmetries with an symmetry array.
#    """
#    M = len(charges)-1
#
#    symmetry = symmetry*np.ones(M)
#
#    if debug:
#        logging.info("{:d} unknowns, {:d} pairwise equality constraints".format(N,M))
#        logging.info("symmetry list ({}):\n{}".format(symmetry.shape,symmetry))
#
#    D = np.atleast_2d(np.zeros((M,N)))
#    q = np.atleast_1d(np.zeros(M))
#    D[:,charges[0]] = 1
#
#    for j in range(M):
#        D[j,charges[j+1]] = -1.0*symmetry[j]
#
#    if debug:
#        logging.info("D ({}):\n{}".format(D.shape,D))
#        logging.info("q ({}):\n{}".format(q.shape,q))
#
#    return D, q


def constructChargegroupConstraints(chargeGroups, N, q=0, debug=False):
    """
    Function to construct D_matrix and q_vector for charge Groups

    Parameters
    ----------
    chargeGroups: list, or 2-D list; ints
        list of atom indices which belong to one charge group with charge q. For
        more than one charge group you can use a two dimensional list (list of
        charge groups) and use a list q, for the charges of the charge groups.
    N: int
        Total number of atoms in your system. Needed to fix the size of the
        D_matrix and q_vector to the system size (number of atoms)
    q: scalar, or list; reals
        describes the total charge of a charge group. Therefore it is a scalar
        for one charge group and a list if more than one charge group is given
    debug: bool
        True:  print debug information
        False: dont print debug information
        default=False

    Returns
    -------
    D_matrix: np.ndarray, dim=2
        D_matrix which carries the constraints of pairwise symmetric charges
    q_vector: np.array
        q_vector carying the total charge of q_1+q_2 (always zero)
    """

    logging.debug("")
    logging.debug("### constructChargegroupConstraints ###")
    logging.debug("")

    M = len(chargeGroups)

    q_vector = np.atleast_1d(q*np.ones(M))

    logging.debug("{:d} unknowns, {:d} pairwise equality constraints".format(N,M))

    D_matrix = np.atleast_2d(np.zeros((M,N)))
    #q = np.atleast_2d(np.zeros(M))

    for j in range(M):
        D_matrix[j,chargeGroups[j]] = 1.0

    logging.debug("D_matrix ({}):\n{}".format(D_matrix.shape, D_matrix))
    logging.debug("q_vector ({}):\n{}".format(q_vector.shape, q_vector))

    return D_matrix, q_vector

def constructTotalChargeConstraint(charge, N, symmetry=1.0, debug=False):
    """
    Function to construct D_matrix and q_vector for the total charge constraint

    Parameters:
    -----------
    charges: float
        Total required charge of system.
    N: int
        Total number of atoms in your system. Needed to fix the size of the
        D_matrix and q_vector to the system size (number of atoms)
    debug: bool
        True:  print debug information
        False: dont print debug information
        default=False

    Return
    ------
    D: np.ndarray, dim=2
        D_matrix filled with ones
    q: np.array
        q_vector carying the (scalar) total charge

    """

    D = np.atleast_2d(np.ones((1, N)))
    q = np.atleast_1d(charge)

    return D,q

# reduce D to matrix of full rank (and accordingly adjust q)
# ATTENTION: contradictory constraints might vanish unnoticed
def construct_D_of_full_rank(D,q):
    D_LI=[D[0]]
    q_LI=[q[0]]
    for i in range(D.shape[0]):
        tmp=[]
        for r in D_LI:
            tmp.append(r)
        tmp.append(D[i]) 
        if np.linalg.matrix_rank(tmp)>len(D_LI):    
            #test if row D[i] is linearly independent from all (row) vectors in D_LI
            D_LI.append(D[i]) 
            #note that matrix_rank does not need to take in a square matrix
            q_LI.append(q[i])
    return np.array(D_LI), np.array(q_LI) 

def read_AtomName_ChargeGroup(file_name, ase2pmd):
    """
    Function to read in csv file of atom names and associated charge groups with 
    the help of pandas. Expects dictionary of ASE indices pointing at atom and 
    residue name tyuples and constructs vector of associated charge
    group. The vectors are all ordered in the same manner. It also returns
    ncgtypes, to reconstruct the charge group if atom names occure more than ones.

    Parameters
    ----------
    file_name: str
        name of the file which is read. File format: one line per atom,
        str, int
        [atom name],[charge group id]
    ase2pmd: dict
        dictionary of ASE atom indices mapped onto tuples of
        ParmEd names and residues

    Return
    ------
    cg2ase: list of list of int
        each sublist represents one charge group and contains ASE indices of 
        all atoms in the group. 
    cg2cgtype: list of int
        Charge groups can reoccur across several residues. 
        Thus every charge group is given an (arbitrary) unique index internally.
        This list maps internal cg indices onto their original (ambiguous) id.
    ncgtypes: int
        the  number of (original) charge groups which are might be used for 
        more than one atom index and thus carry an ambiguous meaning.
    """

    ase2pmd_df = pd.DataFrame(ase2pmd).T
    ase2pmd_df.columns = ['atom','residue']    
    
    # column 0: atom names
    # column 1: charge group number
    pmd2cg_df = pd.read_csv(file_name, sep=',', header=None,  
                            comment='#', names=['atom','cg'])   

    unique_residues = ase2pmd_df['residue'].unique()
    unique_charge_groups = pmd2cg_df['cg'].unique()
        
    ncgtypes = len(unique_charge_groups)
    cg2ase = []
    cg2cgtype = []                                                
    for r in unique_residues:
        # atoms_in_residue = np.where(residues == r)
        for cgtype in unique_charge_groups:
            #np.where( (residues == r ) && ( names == name_cg[:,1] == cg )
            names_in_cg = pmd2cg_df[pmd2cg_df['cg'] == cgtype]['atom']
            cg_sel = ase2pmd_df['atom'].isin(names_in_cg)
            res_sel = (ase2pmd_df['residue'] == r )
            
            new_cg = ase2pmd_df[cg_sel & res_sel]
                                                
            if not new_cg.empty:                                
                cg2ase.append(new_cg.index.values)
                cg2cgtype.append(cgtype)                              
                     
    return cg2ase, cg2cgtype, ncgtypes


def read_ChargeGroup_TotalCharge(file_name):
    """
    Function to read in csv file of charge groups and charges with the help of
    pandas.

    Parameters
    ----------
    file_name: str
        name of the file which is read. File format: 
        int, float
        [ charge group number ], [ charge ]

    Return
    ------
    cg_q: dict of int: float
        dictionary of charge group id ('type') and corresponding total charge.
    """

    cg_q_df = pd.read_csv(file_name, sep=',', header=None, comment='#', 
                          index_col=0, names=['charge'])
    cg_q = cg_q_df['charge'].to_dict()

    return cg_q


# ATTENTION: constructs symmetry constraints across ALL residues:
def read_SameChargedAtoms(file_name, ase2pmd):
    """
    Function to read in csv file of atoms which should have the same charge.
    Automatically enforces same charges for atoms of same name, but
    possibly spread over different residues.

    Parameters
    ----------
    file_name: str
        name of the file which is read. Format: pairwise equalities
        str, str
        [ atom name 1 ], [ atom name 2]
        
    name: list 
        atom names, as ordered in ASE

    Return
    ------
    sym2ase: list of list of int
        each sublist groups atoms of same charge by their ASE indices.
    """
    logging.debug("")
    logging.debug("### read_SameChargedAtoms ###")
    
    ase2pmd_df = pd.DataFrame(ase2pmd).T
    ase2pmd_df.columns = ['atom','residue']    
    
    sca_df = pd.read_csv(file_name, sep=',', header=None,  comment='#')

    sym2ase = []

    # constructs symmetries on atoms of same type across residues
    # this will result in redundant constraints for symmetries 
    # specidfied explicitly in input file for pairwise symmetries
    # for now, we do not care about that
    # TODO: reduce all (possibly recdundant) constraints into 
    # a minimal set of unique ones
    
    logging.debug("")
    logging.debug("Constructing implict symmetries due to atom types.")
    unique_atoms = ase2pmd_df['atom'].unique()
    for a in unique_atoms:
        new_symmetry_group = ase2pmd_df[ ase2pmd_df['atom'] == a ]
        if not new_symmetry_group.empty:                                
            if len(new_symmetry_group.index.values) < 2:
                logging.debug("Apparently, atom type {} only occurs once "
                             "at ASE index {:d}. No symmetry constraint "
                             "created.".format(a,new_symmetry_group.index.values[0]))
            else:
                logging.debug("Add symmetry constraint for atom type {} "
                    "at ASE indices {}.".format(a,new_symmetry_group.index.values))
                sym2ase.append(new_symmetry_group.index.values)   
                
    logging.debug("")
    logging.debug("Constructing explicit symmetries.")
    
    # for i, group in sca_df.iterrows():
    #    sca_sel = []
    #    new_symmetry_group = []
    #    for atom in group:
    #        # Only selects first occurence of atom type on system
    #        # and appends it to symmetry group selection, as all other 
    #        # symmetries accross residues are enforced by implicity 
    #        # name-wise constraints constructed above
    #        sca_sel = ase2pmd_df[ase2pmd_df['atom'] == atom]
    #        if not sca_sel.empty:
    #            new_symmetry_group.append( sca_sel.index[0] )

    #    if new_symmetry_group:
    #        logging.info("Add explicit symmetry constraint for types {} "
    #                "at ASE indices {}.".format(group.values,
    #                            new_symmetry_group ))
    #        sym2ase.append( np.array(new_symmetry_group) )  
    
    for i, group in sca_df.iterrows():
        sca_sel = ase2pmd_df['atom'].isin(group)
        new_symmetry_group = ase2pmd_df[sca_sel]

        if not new_symmetry_group.empty:                                
            logging.debug("Add explicit symmetry constraint for types {} "
                    "at ASE indices {}.".format(group.values,
                        new_symmetry_group.index.values))
            sym2ase.append(new_symmetry_group.index.values)
            
    return sym2ase


def concatenated_constraints(D_matrices, q_vectors):
    """
    Function to concatenate D_matrices and q_vectors to one D_matrix and one
    q_vector by using numpy.hstack and numpy.vstack. The order of D_matrix in
    D_matrices and q_vector in q_vectors should be the same, otherwise the
    constraints are connected wrong.

    Parameters
    ----------
    D_matrices: list of numpy.ndarray (reals)
        list of all D_matrices which should be concateneted
    q_vectors: list of numpy.ndarray (reals)
        list of all q_vectors which should be concateneted

    Return
    ------
    D_matrix: np.ndarray, dim=2
        D_matrix which carries all constraints of the input D_matrices
    q_vector: np.array
        q_vector which carries all constraints of the input q_vectors
    """

    D_matrix = D_matrices[0]
    for d in D_matrices[1:]:
        D_matrix = np.vstack([D_matrix, d])

    q_vector = q_vectors[0]
    for q in q_vectors[1:]:
        q_vector = np.hstack([q_vector, q])

    return D_matrix, q_vector


def read_horton_cost_function(file_name, debug=False):
    """
    Function to read in hortons cost function. You need h5py to read it.
    We read out the three variables A, B, C which characterise the cost function
    by: X^T A X - 2 B X + C, which is the function to minimize.
    Parameters
    ----------
    file_name: str
        file name of the cost function writen by Horton, typically something like
        'xyz.cost.h5' or 'xyz_cost.h5'.

    Return
    ------
    A_horton: 2D numpy.ndarray

    B_horton: numpy.array

    C_horton: float

    N_horton: int
        N_horton is the number of atoms of the structure
    """

    import h5py

    cost_function = h5py.File(file_name)
    cost = cost_function['cost']
    A_horton = cost['A'][:]
    B_horton = cost['B'][:]
    C_horton = cost['C'].value
    N_horton = cost['natom'].value
    if debug:
        logging.info("A: {}".format(A_horton))
        logging.info("B: {}".format(B_horton))
        logging.info("C: {}".format(C_horton))
        logging.info("N: {}".format(N_horton))

    return A_horton, B_horton, C_horton, N_horton


def logResults(X,A,B,C,N):
    """
    Function to log results for debugging purposes.
    
    Parameters
    ----------
    X: np.ndarray, dim=2
        Optimized results as yielded by constrainedMinimize(...)

    A: 2D numpy.ndarray
    B: numpy.array
    C: float
        describe the cost function used for fitting
        
    N: int
        N is the number of atoms of the structure
        
    Return
    ------
    None
    """    
    
    np.set_printoptions(precision=3)
    np.set_printoptions(suppress=True)

    logging.info('charges {}:\n {}\ncharge sum = {}\n'.format( X[:N].T.shape,
                                                        X[:N].T,
                                                        X[:N].T.sum() ))
    logging.info('Lagrange multipliers {}:\n {}'.format( X[N:].T.shape,
                                                  X[N:].T ) )

    ### test the results
    logging.info( 'value of cost function: {}'.format(
        (np.dot(X.T, np.dot(A, X)) - 2*np.dot(B.T, X) - C) ) )
    
# check charge group constraints:
def checkChargeGroups( df, cg2ase, cg2cgtype, cg2q,
    q_cols = ['q','q_unconstrained','q_qtot_constrained',
              'q_cg_qtot_constrained', 'q_sym_qtot_constrained']):
    
    logging.info("")
    logging.info("")
    logging.info("##############################")    
    logging.info("CHARGE GROUP CONSTRAINTS CHECK")            
    logging.info("##############################")    
    logging.info("")
    logging.info("atoms grouped together by their ASE indices:")
    logging.info("{}".format(cg2ase))    
    logging.info("")
    logging.info("desired charge of each group:")
    logging.info("{}".format(cg2q))
    
    for cg_index, ase_indices_in_cg in enumerate(cg2ase):
        logging.info("cg {:d}, type {:d}:".format(cg_index,cg2cgtype[cg_index]))
        for q_col in q_cols:
            q_cg = df.iloc[ase_indices_in_cg][q_col].sum() # select first charge group
            logging.info(
                "    {:>30}:{:8.4f}    absolute error:{:12.4e}".format(
                    q_col,q_cg,q_cg-cg2q[cg_index]))
            
# check symmetry constraints:
def checkSymmetries( df, sym2ase, 
    q_cols = ['q','q_unconstrained','q_qtot_constrained',
              'q_cg_qtot_constrained', 'q_sym_qtot_constrained']):
    
    logging.info("")    
    logging.info("")
    logging.info("##########################")    
    logging.info("SYMMETRY CONSTRAINTS CHECK")            
    logging.info("##########################")    
    logging.info("")
    logging.info("groups of equally charged atoms by their ASE indices:")
    logging.info("{}".format(sym2ase))
    
    for sym_index, ase_indices_in_sym in enumerate(sym2ase):
        #logging.info("cg {:d}, type {:d}:".format(cg_index,cg2cgtype[cg_index]))
        msg = []
        for ase_index in ase_indices_in_sym:
                msg.append("({}, {})".format(
                    df.iloc[ase_index]['atom'], 
                    df.iloc[ase_index]['residue']))
                           
        logging.info("sym {:d}: {}".format(sym_index,"; ".join(msg)))
              
        for q_col in q_cols:
            msg = []
            for ase_index in ase_indices_in_sym:
                msg.append("{:.3f}".format(df.iloc[ase_index][q_col]))
            logging.info("{:>30}: {}".format(q_col,",".join(msg)))    
            
        logging.info("")


def fitESPconstrained(infile_pdb, infile_top, infile_cost_h5, 
    infile_atoms_in_cg_csv, infile_cg_charges_csv, 
    infile_atoms_of_same_charge_csv,
    qtot = 0.0, strip_string=':SOL,CL', 
    implicitHbondingPartners = {'CD4':1,'CD3':1,'CA2':2,'CA3':2,'CB2':2,'CB3':2},
    debug=False, outfile_top = None):
    
    """
    Automizes the whole fitting process from importing Horton's
    cost function over reading constraints from simple text files to
    minimizing, logging and double-checking the results.
    
    Parameters
    ----------
    infile_pdb: str
        PDB file with original (united-atom) molecular structure 
    infile_top: str
        GROMACS topolgy file with original (united-atom) system.
        All #includes shoulb be removed!
    infile_cost_h5: str
        Cost function by HORTON, hdf5 format
    infile_atoms_in_cg_csv: str
        file with atom - charge group assignments in simple 
        "comma separated value" text format, one line per atom:
        str, int
        [atom name],[charge group id]     
    infile_cg_charges_csv: str
        file with charge group - charge assignments in simple 
        "comma separated value" text format, one line per charge group:
        int, float
        [ charge group number ], [ charge (in e) ]
    infile_atoms_of_same_charge_csv: str
        file with pairwise atom symmetry assignments in simple 
        "comma separated value" text format, one line per equality:
        str, str
        [ atom name 1 ], [ atom name 2]
        will have the same charge. Apart from that, all atoms of the same 
        name (but possibly spread over different residues) will have the 
        same charge enforced.
    qtot: float
        The system's total charge
    strip_string: str
        Groups to remove from the initally imported topology in ParmEd.
        ':SOL,CL' by default (solvent and chlorine ions).
    implicitHbondingPartners: dict of str: int
        By default "{'CD4':1,'CD3':1,'CA2':2,'CA3':2,'CB2':2,'CB3':2}"
        Specifies which atoms have (how many) implicit hydrogens around them.
        These hydrogens must equal those used in QM calculations.
    debug: bool
        By default False, uses logging if True
        
        
    Return
    ------
    q: np.ndarray of float, dim=1
        fitted charges, fully constrained
    lambda: np.ndarray of float, dim=1
        Lagrange multipliers
    info_df: pandas.DataFrame
        containing information on the fit in easily accesible pandas dataframe
    cg2ase: list of list of int 
    cg2cgtype: list of int
    cg2q: list of float
    sym2ase: list of list of int
    """    
    
    logging.info("#################")    
    logging.info("fitESPconstrained")            
    logging.info("#################")   

    # A: construct all-atom representation from united-atom structure and topology:
    ua_ase_struct = ase.io.read(infile_pdb)
    ua_pmd_struct = pmd.load_file(infile_pdb)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ua_pmd_top = gromacs.GromacsTopologyFile(infile_top,parametrize=False)
    # throws some warnings on angle types, does not matter for bonding info
    # if error thrown, just try to "reduce" .top as far as possible
    # warnings supressed as shown on
    # https://docs.python.org/2/library/warnings.html

    ua_pmd_top.strip(strip_string) 
        #strip water and electrolyte from system (if not yet done in .top)
    ua_pmd_top.box = ua_pmd_struct.box # Needed because .pdb contains box info
    ua_pmd_top.positions = ua_pmd_struct.positions

    ua_names = [ a.name for a in ua_pmd_top.atoms ]
    ua_residues = [ a.residue.name for a in ua_pmd_top.atoms ]

    aa_ase_struct, aa_pmd_struct, aa_names, aa_residues = \
        insertHbyList(ua_ase_struct,ua_pmd_top,
        implicitHbondingPartners,1.0)

    ua_count = len(ua_ase_struct)     # united atoms structure
    aa_count = len(aa_ase_struct) # all atoms structure

    ua_ase_index = np.arange(ua_count)
    aa_ase_index = np.arange(aa_count)

    aa_atom_residue_list = list(zip(aa_names,aa_residues))
    aa_ase_index = range(aa_count)
    aa_ase2pmd = dict(zip(aa_ase_index,aa_atom_residue_list))
    aa_pmd2ase = dict(zip(aa_atom_residue_list,aa_ase_index))

    ua_atom_residue_list = list(zip(ua_names,ua_residues))
    ua_ase_index = range(ua_count)
    ua_ase2pmd = dict(zip(ua_ase_index,ua_atom_residue_list))
    ua_pmd2ase = dict(zip(ua_atom_residue_list,ua_ase_index))
    
    # TODO: distinction for ua and aa fitting:
    pmd_struct = ua_pmd_struct
    pmd_top = ua_pmd_top
    ase2pmd = ua_ase2pmd
    pmd2ase = ua_pmd2ase
    
    # B: read cost function
    
    A_horton, B_horton, C_horton, N_horton = \
        read_horton_cost_function(file_name = infile_cost_h5)
        
    # C: read constraints files
    
    ### Charge Groups:
    # read in all charge groups and construct the corresponding constraints
    cg2ase, cg2cgtype, ncgtypes = read_AtomName_ChargeGroup(
        file_name = infile_atoms_in_cg_csv, ase2pmd = ase2pmd)
    
    cg_q = read_ChargeGroup_TotalCharge(file_name = infile_cg_charges_csv)
    
    cg2q  = [ cg_q[cg] for cg in cg2cgtype ]
    
    ### Same Charged Atoms
    sym2ase = read_SameChargedAtoms(
        file_name = infile_atoms_of_same_charge_csv, ase2pmd = ase2pmd)
    
    # D: construct constraints matrices
    D_matrix_cg_red, q_vector_cg_red = constructChargegroupConstraints(
        chargeGroups = cg2ase, N = N_horton, q = cg2q, debug = debug)
    
    D_matrix_sym_red, q_vector_sym_red = constructPairwiseSymmetryConstraints(
        charges = sym2ase, N = N_horton, symmetry = 1.0, debug = False)
    
    D_matrix_qtot, q_vector_qtot = constructTotalChargeConstraint(
        charge = qtot, N = N_horton)
    
    D_matrix_all_red, q_vector_all_red = concatenated_constraints(
        D_matrices = [D_matrix_cg_red,D_matrix_sym_red,D_matrix_qtot],
        q_vectors = [q_vector_cg_red,q_vector_sym_red,q_vector_qtot])    
    
    # remove redundant constraints 
    D_matrix_cg, q_vector_cg = construct_D_of_full_rank(
        D_matrix_cg_red, q_vector_cg_red)
    D_matrix_sym, q_vector_sym  = construct_D_of_full_rank(
        D_matrix_sym_red, q_vector_sym_red )
    D_matrix_all, q_vector_all = construct_D_of_full_rank(
        D_matrix_all_red, q_vector_all_red)
    
    if debug:
        D_matrix_cg_qtot_red, q_vector_cg_qtot_red = concatenated_constraints(
            D_matrices = [D_matrix_cg_red,D_matrix_qtot],
            q_vectors = [q_vector_cg_red,q_vector_qtot])        
        D_matrix_cg_qtot, q_vector_cg_qtot = construct_D_of_full_rank(
            D_matrix_cg_qtot_red, q_vector_cg_qtot_red)
        
        
        D_matrix_sym_qtot_red, q_vector_sym_qtot_red = concatenated_constraints(
            D_matrices = [D_matrix_sym_red,D_matrix_qtot],
            q_vectors = [q_vector_sym_red,q_vector_qtot])    
        D_matrix_sym_qtot, q_vector_sym_qtot = construct_D_of_full_rank(
            D_matrix_sym_qtot_red, q_vector_sym_qtot_red)
                
        logging.info("")
        logging.info("")
        logging.info("")
        logging.info("################")    
        logging.info("CONSTRAINTS INFO")            
        logging.info("################")    
        logging.info("")
        logging.info("")
        
        # CG CONSTRAINTS
        rank_D_cg = np.linalg.matrix_rank(D_matrix_cg)
        rank_Dq_cg = np.linalg.matrix_rank(np.hstack((D_matrix_cg,
                                           np.atleast_2d(q_vector_cg).T)))
        logging.info("{:d} CG constraints, rank {:d} "
                     " => {:d} redundant or contradictory constraints.".format(
                         D_matrix_cg_red.shape[0], rank_D_cg, 
                         D_matrix_cg_red.shape[0] - rank_D_cg))
        if rank_D_cg == rank_Dq_cg:
            logging.info("    D rank == [D,q] rank == {:d},"
                         " solvable".format(rank_D_cg))
        else:
            logging.info("    D rank {:d} != [D,q] rank {:d},"
                         " unsolvable".format(rank_D_cg, rank_Dq_cg))

        # CG + QTOT CONSTRAINTS
        rank_D_cg_qtot = np.linalg.matrix_rank(D_matrix_cg_qtot)
        rank_Dq_cg_qtot = np.linalg.matrix_rank(np.hstack((D_matrix_cg_qtot,
                                           np.atleast_2d(q_vector_cg_qtot).T)))
        logging.info("{:d} CG and QTOT constraints, rank {:d} "
                     " => {:d} redundant or contradictory constraints.".format(
                         D_matrix_cg_qtot_red.shape[0], rank_D_cg_qtot, 
                         D_matrix_cg_qtot_red.shape[0] - rank_D_cg_qtot))
        if rank_D_cg_qtot == rank_Dq_cg_qtot:
            logging.info("    D rank == [D,q] rank == {:d},"
                         " solvable".format(rank_D_cg_qtot))
        else:
            logging.info("    D rank {:d} != [D,q] rank {:d},"
                         " unsolvable".format(rank_D_cg_qtot, rank_Dq_cg_qtot))

            
        # SYM CONSTRAINTS
        rank_D_sym = np.linalg.matrix_rank(D_matrix_sym)
        rank_Dq_sym = np.linalg.matrix_rank(np.hstack((D_matrix_sym,
                                           np.atleast_2d(q_vector_sym).T)))
        logging.info("{:d} SYM constraints, rank {:d} "
                     " => {:d} redundant or contradictory constraints.".format(
                         D_matrix_sym_red.shape[0], rank_D_sym, 
                         D_matrix_sym_red.shape[0] - rank_D_sym))
        if rank_D_sym == rank_Dq_sym:
            logging.info("    D rank == [D,q] rank == {:d},"
                         " solvable".format(rank_D_sym))
        else:
            logging.info("    D rank {:d} != [D,q] rank {:d}," 
                         " unsolvable".format(rank_D_sym, rank_Dq_sym))

        # SYM + QTOT CONSTRAINTS
        rank_D_sym_qtot = np.linalg.matrix_rank(D_matrix_sym_qtot)
        rank_Dq_sym_qtot = np.linalg.matrix_rank(np.hstack((D_matrix_sym_qtot,
                                           np.atleast_2d(q_vector_sym_qtot).T)))
        logging.info("{:d} SYM and QTOT constraints, rank {:d} "
                     " => {:d} redundant or contradictory constraints.".format(
                         D_matrix_sym_qtot_red.shape[0], rank_D_sym_qtot, 
                         D_matrix_sym_qtot_red.shape[0] - rank_D_sym_qtot))
        if rank_D_sym_qtot == rank_Dq_sym_qtot:
            logging.info("    D rank == [D,q] rank == {:d},"
                         " solvable".format(rank_D_sym_qtot))
        else:
            logging.info("    D rank {:d} != [D,q] rank {:d},"
                         " unsolvable".format(rank_D_sym_qtot, rank_Dq_sym_qtot))
            
        # ALL CONSTRAINTS
        rank_D_all = np.linalg.matrix_rank(D_matrix_all)
        rank_Dq_all = np.linalg.matrix_rank(np.hstack((D_matrix_all,
                                           np.atleast_2d(q_vector_all).T)))
        logging.info("{:d} ALL constraints, rank {:d} "
                     " => {:d} redundant or contradictory constraint.".format(
                         D_matrix_all_red.shape[0], rank_D_all, 
                         D_matrix_all_red.shape[0] - rank_D_all))
        if rank_D_sym == rank_Dq_sym:
            logging.info("    D rank == [D,q] rank == {:d},"
                         " solvable".format(rank_D_all))
        else:
            logging.info("    D rank {:d} != [D,q] rank {:d},"
                         " unsolvable".format(rank_D_all, rank_Dq_all))
                 
    # E: Minimization 
    
    ### Constrained minimization
    logging.info("")    
    logging.info("########################")    
    logging.info("FULLY CONSTRAINED SYSTEM")            
    logging.info("########################")
    logging.info("")  
    X, A, B = constrainedMinimize(A_matrix = A_horton,
                        b_vector = B_horton,
                        C_scalar = C_horton,
                        D_matrix = D_matrix_all,
                        q_vector = q_vector_all,
                        debug    = debug)
    
    ase2pmd_df = pd.DataFrame(ase2pmd).T
    ase2pmd_df.columns = ['atom','residue']
    ase2pmd_df['q'] = X[:N_horton]

    # additional debug cases
    if debug:     
        ### Unconstrained minimization
        logging.info("")    
        logging.info("####################")    
        logging.info("UNCONSTRAINED SYSTEM")            
        logging.info("####################")
        logging.info("")  
        X_unconstrained, A_unconstrained, B_unconstrained = \
            unconstrainedMinimize(A_matrix = A_horton,
                        b_vector = B_horton,
                        C_scalar = C_horton,
                        debug    = debug)
        
        ### Total charge constraint minimization
        logging.info("")    
        logging.info("####################################")    
        logging.info("SYSTEM WITH TOTAL CHARGE CONSTRAINED")            
        logging.info("####################################")
        logging.info("")  
        X_qtot_constraint, A_qtot_constraint, B_qtot_constraint = \
            constrainedMinimize(A_matrix = A_horton,
                        b_vector = B_horton,
                        C_scalar = C_horton,
                        D_matrix = D_matrix_qtot,
                        q_vector = q_vector_qtot,
                        debug    = debug)
        
        ### Charge group & total charge constraint minimization 
        logging.info("")    
        logging.info("######################################################")    
        logging.info("SYSTEM WITH TOTAL CHARGE AND CHARGE GROUPS CONSTRAINED")            
        logging.info("######################################################")
        logging.info("")  
        X_cg_qtot, A_cg_qtot, B_cg_qtot = \
            constrainedMinimize(A_matrix = A_horton,
                        b_vector = B_horton,
                        C_scalar = C_horton,
                        D_matrix = D_matrix_cg_qtot,
                        q_vector = q_vector_cg_qtot,
                        debug    = debug)
            
        ### Symmetry & total charge constraint minimization
        logging.info("")    
        logging.info("###################################################")    
        logging.info("SYSTEM WITH TOTAL CHARGE AND SYMMETRIES CONSTRAINED")            
        logging.info("###################################################")
        logging.info("")
        X_sym_qtot, A_sym_qtot, B_sym_qtot = \
            constrainedMinimize(A_matrix = A_horton,
                        b_vector = B_horton,
                        C_scalar = C_horton,
                        D_matrix = D_matrix_sym_qtot,
                        q_vector = q_vector_sym_qtot,
                        debug    = debug)
        
        logging.info("")
        logging.info("")
        logging.info("")
        logging.info("#################################")    
        logging.info("RESULTS FOR DIFFERENT CONSTRAINTS")            
        logging.info("#################################")    
        logging.info("")
        logging.info("")
        logging.info("### UNCONSTRAINED ###")
        logResults(X_unconstrained,A_unconstrained,B_unconstrained,C_horton,N_horton)
        
        logging.info("")    
        logging.info("")
        logging.info("### QTOT CONSTRAINED ###")
        logResults(X_qtot_constraint,A_qtot_constraint,
                    B_qtot_constraint,C_horton,N_horton)
        
        logging.info("")
        logging.info("")
        logging.info("### QTOT & CG CONSTRAINED ###")
        logResults(X_cg_qtot,A_cg_qtot,B_cg_qtot,C_horton,N_horton)
        
        logging.info("")
        logging.info("")
        logging.info("### QTOT & SYM CONSTRAINED ###")
        logResults(X_sym_qtot,A_sym_qtot,B_sym_qtot,C_horton,N_horton)


        logging.info("")
        logging.info("")
        logging.info("### FULLY CONSTRAINED ###")
        logResults(X,A,B,C_horton,N_horton)
        
        #ase2pmd_df.columns.append(['q_unconstrained', 'q_qtot_constrained', 'q_qtot_cg_constrained'])
        ase2pmd_df['q_unconstrained'] = X_unconstrained
        ase2pmd_df['q_qtot_constrained'] = X_qtot_constraint[:N_horton]
        ase2pmd_df['q_cg_qtot_constrained'] = X_cg_qtot[:N_horton]
        ase2pmd_df['q_sym_qtot_constrained'] = X_sym_qtot[:N_horton]
        
        checkChargeGroups(ase2pmd_df,cg2ase,cg2cgtype,cg2q)
        checkSymmetries(ase2pmd_df,sym2ase)

 
    #atom_charge_dict = dict(zip(names[0:ua_count],charges))
    if outfile_top:
        for a in pmd_top.atoms:
            a.charge = X[ pmd2ase[(a.name,a.residue.name)] ]
            #pmd_top.write(outfile_top)
            pmd_top.save(outfile_top, overwrite=True)
            
    logging.info("####")    
    logging.info("DONE")            
    logging.info("####")   
    
    return X[:N_horton], X[N_horton:], ase2pmd_df, cg2ase, cg2cgtype, cg2q, sym2ase

### ACTUAL PROGRAM ###
#--------------------#
def main():
    import sys
    import ast
    import argparse
    

    parser = argparse.ArgumentParser(prog='esp-fit-constrained.py',
        description='Estimate charges from a HORTON ESP cost function'
                                     'under arbitrary constraints.')
    parser.add_argument('infile_cost_h5',metavar='cost.h5',
        help='The location of the HORTON cost function in the form '
             '"file.h5:group/cost". This argument must be the same as the '
             'output argument of the script horton-esp-cost.py.')
    parser.add_argument('infile_pdb', metavar='infile.pdb',
        help='PDB file with original (united-atom) molecular structure.')
    parser.add_argument('infile_top', metavar='infile.top',
        help='GROMACS topolgy file with original (united-atom) system. '
        'All #includes shoulb be removed!')
    parser.add_argument('infile_atoms_in_cg_csv', metavar='infile_atoms_in_cg.csv',
        help='file with atom - charge group assignments in simple ' 
        '"comma separated value" text format, one line per atom: '
        'str, int / [atom name],[charge group id]')
    parser.add_argument('infile_cg_charges_csv', metavar='infile_cg_charges.csv',
        help='file with charge group - charge assignments in simple' 
        '"comma separated value" text format, one line per charge group:'
        'int, float / [ charge group number ], [ charge (in e) ]')
    parser.add_argument('infile_atoms_of_same_charge_csv', 
                        metavar='infile_atoms_of_same_charge.csv',
        help='file with pairwise atom symmetry assignments in simple' 
        '"comma separated value" text format, one line per equality:'
        'str, str / [ atom name 1 ], [ atom name 2] will have the same charge. '
        'Apart from that, all atoms of the same name (but possibly spread over '
        'different residues) will have the same charge enforced.')  
    
    parser.add_argument('outfile_csv', metavar='outfile.csv',
        help='Fitted charges will be written to a simple text file.')
    parser.add_argument('outfile_top', nargs='?', metavar='outfile.top', 
        default=None, help="GROMACS .top output file"
        "with updated charges according to given .hdf5")
    
    parser.add_argument('--qtot', '-q', default=0.0, type=float,
        help='The total charge of the system. [default=%(default)s]')
    parser.add_argument('--insertion-rules','-i',
        default="{'CD4':1,'CD3':1,'CA2':2,'CA3':2,'CB2':2,'CB3':2}",
        help="A string representation of a python dictionary, describing how "
        "many implicit hydrogens have been inserted at which atom."
        "Example and default: "
        "{'CD4':1,'CD3':1,'CA2':2,'CA3':2,'CB2':2,'CB3':2}")
    parser.add_argument('-v','--verbose', action='store_true',
        help="Prints a lot of information.")

    args = parser.parse_args()
    
    if args.verbose == True:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.WARNING

    logging.basicConfig(stream=sys.stdout, level=loglevel)  
    logging.info('Using replacement rules "{}"...'.format(args.insertion_rules))
    
    implicitHbondingPartners = ast.literal_eval(args.insertion_rules)
   
    #q, lagrange_multiplier, info_df, cg2ase, cg2cgtype, cg2q, sym2ase
    q, lagrange_multiplier, info_df, cg2ase, cg2cgtype, cg2q, sym2ase = \
    fitESPconstrained(infile_pdb = args.infile_pdb, 
              infile_top = args.infile_top, 
              infile_cost_h5 = args.infile_cost_h5, 
              infile_atoms_in_cg_csv = args.infile_atoms_in_cg_csv, 
              infile_cg_charges_csv = args.infile_cg_charges_csv, 
              infile_atoms_of_same_charge_csv = args.infile_atoms_of_same_charge_csv,
              qtot = args.qtot, strip_string=':SOL,CL', 
              implicitHbondingPartners = implicitHbondingPartners, 
              debug = args.verbose, outfile_top=args.outfile_top)
    
    np.savetxt(args.outfile_csv, q, delimiter=',')    
    
if __name__ == '__main__':
    main()
