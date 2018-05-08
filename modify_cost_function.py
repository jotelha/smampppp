#script which should read in hortons old cost function and add new constraints

# PLAN:
# Write functions which read input files for constraints and use them to
# construct the new cost function. After reading in and construction of the
# cost function, solve it and return the fitted charges.

# ENHANCEMENT:
# Check if the constraints are applied properly or if ther occured an error.
# Check if the constraints are meaningfull or if they exclude each other...

# TODO:
# maybe write a function which can visualize the symmetry constraints
# to proof the input. Highlight atoms which are constrained or sth.
# like this.


import numpy as np
import pandas as pd

### Function definitions ###

def constrainedMinimize(A_matrix, b_vector, C_scalar, D_matrix = None,
                        q_vector = np.array([0]), debug=False):
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

    N = b_vector.shape[0]
    M = q_vector.shape[0]

    npv = float(np.version.version[2:])

    if not isinstance(D_matrix,np.ndarray):
        D_matrix = np.atleast_2d( np.ones(N) )

    if debug:
        print("{:d} unknowns, {:d} equality constraints".format(N,M))
        print("A {}: \n {}".format(A_matrix.shape,A_matrix))
        print("B {}: \n {}".format(b_vector.shape,b_vector))
        print("C {}: \n {}".format(C_scalar.shape,C_scalar))
        print("D {}: \n {}".format(D_matrix.shape,D_matrix))
        print("q {}: \n {}".format(q_vector.shape,q_vector))
    if npv < 13:
        print '\nWARNING:\n Your numpy version {} is old,\n I am falling from '\
            'np.block() back to np.bmat()\n\n'.format(np.version.version)
        A = np.bmat([[ 2*np.atleast_2d(A_matrix), np.atleast_2d(D_matrix).T ],
                     [ np.atleast_2d(D_matrix), np.atleast_2d(np.zeros((M,M)))]])

    else:
        #maybe combine these steps to one np.block operation like for bmat
        A_upper = np.block(
            [ 2*np.atleast_2d(A_matrix), np.atleast_2d(D_matrix).T ])
        A_lower = np.block(
            [ np.atleast_2d(D_matrix), np.atleast_2d(np.zeros((M,M)))])
        A = np.block( [ [ A_upper ], [ A_lower ] ] )

    if debug:
        print("block A ({}): \n {}".format(A.shape,A))

    if npv < 13:
        B = np.bmat( [2*np.atleast_1d(b_vector), np.atleast_1d(q_vector)] ).T

    else:
        B = np.block( [2*np.atleast_1d(b_vector), np.atleast_1d(q_vector)] )

    if debug:
        print("block B ({}): \n {}".format(B.shape,B))

    C = C_scalar

    x = np.linalg.solve(A, B)

    return x, A, B


def constructPairwiseSymmetryConstraints(charges, N, symmetry=1.0, debug=False):
    """
    Function to construct D_matrix and q_vector for a pairwise symmetry

    Parameters:
    -----------
    charges: 1D or 2D numpy.ndarray int
        List of charge indices which should be equal. For more than two charges
        there will be (N-1) constraints for pairwise equal charges. if you give
        a list of lists with equal charges each sublist will be forced to have
        equal charges.
    N: int
        Total number of atoms in your system. Needed to fix the size of the
        D_matrix and q_vector to the system size (number of atoms)
    symmetry: int, -1 or 1
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
        symmetries with an symmetry array.
    """

    charges = np.atleast_2d(charges)
    D = np.ones((1, N))
    q = np.ones((1))
    for charge_list in charges:
        M = len(charge_list)-1

        symmetry = symmetry*np.ones(M)

        if debug:
            print("{:d} unknowns, {:d} pairwise equality constraints".format(N,M))
            print("symmetry list ({}):\n{}".format(symmetry.shape,symmetry))

        D_single = np.atleast_2d(np.zeros((M,N)))
        q_single = np.atleast_1d(np.zeros(M))
        D_single[:,charge_list[0]] = 1

        for j in range(M):
            D_single[j,charges[j+1]] = -1.0*symmetry[j]

        if debug:
            print("D_single ({}):\n{}".format(D_single.shape,D_single))
            print("q_single ({}):\n{}".format(q_single.shape,q_single))

        #add up D_single and q_single to D and q
        D = np.append(D, D_single, axis=0)
        q = np.append(q, q_single, axis=0)

    return D[1:], q[1:]

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
#        print("{:d} unknowns, {:d} pairwise equality constraints".format(N,M))
#        print("symmetry list ({}):\n{}".format(symmetry.shape,symmetry))
#
#    D = np.atleast_2d(np.zeros((M,N)))
#    q = np.atleast_1d(np.zeros(M))
#    D[:,charges[0]] = 1
#
#    for j in range(M):
#        D[j,charges[j+1]] = -1.0*symmetry[j]
#
#    if debug:
#        print("D ({}):\n{}".format(D.shape,D))
#        print("q ({}):\n{}".format(q.shape,q))
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

    M = len(chargeGroups)

    q_vector = np.atleast_1d(q*np.ones(M))

    if debug:
        print("{:d} unknowns, {:d} pairwise equality constraints".format(N,M))

    D_matrix = np.atleast_2d(np.zeros((M,N)))
    #q = np.atleast_2d(np.zeros(M))

    for j in range(M):
        D_matrix[j,chargeGroups[j]] = 1.0

    if debug:
        print("D_matrix ({}):\n{}".format(D_matrix.shape, D_matrix))
        print("q_vector ({}):\n{}".format(q_vector.shape, q_vector))

    return D_matrix, q_vector


def read_Name_Index_ChargeGroup(file_name):
    """
    Function to read in csv file of names and charge groups with the help of
    pandas. It translates the 'pdb' names into ase structure atom indices and
    constructs three vectors of atom names, indices and associated charge
    group. The vectors are all ordered in the same manner. It also returns
    cg_max, to reconstruct the charge group if atom names occure more than ones.

    Parameters
    ----------
    file_name: str
        name of the file which is read

    Return
    ------
    name: numpy.array (str)
        array with atom names
    index: numpy array (int)
        array with atom indices corresponding to the atom names (same order)
    cg: numpy.array (int)
        array of charge groups coressponding to the atom names and indices
        (same order). If atoms with the same name occure more than ones in the
        atom names they will be given different charge groups, by:
        cg = charge group + occurence * (cg_max+1)
        Obviously the associated charge group can be recived by: cg % (cg_max+1)
    cg_max: int
        the maximum number of the old charge groups which are used for more than
        one atom index and thus carry an ambiguous meaning.

    TODO
    ----
        --> read 'index_name' from the structure files (.traj, .pdb, ...) and
            not from a '.txt' file! Use some skript from Johannes/ask him...
    """

    name_cg = pd.read_csv(file_name, sep=',', header=None).values[2:]
    # read names only up to 103 the cutted names are the added H2Os.
    name = np.loadtxt('test_files/name_index.txt', dtype='str')[:103]

    index = np.arange(len(name))
    cg_max = max(np.array(name_cg[:,1],dtype=int)) #maximum charge group number
    cg = []
    for i, n in enumerate(name):
        cg_value = np.where(name_cg[:,0] == n)
        # introduce new cg group for atoms which occure more than ones in the
        # molecule.
        occurence_of_name = len(np.where(name[:i].T == n)[0])
        # cg_group + n*(cg_max*1) to avoid doubly named groups
        cg.append(int(name_cg[cg_value,1][0][0]) + (cg_max+1)*occurence_of_name)

    return name, index, cg, cg_max


def read_ChargeGroup_TotalCharge(file_name):
    """
    Function to read in csv file of charge groups and charges with the help of
    pandas.

    Parameters
    ----------
    file_name: str
        name of the file which is read

    Return
    ------
    cg_q: numpy.ndarray (int)
        array of charge group number and corresponding total charge.
        the charge group is cg_q[0] and the total charge cg_q[1].
    """

    cg_q = pd.read_csv(file_name, sep=',', header=None).values[1:].T
    cg_q = cg_q.astype(dtype=int)

    return cg_q


def read_SameChargedAtoms(file_name):
    """
    Function to read in csv file of atoms which should have the same charge.

    Parameters
    ----------
    file_name: str
        name of the file which is read

    Return
    ------
    sca_by_index: 2D numpy.ndarray (int)
        array of atoms with same charges for different symmetries(+1/-1). The
        same charged atoms (sca) are identified by their indices.
    """

    sca = pd.read_csv(file_name, sep=',', header=[2]).values #same charged atoms

    name = np.loadtxt('test_files/name_index.txt', dtype='str')[:103]

    sca_by_index = [] #same charged atoms ordered by indices of atoms
    for group in sca:
        indices = []
        for atom_name in group:
            indices.append(np.where(name == atom_name)[0])
        sca_by_index.append(np.array(indices).flatten())

    sca_by_index = np.asarray(sca_by_index)

    return sca_by_index


def concatenated_constraints(D_matrices, q_vectors):
    """
    Function to concatenate D_matrices and q_vectors to one D_matrix and one
    q_vector by using numpy.hstack and numpy.vstack. The order of D_matrix in
    D_matrices and q_vector in q_vectors should be the same, otherwise the
    constraints are connected wrong.

    Parameters
    ----------
    D_matrices: numpy.ndarray (reals)
        array of all D_matrices which should be concateneted
    q_vectors: numpy.ndarray (reals)
        array of all q_vectors which should be concateneted

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
        print("A: {}".format(A_horton))
        print("B: {}".format(B_horton))
        print("C: {}".format(C_horton))
        print("N: {}".format(N_horton))

    return A_horton, B_horton, C_horton, N_horton



### ACTUAL PROGRAM ###
#--------------------#

### Read Hortons Cost Function
A_horton, B_horton, C_horton, N_horton = read_horton_cost_function(
    file_name = 'test_files/smamp_esp_cost_ua.h5')

### Charge Groups:
# read in all charge groups and construct the corresponding constraints
name, index, cg, cg_max = read_Name_Index_ChargeGroup(
    file_name = 'test_files/Atom_name_charge_group.csv')
cg_q = read_ChargeGroup_TotalCharge(
    file_name = 'test_files/charge_group_total_charge.csv')

#loop over set of charge groups (each charge group occures only ones)
c_groups = []
charges  = []
for cg_value in set(cg):
    #step with index only necessary if index is not ordered
    indices = index[np.where(cg == cg_value)[0]]
    c_groups.append(indices)
    c = cg_q[1, np.where(cg_q[0] == cg_value%(cg_max+1))][0][0]
    charges.append(c)

#print charges, c_groups

D_matrix_cg, q_vector_cg = constructChargegroupConstraints(
    chargeGroups = c_groups, N = N_horton, q = charges, debug=False)

#print D_matrix, D_matrix.shape
#print q_vector, q_vector.shape


### Same Charged Atoms
sca_by_index = read_SameChargedAtoms(
    file_name='test_files/Atoms_have_the_same_charges.csv')

print sca_by_index

D_matrix_sym, q_vector_sym = constructPairwiseSymmetryConstraints(
    charges = sca_by_index, N = N_horton, symmetry = 1.0, debug = False)

#print D_matrix_sym, D_matrix_sym.shape
#print q_vector_sym, q_vector_sym.shape


### Concatenate The Constraints
D_ms = np.append(D_matrix_cg, D_matrix_sym, axis=0) #D_matrices
q_vs = np.append(q_vector_cg, q_vector_sym)         #q_vectors
D_matrix_all, q_vector_all = concatenated_constraints(D_matrices = D_ms,
                                                      q_vectors = q_vs)

print D_matrix_all, D_matrix_all.shape
print q_vector_all, q_vector_all.shape


### Constrained Minimization
X, A, B = constrainedMinimize(A_matrix = A_horton,
                        b_vector = B_horton,
                        C_scalar = C_horton,
                        D_matrix = D_matrix_all,
                        q_vector = q_vector_all,
                        debug    = False)

print 'Results:'
#prevent scientific notation and make the prints mor readable
np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)
print 'charges {}:\n {}\ncharge sum = {}\n'.format( X[:N_horton].T.shape,
                                                    X[:N_horton].T,
                                                    X[:N_horton].T.sum() )
print 'Lagrange multipliers {}:\n {}'.format( X[N_horton:].T.shape,
                                              X[N_horton:].T )


### test the results
print 'optimized result: ',\
    (np.dot(X.T, np.dot(A, X)) - 2*np.dot(B.T, X) - C_horton)[0,0]
