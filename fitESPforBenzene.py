#!/usr/bin/env python
""" Estimate point charges for Benzene without constraints """

from fitESPconstrained import *
from ase.build import molecule

def main():
    import argparse    

    parser = argparse.ArgumentParser(prog='fitESPforBenzene.py',
        description='Estimate point charges for Benzene without constraints')
    parser.add_argument('infile_cost_h5',metavar='cost.h5',
        help='The location of the HORTON cost function in the form '
             '"file.h5:group/cost". This argument must be the same as the '
             'output argument of the script horton-esp-cost.py.')
    args = parser.parse_args()

    A_horton, B_horton, C_horton, N_horton = \
        read_horton_cost_function(file_name = args.infile_cost_h5)    

    X_unconstrained, A_unconstrained, B_unconstrained = \
            unconstrainedMinimize(A_matrix = A_horton,
                        b_vector = B_horton,
                        C_scalar = C_horton)

    benzene = molecule('C6H6')
    for i, a in enumerate(benzene):
        a.charge = X_unconstrained[i]
    for a in benzene:
        print("{}, {: .3f}".format(a.symbol,a. charge))

if __name__ == '__main__':
    main()


