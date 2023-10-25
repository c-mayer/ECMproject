#!/usr/bin/env python

"""Python script to validate results of the developed script with the results of the predecessor ecmtool.

@author: Christian Mayer
"""

import numpy as np
import argparse


def read_in_matrix(file, fileformat, objective=True):
    """Reads in matrix of the given ecmtool result file."""
    if fileformat == 'csv':
        delimiter = ','
    if fileformat == 'txt':
        delimiter = '\t'
    matrix = np.genfromtxt(file, delimiter=delimiter)
    if np.isnan(matrix[0]).all(): # proofs if all elements of the first row are nan --> if header in matrix then True
        matrix = matrix[1:] # take out header
    if objective and "ecmtool" in file: # proofs if objective function should be taken out and if file is an ecmtool result file
        matrix = matrix[:,:-1]
    return matrix


def normalize_array(array, column):
    """Takes an array and normalizes each row with the corresponding absolute value of the given column (column as integer)."""
    return np.array([row / abs(row[column]) for row in array])


def proof_for_uniqueness(matrix, tolerance=2e-5, matrix_non_unique=False):
    """Proofs if all ECMs of matrix are only existing once. If matrix_non_unique=True, returns all non unique ECMs."""
    ub = matrix + tolerance # upper border
    lb = matrix - tolerance # lower border
    res = np.logical_and((ub[:, None] >= matrix), (lb[:, None] <= matrix)).all(-1) # kind of identity matrix if row just matches with itself
    if np.array_equal(np.array(res, dtype=int), np.identity(np.shape(matrix)[0], dtype=int)): # proofs if res is an identity matrix
        print("All ECMs are unique.")
    # creates list with all indices of doulbe entries
    # counts all True and False of ECM --> np.unique(ECM, return_counts=True)
    # proofs if amount of True is greater than one (shouldnt be if every entry just matches with itself)
    index_list = [i for i, ECM in enumerate(res) if np.unique(ECM, return_counts=True)[1][1] > 1]
    if len(index_list) > 0:
        print(f"There are {len(index_list)} ECMs not unique.")
        if matrix_non_unique:
            non_unique_matrix = matrix[index_list]
            print(non_unique_matrix)
            return non_unique_matrix
    

def validate_fast(matrix_1, matrix_2, tolerance=1e-5):
    """Compares rows of matrix_1 with rows of matrix_2. Returns numpy array with all non_matching ECMs ."""
    ub = matrix_1 + tolerance # upper border
    lb = matrix_1 - tolerance # lower border
    res = np.logical_and((ub[:, None] >= matrix_2), (lb[:, None] <= matrix_2)).all(-1).any(-1)
    # res is vector with bool value for each row of matrix_1 --> True if row also in matrix_2
    non_matching = matrix_1[np.logical_not(res)]
    if res.all(): # all values of vector are true
        print(f"All ECMs are the same.")
    else:
        print(f"{len(non_matching)} ECMs are not matching.")
    return non_matching


def validate_verbose(matrix_1, matrix_2, tolerance=1e-5):
    """Proofs ECMs on appearance in both given matrices. All ECMs of both matrices which dont match are returned as tuple (non_1, non_2)."""
    non_ECM_1 = []
    indices = []
    # comparison of the rows of matrices
    for i_1, row in enumerate(matrix_1):
        success = 0 # counts how often ECM in matrix_1 matches with ECM in matrix_2
        ob = row + tolerance
        lb = row - tolerance
        for i_2, comparison in enumerate(matrix_2):
            if (ob > comparison).all() and (lb < comparison).all():
                success += 1
                indices.append(i_2) # remember which rows of matrix_2 matched with rows of matrix_1     
        # output matrix_1
        if success == 0:
            non_ECM_1.append(row)
        elif success > 1:
            print(f"We have a one to many problem!")
            print(f"ECM with index {i_1} of matrix_1 was matched by {success} ECMs of matrix_2.")

    if len(non_ECM_1) > 0:
        print(f"{len(non_ECM_1)} ECMs of matrix_1 are not the same.")
    else:
        print(f"All ECMs of matrix_1 are the same.")
    non_matching_1 = np.array(non_ECM_1)

    # output matrix_2
    # create dictionary with indices of matrix_2
    ind_dict = {}
    for key in range(np.shape(matrix_2)[0]):
        ind_dict[key] = 0
    # count how often each index in indices list
    for i in indices:
        ind_dict[i] += 1
    # loop through dictionary
    n = 0
    fail = 0
    fail_indices = [] # list wit all indices of matrix_2, not matched matrix_1
    for key, value in ind_dict.items():
        if value > 1:
            print(f"We have a many to one problem!")
            print(f"ECM with index {key} of matrix_2 was matched by {value} ECMs of matrix_1.")
        elif value == 0:
            fail += 1
            fail_indices.append(key)
    # output for matrix_2
    if fail > 0:
        print(f"{fail} ECMs of matrix_2 are not the same.")
    else:
        print(f"All ECMs of matrix_2 are the same.")
    non_matching_2 = matrix_2[fail_indices]
    
    # give matrices as output
    return (non_matching_1, non_matching_2)


def calculate_percentages(matrix, non_matching_matrix):
    """Takes original matrix and matrix of non-matches from validate function and calculates percentages of outcomes."""
    n_total = len(matrix)
    n_non = len(non_matching_matrix)
    n = n_total - n_non
    print(f"From {n_total} ECMs:")
    print(f"{n} ECMs are matching.")
    print(f"{(n/n_total)*100} % of ECMs are matching.")
    print(f"{n_non} ECMs are not matching.")
    print(f"{(n_non/n_total)*100} % of ECMs are not matching.")


if __name__ == '__main__':
    
    ### argparse ###
    parser = argparse.ArgumentParser(description="Python script to validate results of the developed script with the results of the predecessor ecmtool. @author: Christian Mayer", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # create group for required arguments
    parser_req = parser.add_argument_group('required arguments')

    parser_req.add_argument('-f1', '--file1',
                            help='Enter input result file.',
                            type=str,
                            metavar='PATH_TO_FILE',
                            action='store', 
                            required=True)
    parser_req.add_argument('-f2', '--file2',
                            help='Enter input result file.',
                            type=str,
                            metavar='PATH_TO_FILE',
                            action='store', 
                            required=True)

    # optional arguments
    parser.add_argument('-ff1', '--fileformat1',
                        help='Give fileformat for file 1. (choices: csv or txt)',
                        type=str,
                        metavar='STR',
                        choices=['csv', 'txt'],
                        action='store',
                        default='csv')
    parser.add_argument('-ff2', '--fileformat2',
                        help='Give fileformat for file 2. (choices: csv or txt)',
                        type=str,
                        metavar='STR',
                        choices=['csv', 'txt'],
                        action='store',
                        default='csv')
    parser.add_argument('-n', '--norm_column',
                        help='Give column on which both input files get normalized. Column should not contain any zeros.',
                        type=int,
                        metavar='INT',
                        action='store',
                        default=0)
    parser.add_argument('-v', '--verbose',
                        help='Results will be calculated with slower but more detailed function. Output will give more information.',
                        action='store_true')
    parser.add_argument('-ob', '--objective',
                        help='Set "True", if objective column exists in ecmtool results.',
                        type=bool,
                        metavar='BOOL',
                        action='store',
                        default=True)
    parser.add_argument('-tol_val', '--tolerance_validation',
                        help='Give tolerance, when two numbers are seen as the same for validation.',
                        type=float,
                        metavar='FLOAT',
                        action='store',
                        default=1e-5)
    parser.add_argument('-tol_un', '--tolerance_unique',
                        help='Give tolerance, when two numbers are seen as the same for proof of uniqueness. Useful value: double of the validation tolerance',
                        type=float,
                        metavar='FLOAT',
                        action='store',
                        default=2e-5)
    parser.add_argument('-mnu', '--matrix_non_unique',
                        help='Set "True", if you want a matrix with all non unique entries as additional output.',
                        type=bool,
                        metavar='BOOL',
                        action='store',
                        default=False)
    
    args = parser.parse_args()

    ### main ###
    # load in matrices
    matrix_1 = read_in_matrix(args.file1, args.fileformat1, objective=args.objective)
    matrix_2 = read_in_matrix(args.file2, args.fileformat2, objective=args.objective)

    # normalize column values rowwise on column <norm_column> of both matrices
    matrix_1 = normalize_array(matrix_1, args.norm_column)
    matrix_2 = normalize_array(matrix_2, args.norm_column)

    # proof each matrix on unique entries
    print(f"{args.file1}:")
    print(f"shape of matrix: {np.shape(matrix_1)}")
    proof_for_uniqueness(matrix_1, tolerance=args.tolerance_unique, matrix_non_unique=args.matrix_non_unique)
    print(f"{args.file2}:")
    print(f"shape of matrix: {np.shape(matrix_2)}")
    proof_for_uniqueness(matrix_2, tolerance=args.tolerance_unique, matrix_non_unique=args.matrix_non_unique)

    # compare ECMs of both files
    if args.verbose:
        non_matching_1, non_matching_2 = validate_verbose(matrix_1, matrix_2, tolerance=args.tolerance_validation)
    else:
        print(f"{args.file1}:")
        non_matching_1 = validate_fast(matrix_1, matrix_2, tolerance=args.tolerance_validation)
        print(f"{args.file2}:")
        non_matching_2 = validate_fast(matrix_2, matrix_1, tolerance=args.tolerance_validation)

    # print percentages of file_1 compared to file_2
    print(f"{args.file1}:")
    calculate_percentages(matrix_1, non_matching_1)

    # print percentages of file_2 compared to file_1
    print(f"{args.file2}:")
    calculate_percentages(matrix_2, non_matching_2)
    