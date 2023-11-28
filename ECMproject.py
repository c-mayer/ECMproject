#!/usr/bin/env python

"""
Python script to calculate ECM's from given metabolic model. Reads in metabolic model in xml format and creates list of all ECMs as csv file in default settings. If you also want intermediate result for development, set developer flag. See help page (-h) for further information.

@author: Christian Mayer, Marcus Holzer
"""

### import statements ###    
import argparse
import cobra
import numpy as np
import pandas as pd
from cobra.io import read_sbml_model
import io
from contextlib import redirect_stderr
import os
import subprocess
import time as ti
from fractions import Fraction
import json
import sys
import multiprocessing as mp
from itertools import islice
import gzip
from functools import partial


### functions ###

def read_model(input_filename):
    """
    Reads metabolic model from sbml file using the 'cobra.io.read_sbml_model' functions. Reads io string during reading
    model and catches exchange reactions added by cobrapy. These reactions are going to be removed from the model again.
    :param input_filename: sbml file with cobrapy compatible metabolic model
    :return:
        - model - cobrapy model
        - reas_added - list of reactions added by cobrapy
    """
    
    stream = io.StringIO()
    with redirect_stderr(stream):
        model = cobra.io.read_sbml_model(input_filename)
    console_output = stream.getvalue()

    reas_added = []
    for text in console_output.split('\n'):
        if text.startswith("Adding exchange reaction"):
            tmp = text.split()
            reas_added.append(tmp[3])

    return model, reas_added


def rm_reactions(model, rea_list):
    """
    Removes exchange reaction that were added by cobrapy.
    :param model: cobrapy model
    :param list rea_list: list of reaction names that will be removed
    :return:
        - model - altered cobrapy model
    """
    model.remove_reactions(rea_list)
    return model


def indicate_exchange(model, extern_compartments):
    """Iterates over each reaction in the cobra model and adds an reaction.exchange attribute (value: True or False) to each reaction.
    This attribute is used later in the code to identify exchange reactions.
    :param model: cobrapy model
    :param list extern_compartments: list of extern compartments created by function get_extern_compartments
    :return: None
    """
    # if no SBO term and products and reactands on reaction --> cobrapy cant detect it as extern
    # we take each reaction as external if:
    # SBO:0000627
    # OR
    # reaction id starts with 'R_EX_' or 'EX_'
    # OR
    # reaction is one sided (reaction.boundary = True) and no demand or sink (we proof for demand and sink before we proof for boundary)
    # if exchange reaction, than reaction.exchange = True
    for reaction in model.reactions:
        # reaction id starts with 'R_EX_' or 'EX_'
        if 'EX_' in reaction.id[:5]:
            print(f'exchange: {reaction.id}')
            reaction.exchange = True
            continue

        # unfortunately, demand reactions do have the same SBO Term as exchange reactions
        # demand reactions differ from exchange reactions in case of exchange reactions work with metabolites of external compartment
        # we could use that (if exchange sbo term --> proof if one metabolite is in external compartment):
        if 'sbo' in reaction.annotation:
            if reaction.annotation['sbo'] == 'SBO:0000627':
                extern = False
                for metabolite in reaction.metabolites:
                    # compartment name needs to be given and contain 'ex' or compartment id needs to contain 'e'
                    if (metabolite.compartment in extern_compartments and model.compartments[metabolite.compartment]) or 'e' in metabolite.compartment:
                        extern = True
                if extern:
                    # add reaction to exchange reactions
                    print(f'exchange: {reaction.id}')
                    reaction.exchange = True
                    continue


        # kill all demand reactions by id prefix 'R_DM_'
        if 'DM_' in reaction.id[:5]:
            print(f'demand: {reaction.id}')
            reaction.exchange = False
            continue
        # kill all sink reactions by SBO Term
        if 'sbo' in reaction.annotation:
            if reaction.annotation['sbo'] == 'SBO:0000632':
                print(f'sink: {reaction.id}')
                reaction.exchange = False
                continue
        # kill all sink reactions by id prefix 'R_SINK_'
        if 'SINK_' in reaction.id[:6]:
            print(f'sink: {reaction.id}')
            reaction.exchange = False
            continue

        # proof if one sided
        if reaction.boundary:
            # since we already ended iteration if criterium was for sink or demand, all remaining reaction.boundary = True reactions hopefully are exchange reactions
            print(f'exchange: {reaction.id}')
            reaction.exchange = True
            continue
        # if no exchange reaction (iteration reaches this point)
        reaction.exchange = False


def mark_exchange(model, extern):
    """
    Takes all reaction ids of the model and compares them with the given list of extern reactions. If matched, exchange attribute of reaction is set to True.
    :param model: cobrapy model
    :param list extern: list of reaction ids which will be compared
    :return: None
    """
    for reaction in model.reactions:
        if reaction.id in extern:
            reaction.exchange = True
        else:
            reaction.exchange = False


def write_reaction_direction(model):
    """Loops through all reactions of a given cobrapy model and writes an attribute "direction" depending on the reaction borders."""
    for reaction in model.reactions:
        if reaction.lower_bound < 0 and reaction.upper_bound > 0:
            reaction.direction = "both"
        elif reaction.upper_bound > 0 and reaction.lower_bound >= 0:
            reaction.direction = "forward"
        elif reaction.lower_bound < 0 and reaction.upper_bound <= 0:
            reaction.direction = "backward"
        else:
            reaction.direction = False
            print(f'Reaction {reaction.id} has no feasible boundaries!')


def correct_stochiometric_matrix(model, smatrix):
    """Changes algebraic sign of each column in stochiometric matrix with backward reaction."""
    for i, reaction in enumerate(model.reactions):
        if reaction.direction == "backward":
            smatrix[:,i] = smatrix[:,i] * -1
    return smatrix


def ex_reaction_id_and_index(model):
    """
    Creates a list of tuples from a given cobrapy model. Each tuple consists of exchange reaction id and index in order of all reactions. Needed for mplrs project.
    """
    return [(reaction.id, i) for i, reaction in enumerate(model.reactions) if reaction.exchange is True]


def get_extern_compartments(model):
    """Writes all compartments in list which have "ex" in their name or no name to get all extern compartments of a given cobrapy model. Takes also all compartments with no name, because metabolites get filtered after that for id in oter function."""
    # we list all compartment ids of the model which have "ex" in their name or no name
    return [compartment for compartment in model.compartments if "ex" in model.compartments[compartment].lower() or not model.compartments[compartment]]


def ex_metabolites_id_and_index(model, extern_compartments):
    """Creates for each extern metabolite which is attched to an external reaction a tuple. Each tuple consists of metabolite id and index in order of the model. 
    Needed for header of final results and for slicing stochiometric matrix."""
    # want to get all external metabolites
    # mplrs kicked out all internal ones already but our header does not know which
    # so we need to kick again
    new_ex_metabolites = []
    for reaction in model.reactions:
        if reaction.exchange:
            if reaction.reactants and reaction.products:
                # if exchange reactions are both sided in the model, we search only for metabolites with extern compartment
                # id of extern compartment needs to contain an "e" and compartment name needs "ex" or nothing, otherwise it will not work
                new_ex_metabolites += [metabolite.id for metabolite in reaction.metabolites.keys() if "e" in metabolite.compartment and metabolite.compartment in extern_compartments]
            else:
                new_ex_metabolites += [metabolite.id for metabolite in reaction.metabolites.keys()]

    # in case some external reactions share the same external metabolite and bring it in right order
    # right order is not necessary for ex_metabolites
    new_ex_metabolites = set(new_ex_metabolites)
    return [(metabolite.id, index) for index, metabolite in enumerate(model.metabolites) if metabolite.id in new_ex_metabolites]


def split_extern_reversible_reactions(model):
    """Splitting all reversible extern reactions of a cobrapy model into two cobrapy reaction objects."""
    for reaction in model.reactions:
        if reaction.reversibility and reaction.exchange: # since splitted reactions are not reversible, they dont get splitted again
            # create backward irreversible reaction from reversible reaction
            backward_reaction = cobra.Reaction(reaction.id + "_b")
            backward_reaction.name = reaction.name # reaction name is the same by purpose (remerging if name is the same)
            backward_reaction.subsystem = reaction.subsystem
            backward_reaction.lower_bound = 0.  # make it irreversible
            backward_reaction.exchange = True
            if 'sbo' in reaction.annotation:
                backward_reaction.annotation['sbo'] = reaction.annotation['sbo']
            # add reaction to model
            model.add_reactions([backward_reaction])
            # add metabolites to reaction
            metabolite_dict = reaction.metabolites
            for object in metabolite_dict:
                backward_reaction.add_metabolites({object: (metabolite_dict[object] * -1)})

            # alter forward reaction to split
            reaction.id = reaction.id + "_f"
            #reaction.name = reaction.name
            reaction.lower_bound = 0


def as_fraction(number, approximation=None):
    """
    Takes integers or floats and converts them to rational strings.
    
    Parameters
    ----------
    number : int, float
        Give number for conversion to a rational.
    approximation : int default None
        Sets the border for number a denominator can have.
        If border of 1e12 is set, you get a rational unequal to zero for numbers down to 1e-12 and zero for numbers lower than that.
        Approximation default is 1e6. Keep attention. Less strict approximation borders can lead to huge numbers which slow down the calculation intensely.
    
    Returns
    -------
    str
        A rational or integer number as a string.
    
    Examples
    --------
    >>> into_fractions(5)
    '5'
    >>> into_fractions(0)
    '0'
    >>> into_fractions(0.5)
    '1/2'
    >>> into_fractions(1/6)
    '1/6'
    >>> into_fractions(1/6, approximation=1e20)
    '6004799503160661/36028797018963968'
    >>> into_fractions(1/1e6)
    '1/1000000'
    >>> into_fractions(1/1e7)
    '0'
    """
    if approximation:    
        return str(Fraction(number).limit_denominator(int(float(approximation))))
    else:
        return str(Fraction(number).limit_denominator())


def write_h_representation(smatrix, model, tmp_dir, core_name, approximation=None):
    """
    Takes stochiometric matrix of cobra model and converts it into an .ine file with H-representation.
    
    :param smatrix: stochiometric matrix created by cobrapy
    :param model: takes cobrapy model of the loaded network
    :param core_name: filename; created file will be named after core_name + .ine
    :param approximation: number; if value is given, floats will get "rounded" to fractions with max this number as denominator (e.g. 1e6)
    :return: None
    """

    # dimension variables
    stoich = len(smatrix) # n_rows of stoichiometric matrix
    n_reactions = len(model.reactions) # n_columns of matrix

    # calculate number of reversible and irreversible reactions
    n_rev_reactions = 0
    for reaction in model.reactions:
        if reaction.reversibility:
            n_rev_reactions += 1

    n_irrev_reactions = n_reactions - n_rev_reactions


    ### create .ine file ###
    file = open(tmp_dir + core_name + ".ine", "w")

    # write header ine file
    file.write("* " + core_name + "\n")
    file.write("H-representation" + "\n")

    # linearity
    file.write("linearity " + str(stoich))
    for i in range(1, stoich + 1):
        file.write(" " + str(i))
    file.write("\n")

    # begin of matrix
    file.write("begin" + "\n")
    # write shape of matrix and rational
    file.write(str(stoich + n_irrev_reactions) + " " + str(len(smatrix[0]) + 1) + " rational \n")

    # write matrix into file
    # writes stochiometric matrix into file
    for line in smatrix:
        file.write(str(0) + " ") # b column
        for val in line:
            file.write(as_fraction(val, approximation=approximation) + " ")
        file.write("\n")

    # create Irreversibility constrain matrix
    for index, reaction in enumerate(model.reactions):
        if not reaction.reversibility:
            file.write(str(0) + " ") # b column
            for index2 in range(0, n_reactions):
                if index == index2:
                    file.write(str(1) + " ")
                else:
                    file.write(str(0) + " ")
            file.write("\n")

    file.write("end" + "\n")

    file.close()


def redund(n_processes, path_mplrs, tmp_dir, core_name, verbose=True):
    """Performs mplrs redund on given .ine file."""
    #original_cmd = "mpirun -np 3 /opt/lrslib/v072/mplrs -redund ./h_representation.ine > input_postredund.ine"
    cmd = ["mpirun", "-np", str(n_processes), path_mplrs, tmp_dir + core_name + ".ine", tmp_dir + core_name + "_postredund.ine", "-redund"]
    if verbose:
        subprocess.run(cmd)
    else:
        subprocess.run(cmd, stdout=subprocess.DEVNULL)


def rewrite_postredund_file(tmp_dir, core_name):
    """Rewrite postredund.ine file."""
    # delete last seven rows
    file = open(tmp_dir + core_name + "_postredund.ine", "r+")
    lines = file.readlines() # save all lines in object
    file.seek(0) # go to start of file with pointer
    file.truncate() # delete all lines from pointer downwards (delete whole file)
    file.writelines(lines[:-7])
    file.close()


def write_project_line(ex_reactions, tmp_dir, core_name):
    """Write mplrs project line to postredund.ine file."""
    file = open(tmp_dir + core_name + "_postredund.ine", "a")
    # write project line into file
    file.write("project " + str(len(ex_reactions)))
    for reaction, i in ex_reactions:
        file.write(" " + str(i + 1))
    file.write("\n")
    file.close()


def mplrs_project(smatrix, ex_reactions, n_processes, path_mplrs, tmp_dir, core_name, rows=60, lastp=10, lastrows=10, verbose=True):
    """Performs mplrs project on postredund.ine file."""
    # calculate n times
    # if we want the last H-representation
    n = len(smatrix[0]) - len(ex_reactions)
    
    print(f'Start projection of postredund file.')
    # mplrs project
    cmd = ["mpirun", "-np", str(n_processes), path_mplrs, tmp_dir + core_name + "_postredund.ine", tmp_dir + core_name + "_h.projected", "-rows", str(rows),  "-lastp", str(lastp), "-lastrows", str(lastrows)]
    for i in range(0, n):
        print(f'Run {i + 1}/{n}')
        if verbose:
            subprocess.run(cmd)
        else:
            subprocess.run(cmd, stdout=subprocess.DEVNULL)
        os.replace(tmp_dir + core_name + "_h.projected", tmp_dir + core_name + "_postredund.ine")
    os.replace(tmp_dir + core_name + "_postredund.ine", tmp_dir + core_name + "_h.projected")


def mfel_project(mfel_file, n_processes, tmp_dir, core_name, rows=20, lastp=20, lastrows=5, verbose=True):
    """Performs mplrs project on postredund.ine file with mfel.tcsh script. Its a modified fel script of the creators of mplrs."""
    print(f'Start projection of postredund file.')
    # give execute file permission
    cmd = ["chmod", "+x", mfel_file]
    subprocess.run(cmd)
    
    # run mfel file
    # we need to run this with all first three positions, otherwise it wont work properly ($1, $2, $3 --> $# >= 3)
    cmd = [mfel_file, tmp_dir + core_name + "_postredund.ine", tmp_dir + core_name + '_h.projected', str(n_processes), str(rows), str(lastp), str(lastrows)] # default cores are all
    if verbose:
        subprocess.run(cmd)
    else:
        subprocess.run(cmd, stdout=subprocess.DEVNULL)


def mplrs_conversion(n_processes, path_mplrs, tmp_dir, core_name, verbose=True):
    """Convert H-representation to V-representation via mplrs."""
    cmd = ["mpirun", "-np", str(n_processes), path_mplrs, tmp_dir + core_name + "_h.projected", tmp_dir + core_name + ".projected"]
    print(f'Run conversion')
    if verbose:
        subprocess.run(cmd)
    else:
        subprocess.run(cmd, stdout=subprocess.DEVNULL)


def slice_stochio_matrix(smatrix, ex_reactions, ex_metabolites):
    """
    Slice stochiometric matrix to a form where it can be calculated with our merged V-representation matrix.
    Takes list with tuples of extern reaction (reaction.id, index) and list with tuples of metabolites attached to extern reactions (metabolite.id, index) and slices stochiometric matrix with these indices.
    Returns sliced stochiometric matrix.
    """
    # get indices for slicing stochiometric matrix
    reaction_indices = [i for reaction, i in ex_reactions]
    metabolite_indices = [i for reaction, i in ex_metabolites]
    # new stochiometric matrix with m * r_ex (m = metabolite; r_ex = exchange reaction)
    # first: take out all non-external reactions
    smatrix = smatrix[:,reaction_indices]
    # second: take out all metabolites not attached to any extern reaction
    smatrix = smatrix[metabolite_indices]
    # alternative: smatrix = smatrix[~np.all(smatrix == 0, axis=1)]
    return smatrix

    
def merge_model(model):
    """Merges together splitted reactions from function split_extern_reversible_reactions in given cobrapy model and returns indices of reaction pairs to merge on pre-ECMs."""
    # create a "list" of the exchange reactions with new indices for after the projection
    ex_reactions_projected = []
    count = 0
    for reaction in model.reactions:
        if reaction.exchange:
            ex_reactions_projected.append((count, reaction))
            count += 1
    # use the new list of exchange reactions and find the reactions that were split in the preprocessing step
    # create a list of reaction pair indices for the reactions that will be merged
    # use the merge_reactions function to merge the plit reactions
    reaction_pair_index = []
    for index_1, reaction_1 in ex_reactions_projected:
        if reaction_1.id[-2:] == '_f':
            for index_2, reaction_2 in ex_reactions_projected:
                if reaction_1.name == reaction_2.name and reaction_1.id != reaction_2.id:
                    reaction_pair_index.append((index_1, index_2))
                    merge_reactions(reaction_1, reaction_2)
    return reaction_pair_index


def merge_reactions(reaction_f, reaction_b):
    """Merges two reaction objects together. Rewrites the forward reaction to original state and deletes backward reaction."""
    # rewrite forward reaction
    reaction_f.id = reaction_f.id[:-2]
    reaction_f.lower_bound = -1000
    # delete backward reaction
    model.remove_reactions([reaction_b])


def write_output_header(ex_metabolites, separator=','):
    """Creates header for the output file."""
    header = ''
    count = 0
    stop = len(ex_metabolites)
    for metabolite, index in ex_metabolites:
        header += metabolite
        count += 1
        if count == stop:
            break
        header += separator
    return header


def convert_chunk(chunk):
    """Parse chunk from V-representation and load into np_array."""
    array = [line.strip().split()[1:] for line in chunk if line.strip().startswith('0')] # takes all lines which start with 0 and writes it splitted as elements into list (first column not taken)
    return np.array(array, dtype=float) # choose floats instead of integer


def merge_V_representation(reaction_pair_index, array, n_reactions):
    """Merges array on given column pair and takes out all rows just containing zeros after merging."""
    for index_pair in reaction_pair_index:
        index_1, index_2 = index_pair
        merged_reaction = array[:,index_1] - array[:,index_2]
        array[:,index_1] = merged_reaction
    # since all backward reations at the end:
    array = array[:,:(n_reactions)]
    # removing all rows just containing 0
    array = array[~np.all(array == 0, axis=1)]
    return array


def get_ECMs(array, smatrix):
    """
    Creation of ECMs by multiplying transposed stochiometric matrix with pre-ECM matrix (merged matrix from V-representation).
    Input array of merged V-representation and sliced stochiometric matrix.
    Outputs an array where each column shows an external metabolite and each row is an ECM.
    """
    # array with reaction rates for all exchange reactions with e * r_ex (e = one pre-ecm; r_ex = exchange reaction)
    # new stochiometric matrix with m * r_ex (m = metabolite; r_ex = exchange reaction)
    # m * r_ex transposed = r_ex * m
    # e * r_ex multiplied r_ex * m --> e * m
    ECMs = np.matmul(array, (smatrix.T * -1))
    return ECMs


def normalize_array(array):
    """Normalize given array on column with the maximum value."""
    max_list = np.amax(array, axis=1)
    return np.divide(array, max_list[:, np.newaxis])


def postprocessing(ECM_queue, count_queue, reaction_pair_index, n_reactions, smatrix, outputfile, separator, decimals, gzipped, parallel, control, chunk):
    """Receives chunks and converts them to arrrays with ECMs. Outputs ECMs in csv, txt or console if gzipped is False. Compresses ECMs and writes them to queue if gzipped is True."""
    ECM = convert_chunk(chunk)
    # if cunksize is very low, it can happen, that no ECM is in chunk (especially at start where a lot of header lines are)
    if len(ECM) == 0:
        return None
    ECM = merge_V_representation(reaction_pair_index, ECM, n_reactions)
    ECM = get_ECMs(ECM, smatrix)
    ECM = normalize_array(ECM)
    #if decimals <=4:
    #    np.ECM.astype(dtype=np.float32) # to reduce memory consumption
    # count number of ECMs
    count = np.shape(ECM)[0]
    if gzipped and not parallel:
        with gzip.open(outputfile, 'a') as output_file:
            np.savetxt(output_file, ECM, delimiter=separator, fmt=f'%1.{decimals}f')
        return count
    elif gzipped and parallel:
        ECM_queue.put(ECM)
        count_queue.put(count)
        if control:
            control.release()
    else:
        with open(outputfile, 'a') as output_file:
            np.savetxt(output_file, ECM, delimiter=separator, fmt=f'%1.{decimals}f')
        if parallel and control:
            control.release()
        if not parallel:
            return count


def count_queue_processing(count_queue, ECM_count, chunksize):
    """Receives count_queue and multiprocessing.Value object. Calculates total number of ECMs and returns it by writitng to multiprocessing.Value object."""
    # get number of all ECMs
    add = 0
    count = 0
    while True:
        while not count_queue.empty():
            add = count_queue.get()
            count += add
        if count_queue.empty():
            if add == 0:
                ti.sleep(3) # initial waiting time not depending on chunksize
            if chunksize <= 10000:
                ti.sleep(1)
            if chunksize > 10000 <= 100000:
                ti.sleep(chunksize * 0.00014)
            if chunksize > 100000:
                ti.sleep(chunksize * 0.0002)
            if count_queue.empty():
                break
    # write total number of ECMs into mp.object
    ECM_count.value = count


def ECM_queue_processing(ECM_queue, header, chunksize, outputfile):
    """Receives ECM_queue and writes it to a .gz file."""
    # write output serial from queue to .gz
    add = 0
    with gzip.open(outputfile, 'wb') as output_file:
        output_file.write((header + '\n').encode())
        output_file.close()
    with gzip.open(outputfile, 'a') as output_file:
        while True:
            while not ECM_queue.empty():
                ECM = ECM_queue.get()
                add = np.shape(ECM)[0]
                np.savetxt(output_file, ECM, delimiter=separator, fmt=f'%1.{decimals}f')
            if ECM_queue.empty():
                if add == 0:
                    ti.sleep(3) # initial waiting time not depending on chunksize
            if chunksize <= 10000:
                ti.sleep(1)
                if ECM_queue.empty():
                    break
            if chunksize > 10000 <= 100000:
                ti.sleep(chunksize * 0.00014)
                if ECM_queue.empty():
                    break
            if chunksize > 100000:
                ti.sleep(chunksize * 0.0002)
                if ECM_queue.empty():
                    break
    output_file.close()


if __name__ == '__main__':

    start = ti.time()
    
    ### argparse ###
    parser = argparse.ArgumentParser(description="Python script to calculate ECM's from given metabolic model.\n@author: Christian Mayer, Marcus Holzer, Bianca Buchner", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # create group for required arguments
    parser_req = parser.add_argument_group('required arguments')

    parser_req.add_argument('-f', '--file',
                            help='Enter input sbml-File.',
                            type=str,
                            metavar='PATH_TO_FILE',
                            action='store', 
                            required=True)
    parser_req.add_argument('-m', '--model_name',
                            help='Enter name of the model. (not filename) If no --result_name is given, model name is used to name result files.',
                            type=str,
                            metavar='STR',
                            action='store',
                            required=True)

    # optional arguments
    parser.add_argument('-n', '--n_processes',
                        help='Give number of processes mplrs will use.',
                        type=int,
                        metavar='INT',
                        action='store',
                        default=3)
    parser.add_argument('-ch', '--chunksize',
                        help='Give the size of chunks by number of lines per chunk from V-representation for postprocessing via multiprocessing. Lower chunksizes (e.g. 50000) are preferable on lower sized models.',
                        type=int,
                        metavar='INT',
                        action='store',
                        default=100000)
    parser.add_argument('-p', '--parallel',
                        help='If flag is given, postprocessing will be performed parallel with multiprocessing. This has high RAM usage on big models but postprocessing step is much faster.',
                        action='store_true')
    parser.add_argument('-gz', '--gzipped',
                        help='If flag is given, result file will be .gz file. In parallel mode, creation of a .gz file will take significantly more time.',
                        action='store_true')
    parser.add_argument('-t', '--time',
                        help='If flag is given as an option, script creates an additional json file in outpath with time measurements of different working steps. Note, that total time doesnt measure time for import statements.',
                        action='store_true')
    parser.add_argument('--pool',
                        help='If flag is given, postprocessing will be performed with mp.Pool instead of mp.Process.',
                        action='store_true')
    parser.add_argument('-op', '--only_projection',
                        help='If flag is given, only H-projection will be created in tmp directory. Stops after projection step.',
                        action='store_true')
    parser.add_argument('-fv', '--fluxvertices',
                        help='If flag is given, only V-representation will be created in tmp directory. Stops after vertex enumeration step.',
                        action='store_true')
    parser.add_argument('-po', '--only_postprocessing',
                        help='Give V-representation as input. Only postprocessing and certain steps of preprocessing which are necessary will be performed. Usefull, if comprehensive V-representation of a model already exists.',
                        type=str,
                        metavar='FILE',
                        action='store',
                        default=None)
    parser.add_argument('-mp', '--mplrs',
                        help='Path to mplrs file.',
                        type=str,
                        metavar='FILE',
                        action='store',
                        default='mplrs')
    parser.add_argument('-o', '--outpath',
                        help='Directory, where results shall be saved.',
                        type=str,
                        metavar='PATH',
                        action='store',
                        default='./')
    parser.add_argument('-rn', '--result_name',
                        help='File name of result. If this option is not given, name of result file will be outpath + model_name.',
                        type=str,
                        metavar='STR',
                        action='store',
                        default=None)
    parser.add_argument('-sep', '--separator',
                        help='Option, which separator output will have. If no --result_name is given, result file will end with .csv for [,;] and .txt for everything else.',
                        type=str,
                        metavar='STR',
                        action='store',
                        default=',')
    parser.add_argument('-d', '--decimals',
                        help='Give the number of decimals, the resulting ECMs will have.',
                        type=int,
                        metavar='INT',
                        action='store',
                        default=4)
    parser.add_argument('-ex', '--extern',
                        help='Give reaction ids as input (e.g. -ex R1 R2). Just the given reactions are marked as extern and will show up in the ECMs. Use only for real extern reactions, not internal ones! Beware, if the chosen reactions are not all extern reactions of the model, results are unbalanced.',
                        type=str,
                        nargs='+',
                        metavar='R1 R2 R3 ...',
                        action='store',
                        default=None)
    parser.add_argument('-tmp', '--tmppath',
                        help='Directory, where tmp files get stored.',
                        type=str,
                        metavar='PATH',
                        action='store',
                        default='./tmp/')
    parser.add_argument('-dv', '--developer',
                        help='Intermediate files in tmp directory get not deleted. H- and V-representation are maintained.',
                        action='store_true')
    parser.add_argument('-v', '--verbose',
                        help='If flag is given, mplrs will show the whole output.',
                        action='store_true')
    parser.add_argument('-ap', '--approximation',
                        help='Approximation for float numbers converted to rationals for .ine files. The less strict the border, the higher rationals can get. Keep attention, large numbers can lead to intense performance issues. If approximation is set to None (default), border is 1e6. Floats will get "rounded" to fractions with max this number as denominator. (e.g.: 1e3 means denominator can be max 1000 --> 1/1000: 0.001)',
                        type=str,
                        metavar='STR',
                        action='store',
                        default='1e06')
    parser.add_argument('-mf', '--mfel',
                        help='Mplrs project will be performed with mfel file if flag is given as an option.',
                        action='store_true')
    parser.add_argument('-ro', '--rows',
                        help='The number of rows per job for the mplrs algorithm.',
                        type=int,
                        metavar='INT',
                        action='store',
                        default=20)
    parser.add_argument('-lr', '--lastrows',
                        help='The number of rows for the last lastp jobs for the mplrs algorithm.',
                        type=int,
                        metavar='INT',
                        action='store',
                        default=20)
    parser.add_argument('-lp', '--lastp',
                        help='Give the percentage of processes, which get used for lastrows of mplrs algorithm.',
                        type=int,
                        metavar='INT',
                        action='store',
                        default=10)
    
    
    args = parser.parse_args()

    time_initial_setup_start = ti.time()
    
    # process id
    print(f'Process ID: {os.getpid()}')
    
    ### set names
    sbmlfile = args.file
    core_name = args.model_name
    path_mplrs = args.mplrs
    n_processes = args.n_processes
    separator = args.separator
    decimals = args.decimals
    approximation = args.approximation
    rows = args.rows
    lastp = args.lastp
    lastrows = args.lastrows
    chunksize = args.chunksize
    developer = args.developer
    verbose = args.verbose
    gzipped = args.gzipped
    only_projection = args.only_projection
    fluxvertices = args.fluxvertices
    only_postprocessing = args.only_postprocessing
    parallel = args.parallel
    extern = args.extern
    pool_switch = args.pool
    
    if pool_switch and not parallel:
        raise Exception('The --pool option can only be performed if --parallel option is enabled.')
    
    # set temporary directory
    tmp_dir = args.tmppath + '/'
    # set result directory
    outpath = args.outpath #'./'
    
    # test if outpath directory exists
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
        print(f'Created directory {outpath}.')
    
    # test if path for tmp file exists
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)
        print(f'Created directory {tmp_dir}.')
    
    # get path, where this file is saved
    dir_path = os.path.dirname(os.path.realpath(__file__))
    # mfel file always has to be in directory ./mplrs_scripts/
    mfel_file = dir_path + '/mplrs_scripts/mfel.tcsh'
    
    # set input and output file
    if not only_postprocessing:
        v_representation = tmp_dir + core_name + '.projected'
    else:
        v_representation = only_postprocessing
    
    if args.result_name:
        outputfile = args.result_name
    else:
        if (separator == ',') or (separator == ';'):
            end_of_file = 'csv'
        else:
            end_of_file = 'txt'

        if gzipped:
            outputfile = outpath + core_name + '_ecms.' + end_of_file + '.gz'
        else:
            outputfile = outpath + core_name + '_ecms.' + end_of_file

    ### preprocessing
    
    # create temporary directory
    try:
        os.mkdir(tmp_dir)
    except FileExistsError as error:
        print('Directory for temporary files already exists.')
    
    # read in model via cobra
    model, added_reas = read_model(sbmlfile)
    # remove added reactions from cobra if added
    model = rm_reactions(model, added_reas)
    print(f'Model: {sbmlfile} succesfully parsed with cobra.')
    
    # get extern compartments (also compartments with no name are taken but filtered out later)
    extern_compartments = get_extern_compartments(model)
    
    # add an exchange attribute to all reactions
    if extern:
        # if extern reactions are given manually
        mark_exchange(model, extern)
    else:
        # proof all reactions if they are exchange reactions
        indicate_exchange(model, extern_compartments)
    
    # split all reversible reactions given by new_ex_reactions
    split_extern_reversible_reactions(model)
    
    # create ex_reactions with splitted reactions
    ex_reactions = ex_reaction_id_and_index(model)
    
    # write direction of each reaction as an attribute
    write_reaction_direction(model)
    
    if not only_postprocessing:
        # create stochiometric matrix with splitted reactions
        smatrix = cobra.util.array.create_stoichiometric_matrix(model)
        # correct stochiometric matrix for backward reactions
        smatrix = correct_stochiometric_matrix(model, smatrix)

        ### create H- and V-representation
        # write H-representation to .ine file from cobra model
        write_h_representation(smatrix, model, tmp_dir, core_name, approximation=approximation)

        # perform mplrs redund on file; create .postredund file
        redund(n_processes, path_mplrs, tmp_dir, core_name, verbose=verbose)
        # rewrite postredund file
        rewrite_postredund_file(tmp_dir, core_name) # deletes last rows of .postredund files
        write_project_line(ex_reactions, tmp_dir, core_name) # writes project line into .postredund file
        print('H-representation created.')

        time_initial_setup_end = ti.time()
        print(f'Time for creation of postredund file (inputfile for projection): {round((time_initial_setup_end - time_initial_setup_start), 5)} seconds.')

        time_projection_start = ti.time()

        # projection
        if args.mfel:
            mfel_project(mfel_file, n_processes, tmp_dir, core_name, rows=rows, lastp=lastp, lastrows=lastrows, verbose=verbose)
        else:
            mplrs_project(smatrix, ex_reactions, n_processes, path_mplrs, tmp_dir, core_name, rows=rows, lastp=lastp, lastrows=lastrows, verbose=verbose) # creates _h.projected file, .postredund file disappears

        time_projection_end = ti.time()
        print(f'Projection ran in {round((time_projection_end - time_projection_start), 5)} seconds.')

        if only_projection:
            end = ti.time()
            print(f'Ran in {round(end - start, 5)} seconds.')
            if args.time:
                # export times to file
                time_total = end - start
                time_initial_setup = time_initial_setup_end - time_initial_setup_start
                time_projection = time_projection_end - time_projection_start
                times_dict = {'time_total': time_total, 'time_initial_setup': time_initial_setup, 'time_projection': time_projection}
                print(f'Export times and memory to {outpath + core_name + "_times.json"}.')
                json.dump(times_dict, open(outpath + core_name + '_times.json','w'))

            exit()

        # enumeration
        time_enumeration_start = ti.time()
        mplrs_conversion(n_processes, path_mplrs, tmp_dir, core_name, verbose=verbose) # converts H to V-representation, .projected file is V-representation
        time_enumeration_end = ti.time()
        print(f'Conversion ran in {round((time_enumeration_end - time_enumeration_start), 5)} seconds.')  

        print('V-representation created.')

        if fluxvertices:
            end = ti.time()
            print(f'Ran in {round(end - start, 5)} seconds.')
            if args.time:
                # export times to file
                time_total = end - start
                time_initial_setup = time_initial_setup_end - time_initial_setup_start
                time_projection = time_projection_end - time_projection_start
                time_enumeration = time_enumeration_end - time_enumeration_start
                times_dict = {'time_total': time_total, 'time_initial_setup': time_initial_setup, 'time_projection': time_projection, 'time_enumeration': time_enumeration}
                print(f'Export times and memory to {outpath + core_name + "_times.json"}.')
                json.dump(times_dict, open(outpath + core_name + '_times.json','w'))
            exit()
    
    time_postprocessing_start = ti.time()
    
    ## prepare stochiometric matrix
    # merge forward and backward reaction
    reaction_pair_index = merge_model(model)
    # create stochiometric matrix again with merged reactions (needed for calculation of ECMs)
    # stochiometric matrix with m * r
    smatrix = cobra.util.array.create_stoichiometric_matrix(model)
    # correct stochiometric matrix for irreversible backward reactions
    smatrix = correct_stochiometric_matrix(model, smatrix)
    # create ex_reactions again with merged reactions
    ex_reactions = ex_reaction_id_and_index(model)
    # get all metabolites attached to external reactions and their indices
    ex_metabolites = ex_metabolites_id_and_index(model, extern_compartments)
    # slice stochiometric matrix to just get metabolites attached to extern reactions and extern reactions
    smatrix = slice_stochio_matrix(smatrix, ex_reactions, ex_metabolites)
    
    ### create ECMs from V-representation

    # get number of extern reactions after merging the model
    n_reactions = len(ex_reactions)
    
    if parallel:
        print(f'Performing postprocessing in multiprocessing mode.')
    
    ### output ECMs
    print(f'Reading in V-representation and output results.')
    ## header
    # create header for result
    header = write_output_header(ex_metabolites, separator=separator)
    # write header if not gzipped
    # since all output for .gz files have to derive from the same process, gzipped header is not written here (lookup in queue_processing function)
    if not gzipped:
        output_file = open(outputfile, "w+")
        output_file.write(header + '\n')
        output_file.close()
    if gzipped and not parallel:
        output_file = gzip.open(outputfile, 'wb')
        output_file.write((header + '\n').encode())
        output_file.close()
  
    ## convert pre_ECM chunks to ECM chunks and write them to output via multiprocessing or single processing
    
    if pool_switch:
        control = False
    else:
        # defining semaphore to control maximal number of parllel processes
        control = mp.Semaphore(n_processes)
    
    
    with mp.Manager() as manager:
        # write queue for ECMs and counts
        count_queue = manager.Queue() # stays empty if not gzipped
        ECM_queue = manager.Queue() # stays empty if not gzipped
        
        # setting function object with arguments
        postprocessing_part = partial(postprocessing, ECM_queue, count_queue, reaction_pair_index, n_reactions, smatrix, outputfile, separator, decimals, gzipped, parallel, control)
        
        # set up name for total number of ECMs
        if parallel and gzipped:
            # write mp value for ECM_count
            ECM_count = manager.Value('i', 0)
        if not parallel:
            ECM_count = 0
        
        with mp.Pool(n_processes) as pool:
            # open v_representation
            with open(v_representation, "r") as inputfile:
                n = 0
                while True:
                    n += 1
                    # take n_lines (chunksize) from file generator object, put them into a list (array)
                    v_chunk = list(islice(inputfile, chunksize))
                    
                    if not v_chunk:               
                        if parallel:
                            print(f'Total number of executed processes for postprocessing: {n-1}')
                            if pool_switch: # the more pythonic and "stable" approach but needs a lot of RAM on big vertex enumeration files
                                pool.close()
                                pool.join()
                            else:
                                if chunksize < 50000:
                                    ti.sleep(0.5 + (chunksize * 0.00001)) # sleep is needed to get all ECMs of the last chunks since last process sometimes finishes earlier than processes before
                                else:
                                    ti.sleep(chunksize * 0.00005)
                                process.join()
                                
                            if gzipped:
                                count_queue_process.join()
                                ECM_queue_process.join()
                                ECM_count = ECM_count.value # get total number of ECMs
                        break
                    
                    if parallel and not pool_switch:
                        control.acquire()
                        # process postprocessing with cores from pool
                        process = mp.Process(target=postprocessing, args=[ECM_queue, count_queue, reaction_pair_index, n_reactions, smatrix, outputfile, separator, decimals, gzipped, parallel, control, v_chunk])
                        process.start() # start the child process with chunk
                    elif parallel and pool_switch:
                        pool.apply_async(postprocessing_part, (v_chunk, ))
                    else:
                        count = postprocessing(ECM_queue, count_queue, reaction_pair_index, n_reactions, smatrix, outputfile, separator, decimals, gzipped, parallel, control, v_chunk)
                        ECM_count += count

                    if n == 1 and parallel and gzipped:
                        # add count output
                        count_queue_process = mp.Process(target=count_queue_processing, args=[count_queue, ECM_count, chunksize])
                        count_queue_process.start()
                        # write ECM output to gzip file
                        # has to be done in parallel since queue only has a certain number of entries it can hold at a time
                        # it is inportant, that only the same process writes to .gz file --> so it got serialized here
                        ECM_queue_process = mp.Process(target=ECM_queue_processing, args=[ECM_queue, header, chunksize, outputfile])
                        ECM_queue_process.start()
    
    # calculate number of ECMs in parallel mode
    # got seperated from count_queue, since this queue does not give reliable ECM counts all the time
    if parallel and not gzipped: # gzipped nevertheless works with this queue since we cant ask wc -l on gzip file to get unzipped number of lines
        cmd = ["wc", "-l", outputfile]
        ECM_count_process = subprocess.run(cmd, capture_output=True, text=True)
        ECM_count = ECM_count_process.stdout
        ECM_count = int(ECM_count.split()[0]) -1
    
    
    ## remove tmp directory or just files created there
    if not developer:
        # get all paths of files of tmp directory as list
        all_tmp_filenames = [list(os.walk(tmp_dir))[0][0] + file for file in list(os.walk(tmp_dir))[0][2]]
        # delete all temporary files created
        for file in all_tmp_filenames:
            if core_name in file:
                os.remove(file)

        # delete directory
        try:
            os.rmdir(tmp_dir)
        except OSError as error:
            print('Directory not empty. Only deleted temporary files created in this process and not whole directory.')
    
    
    time_postprocessing_end = ti.time()
    print(f'Time for creation of result (after conversion): {round((time_postprocessing_end - time_postprocessing_start), 5)} seconds.')
    
    print(f'{ECM_count} ECMs are found.')
    
    end = ti.time()
    print(f'Ran in {round(end - start, 5)} seconds.')
    
    if args.time:       
        # export times to file
        time_total = end - start
        if not only_postprocessing:
            time_initial_setup = time_initial_setup_end - time_initial_setup_start
            time_projection = time_projection_end - time_projection_start
            time_enumeration = time_enumeration_end - time_enumeration_start
        else:
            time_initial_setup = 0
            time_projection = 0
            time_enumeration = 0
        time_postprocessing = time_postprocessing_end - time_postprocessing_start
        times_dict = {'time_total': time_total, 'time_initial_setup': time_initial_setup, 'time_projection': time_projection, 'time_enumeration': time_enumeration, 'time_postprocessing': time_postprocessing, 'ECM_count': ECM_count}
        print(f'Export times and memory to {outpath + core_name + "_times.json"}')
        json.dump(times_dict, open(outpath + core_name + '_times.json','w'))           
