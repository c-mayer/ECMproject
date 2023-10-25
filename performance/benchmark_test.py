#!/usr/bin/env python

"""Python script to benchmark clock time of own implementation with clock time of the predecessor ecmtool and mfel.

Requirements to work when working with all three tool options: Needs environment with following programs:
                      1) installed tsch for mfel scrip (conda install -c conda-forge tcsh)
                      2) installed ecmtool (git clone https://github.com/SystemsBioinformatics/ecmtool.git)
                         python 3.8
                         ecmtool requirements into ecmtool directory (pip install -r requirements.txt --upgrade -t .)
                         ecmtool lib (conda install -c conda-forge lrslib=70.a)
                      3) cobra (conda install -c conda-forge cobra)
                         and my own code

@author: Christian Mayer
"""

# activate conda channel "all" to execute
import subprocess
from time import time
import json
from collections import defaultdict
import statistics as stat
import argparse
import os


def get_json_times(time_file, sbmlfile, own_tool, ecmtool, gnutime, path2mplrs, tmp_path, core_name, script, n_cores, n_repeats, chunksize):
    '''Runs code multiple times with n_cores and N-repeats and creates a dictionary which stores all times in seconds gathered from json files in dictionary with cores as keys. Each key itself opens a new dictionary with different kind of times as keys and list of repeats as values.
    Structure of dictionary: dict = {cores: {times: {'all_times': []}}}'''
    result_dict = defaultdict(dict)
    for cores in n_cores:
        print(f'Run with {cores} cores.')
        for i in range(n_repeats):
            print(f'Repeat {i + 1}')
            # run script
            if script == 'own':
                cmd = [gnutime, '-v', own_tool, '-f', sbmlfile, '-m', core_name, '-n', str(cores), '-o', tmp_path, '-mp', path2mplrs, '--time', '--chunksize', str(chunksize)] # default cores are 3 and chunksize 100 000
            elif script == 'ecmtool':
                cmd = [gnutime, '-v', 'python', ecmtool, '--model_path', sbmlfile, '--out_path', tmp_path + core_name + "_ecmtool.csv", '--path2mplrs', path2mplrs, '--processes', str(cores)] # default cores are 3
            elif script == 'mfel':
                cmd = [gnutime, '-v', own_tool, '-f', sbmlfile, '-m', core_name, '-n', str(cores), '-o', tmp_path, '-mp', path2mplrs, '--time', '--mfel', '--chunksize', str(chunksize), '-p'] # default cores are 3 and chunksize 100 000
            #output = subprocess.run(cmd, stderr=subprocess.PIPE) # capture_output=False leads to output written to console
            #print(output)
            
            process = subprocess.Popen(cmd, stderr=subprocess.PIPE)
            process.wait()

            output = process.stderr.read().decode() # output of gnutime goes into stderr
            output = output.split('\n')
            
            # get max rss and swap from process
            for line in output:
                if line.strip().startswith('Maximum resident set size'):
                    rss = int(line.split(':')[1])
                if line.strip().startswith('Swaps'):
                    swaps = int(line.split(':')[1])
            
            # get times dictionary from json file
            dictionary = json.load(open(time_file))
            if i == 0: # first iteration
                # all keys of input dict
                for key in dictionary.keys():
                    result_dict[cores][key] = {'all_values': []}
                # memory keys
                result_dict[cores]['max_rss'] = {'all_values': []}
                result_dict[cores]['swap'] = {'all_values': []}

            # all values of input dict
            for key, value in dictionary.items():
                result_dict[cores][key]['all_values'].append(value)
            # memory values
            result_dict[cores]['max_rss']['all_values'].append(rss)
            result_dict[cores]['swap']['all_values'].append(swaps)
            
    return result_dict


def get_mean_stdev(result_dict, stdev=True):
    '''Modifies a dictionary created by the function "get_json_times" with all mean times and standard deviation in seconds.
    Structure of dictionary: dict = {cores: {times: {mean: int, stdev: int}}}'''
    for key_1 in result_dict.keys(): # first key for cores
        for key_2 in result_dict[key_1].keys(): # second key for kind of time
            time_list = result_dict[key_1][key_2]['all_values'] # third key for all_times, extract list from there
            # calculate mean
            mean_time = stat.mean(time_list)
            # add mean
            result_dict[key_1][key_2]['mean'] = mean_time
            if stdev:
                # calculate standard deviation
                stdev_time = stat.stdev(time_list)
                # add stdev
                result_dict[key_1][key_2]['stdev'] = stdev_time
            
    return result_dict


def total_time_means_to_csv(result_dict, n_cores, filename, core_name):
    'Takes total time means of result dictionary and writes them into csv file. Just can do it for one model at a time.'
    if not os.path.exists(filename):
        # create header
        header = ['models'] + [str(cores) for cores in n_cores]
        with open(filename , 'w') as file:
            file.write(','.join(header) + '\n')
    with open(filename , 'a') as file:
        line = [core_name] + [str(result_dict[cores]['time_total']['mean']) for cores in n_cores]
        file.write(','.join(line) + '\n')


def all_means_to_csv(result_dict, n_cores, filename):
    'Creates a csv file with all means for the model in progress.'
    columns = ['cores']
    index = 0
    with open(filename, 'w') as file:
        for cores in n_cores:
            row = [cores]
            for key in result_dict[cores]:
                if index == 0:
                    columns.append(str(key))
                val = result_dict[cores][key]['mean']
                row.append(str(val))
            if index == 0:
                file.write(','.join(columns) + '\n')
            file.write(','.join(row) + '\n')
            index += 1
        

if __name__ == '__main__':
    
    ### argparse ###
    parser = argparse.ArgumentParser(description='Python script to benchmark clock time of own implementation with clock time of the predecessor ecmtool and mfel. For requirements to work see docstring of file. @author: Christian Mayer', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # create group for required arguments
    parser_req = parser.add_argument_group('required arguments')

    parser_req.add_argument('-f', '--file',
                            help='Enter input sbml-File.',
                            type=str,
                            metavar='PATH_TO_FILE',
                            action='store', 
                            required=True)
    parser_req.add_argument('-m', '--model_name',
                            help='Enter name of the model. (not filename)',
                            type=str,
                            metavar='STR',
                            action='store',
                            required=True)
    parser_req.add_argument('-t', '--tools',
                            help='Give list of tools which will be used. Three entries are posible: "own", "ecmtool", "mfel". At least one of these entries must be given. e.g. if you want own and ecmtool: own,ecmtool',
                            type=str,
                            metavar='STR',
                            action='store',
                            required=True)

    # optional arguments
    
    parser.add_argument('-n', '--n_cores',
                        help='Give list of numbers of cores mplrs will use. e.g. if you want 10 and 40 cores: 10,40',
                        type=str,
                        metavar='STR',
                        action='store',
                        default='10,20,40,60')
    parser.add_argument('-r', '--repeats',
                        help='Give the number of repeats.',
                        type=int,
                        metavar='INT',
                        action='store',
                        default=3)
    parser.add_argument('-ch', '--chunksize',
                        help='Give the size of chunks by number of lines per chunk from V-representation for postprocessing via multiprocessing.',
                        type=int,
                        metavar='INT',
                        action='store',
                        default=100000)
    parser.add_argument('-ns', '--no_stdev',
                        help='If given, no standard deviation will be calculated. Useful if repeat is set to one.',
                        action='store_true')
    parser.add_argument('-op', '--outpath',
                        help='Directory, where output shall be saved.',
                        type=str,
                        metavar='PATH',
                        action='store',
                        default='./')
    parser.add_argument('-mp', '--mplrs',
                        help='Path to mplrs file.',
                        type=str,
                        metavar='FILE',
                        action='store',
                        default='/opt/lrslib/v072/mplrs')
    parser.add_argument('-ot', '--own_tool',
                        help='Path to my own new mplrs project execution file.',
                        type=str,
                        metavar='FILE',
                        action='store',
                        default='../ECMproject.py')
    parser.add_argument('-et', '--ecmtool',
                        help='Path to ecmtool main_modified.py file.',
                        type=str,
                        metavar='FILE',
                        action='store',
                        default='../../../ecmtool_13_07_modified/main.py')
    parser.add_argument('-gt', '--gnutime',
                        help='Path to "time" file of gnutime for memory measurement.',
                        type=str,
                        metavar='FILE',
                        action='store',
                        default='./gtime/bin/time')
    parser.add_argument('-v', '--verbose',
                        help='Result dictionary will be printed also as console output.',
                        action='store_true')
    
    args = parser.parse_args()
    
    start = time()
    
    # set variables
    sbmlfile = args.file
    core_name = args.model_name
    tools = args.tools.split(',')
    n_cores = args.n_cores.split(',')
    n_repeats = args.repeats
    own_tool = args.own_tool
    ecmtool = args.ecmtool
    path2mplrs = args.mplrs
    outpath = args.outpath
    chunksize = args.chunksize
    gnutime = args.gnutime

    
    if n_repeats < 1:
        raise Exception('Given number of repeats has to be at least 1.')
    
    if args.no_stdev or n_repeats <= 1:
        stdev = False
    else:
        stdev = True
    
    # proof which tools are given
    print(f'Choosen tools: {tools}')
    if 'own' in tools:
        own_switch = True
    else:
        own_switch = False
    
    if 'ecmtool' in tools:
        ecmtool_switch = True
    else:
        ecmtool_switch = False
    
    if 'mfel' in tools:
        mfel_switch = True
    else:
        mfel_switch = False
    
    print(f'Choosen cores: {n_cores}')
    
    # set up architecture for further performance measurement
    dir_path = os.path.dirname(os.path.realpath(__file__))
    tmp_path = dir_path + '/results/'
    # create temporary directory
    if not os.path.exists(tmp_path):
        os.mkdir(tmp_path)
    
    if own_switch:
        ##### own code #####
        print('=== Test for own code. ===')
        script = 'own'
        # json file
        time_file = './results/' + core_name + '_times.json'
        result_own = get_json_times(time_file, sbmlfile, own_tool, ecmtool, gnutime, path2mplrs, tmp_path, core_name, script, n_cores, n_repeats, chunksize)
        result_own = get_mean_stdev(result_own, stdev=stdev)
        json.dump(result_own, open(outpath + core_name + '_mplrs_project_result_times.json','w'))
        print(f'{outpath + core_name + "_mplrs_project_result_times.json"} created.')
        # csv files
        csv_file = outpath + '/' + 'all_total_times_own.csv'
        total_time_means_to_csv(result_own, n_cores, csv_file, core_name)
        print(f'{csv_file} created')
        csv_file = outpath + '/' + core_name + '_own.csv'
        all_means_to_csv(result_own, n_cores, csv_file)
        print(f'{csv_file} created')
        
        if args.verbose:
            print('=== Result own code. ===')
            print(result_own)

        ## times ##
        #total
        #initial_setup
        #projection
        #enumeration
        #postprocessing


    if ecmtool_switch:
        ##### ecmtool #####
        print('=== Test for ecmtool. ===')
        script = 'ecmtool'
        # json file
        time_file = './results/' + core_name + '_ecmtool_times.json'
        result_ecmtool = get_json_times(time_file, sbmlfile, own_tool, ecmtool, gnutime, path2mplrs, tmp_path, core_name, script, n_cores, n_repeats, chunksize)
        result_ecmtool = get_mean_stdev(result_ecmtool, stdev=stdev)
        json.dump(result_ecmtool, open(outpath + core_name + '_ecmtool_result_times.json','w'))
        print(f'{outpath + core_name + "_ecmtool_result_times.json"} created.')
        # csv files
        csv_file = outpath + '/' + 'all_total_times_ecmtool.csv'
        total_time_means_to_csv(result_ecmtool, n_cores, csv_file, core_name)
        print(f'{csv_file} created')
        csv_file = outpath + '/' + core_name + '_ecmtool.csv'
        all_means_to_csv(result_ecmtool, n_cores, csv_file)
        print(f'{csv_file} created')
        
        if args.verbose:
            print('=== Result ecmtool. ===')
            print(result_ecmtool)

        ## times ##
        #total
        #preprocessing
        #first vertex enumeration
        #intermediate processing
        #second vertex enumeration
        #postprocessing


    if mfel_switch:
        ###### own code with mfel script in tcsh #####
        print('=== Test for mfel. ===')
        script = 'mfel'
        # json file
        time_file = './results/' + core_name + '_times.json'
        result_mfel = get_json_times(time_file, sbmlfile, own_tool, ecmtool, gnutime, path2mplrs, tmp_path, core_name, script, n_cores, n_repeats, chunksize)
        result_mfel = get_mean_stdev(result_mfel, stdev=stdev)
        json.dump(result_mfel, open(outpath + core_name + '_mfel_result_times.json','w'))
        print(f'{outpath + core_name + "_mfel_result_times.json"} created.')
        # csv files
        csv_file = outpath + '/' + 'all_total_times_mfel.csv'
        total_time_means_to_csv(result_mfel, n_cores, csv_file, core_name)
        print(f'{csv_file} created')
        csv_file = outpath + '/' + core_name + '_mfel.csv'
        all_means_to_csv(result_mfel, n_cores, csv_file)
        print(f'{csv_file} created')
        
        if args.verbose:
            print('=== Result mfel. ===')
            print(result_mfel)

        ## times ##
        #total
        #initial_setup
        #projection
        #enumeration
        #postprocessing
    
    # remove temporary directory
    cmd = ['rm', '-rf', tmp_path]
    subprocess.run(cmd, capture_output=False)
    
    end = time()
    print(f'Ran in {round(end - start, 5)} seconds.')
