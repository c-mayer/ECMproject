#!/usr/bin/env python

"""Python script to plot times of all sections of the tools ECMproject and optionally ecmtool.

@author: Christian Mayer
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import container
import statistics as stat
import json
from time import time
import os

def load_time_results(filename):
    '''Loads in dictionary of time result file.'''
    return json.load(open(filename))

def get_plotting_input(dictionary, lookup_time):
    core_list = [int(cores) for cores in dictionary.keys()]
    mean_list = [dictionary[str(cores)][lookup_time]['mean'] for cores in core_list]
    try:
        stdev_list = [dictionary[str(cores)][lookup_time]['stdev'] for cores in core_list]
    except KeyError:
        stdev_list = [0 for cores in core_list]
    return (core_list, mean_list, stdev_list)


if __name__ == '__main__':
    
    ### set variables ###
    
    # set up architecture for further performance measurement
    dir_path = os.path.dirname(os.path.realpath(__file__))
    outpath = dir_path + '/figures/processes/'
    # create figures directory
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    
    inpath_ECMproject = dir_path + '/times/ECMproject/process/'
    inpath_ecmtool = dir_path + '/times/ecmtool/'
    
    # Switches
    log_scale = True
    ecmtool = False
    times_ECMproject =  ['time_total', 'time_initial_setup', 'time_projection', 'time_enumeration', 'time_postprocessing']
    times_ecmtool = ['time_total', 'time_preprocessing', 'time_first_vertex_enumeration', 'time_intermediate_processing', 'time_second_vertex_enumeration', 'time_postprocessing']
    core_name = 'mmsyn_sm15'
    n_processes = [10, 20, 30, 40, 50, 60]
    lowest_n_processes = min(n_processes)
    stdev = False
    markers = ['o', 's', 'x', '^', 'v', '>', '<']
    
    # Input File
    ECMproject_json = inpath_ECMproject + '/' + core_name + '_mplrs_project_result_times.json'
    if ecmtool:
        ecmtool_json = inpath_ecmtool + '/' + core_name + '_ecmtool_result_times.json'
    
    ### create plot
    fig, ax = plt.subplots()   
    
    # loading in results and plotting
    if not ecmtool:
        #plt.title('ECMproject ' + core_name)
        # ECMproject
        dictionary_ECMproject = load_time_results(ECMproject_json)
        for index, lookup_time in enumerate(times_ECMproject):
            cores_ECMproject, means_ECMproject, stdevs_ECMproject = get_plotting_input(dictionary_ECMproject, lookup_time)
            
            if stdev:
                plt.errorbar(cores_ECMproject, means_ECMproject, stdevs_ECMproject, label=' '.join(lookup_time.split('_')[1:]), linestyle='-', markersize=4, marker=markers[index], capsize=3)
            else:
                plt.plot(cores_ECMproject, means_ECMproject, label=' '.join(lookup_time.split('_')[1:]), linestyle='-', markersize=4, marker=markers[index])
        
        
    else:
        #plt.title('ecmtool ' + core_name)
        # ecmtool
        dictionary_ecmtool = load_time_results(ecmtool_json)
        for index, lookup_time in enumerate(times_ecmtool):
            cores_ecmtool, means_ecmtool, stdevs_ecmtool = get_plotting_input(dictionary_ecmtool, lookup_time)
            
            if stdev:
                plt.errorbar(cores_ecmtool, means_ecmtool, stdevs_ecmtool, label=' '.join(lookup_time.split('_')[1:]), linestyle='-', markersize=4, marker=markers[index], capsize=3)
            else:
                plt.plot(cores_ecmtool, means_ecmtool, label=' '.join(lookup_time.split('_')[1:]), linestyle='-', markersize=4, marker=markers[index])
    
    
    # legend
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles] # take out errorbars of legend objects if object stems from plt.errorbar
    plt.legend(handles, labels, loc="center", bbox_to_anchor=(0.5, -0.18), ncols=3)
    
    # axes
    plt.margins(x=0.03, tight=True)
    ax.tick_params(axis='both', direction='in', which='both')
    
    # x-axis
    plt.xlabel('number of processes')
    plt.xticks(n_processes, n_processes)
    
    # y axis
    ax_label = 'wall clock time [s]'
    plt.ylabel(ax_label)
    
    if log_scale:
        plt.yscale('log')
        plt.text(lowest_n_processes, 60 + 5, '1 minute')
        plt.text(lowest_n_processes, 3600 + 200, '1 hour')
        #plt.text(lowest_n_processes, 86400 + 2000, '1 day')
    else:
        plt.yscale('linear')
        y_max = 45000
        plt.ylim(0,y_max)
        interval = y_max/100
        #plt.text(lowest_n_processes, 60 + interval, '1 minute')
        plt.text(lowest_n_processes, 3600 + interval, '1 hour')
        plt.text(lowest_n_processes, 43200 + interval, '1/2 day')
    
    # grid
    if log_scale:
        ax.axhline(60, linestyle='--', color='grey', linewidth=0.5)
        ax.axhline(3600, linestyle='--', color='grey', linewidth=0.5)
        #ax.axhline(86400, linestyle='--', color='grey', linewidth=0.5)
    else:
        #ax.axhline(60, linestyle='--', color='grey', linewidth=0.5)
        ax.axhline(3600, linestyle='--', color='grey', linewidth=0.5)
        ax.axhline(43200, linestyle='--', color='grey', linewidth=0.5)
    
    # save and close file
    if not ecmtool:
        tool = 'ECMproject'
    else:
        tool = 'ecmtool'
    
    if log_scale:
        plt.savefig(outpath + tool + '_' + core_name + '_log_time_sections_processes' + '.png', bbox_inches='tight')
    else:
        plt.savefig(outpath + tool + '_' + core_name + '_time_sections_processes' + '.png', bbox_inches='tight')
    plt.close()
    