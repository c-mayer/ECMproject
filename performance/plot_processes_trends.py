#!/usr/bin/env python

"""Python script to calculate sum of times and plot times of different models for ECMproject and optionally for ecmtool.

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
    log_scale = False
    ecmtool = True
    memory = False
    time = 'time_total'
    core_name = 'mmsyn_sm15'
    n_processes = [10, 20, 30, 40, 50, 60]
    lowest_n_processes = min(n_processes)
    stdev = False
    
    # Input File
    own_tool_json = inpath_ECMproject + '/' + core_name + '_mplrs_project_result_times.json'
    if ecmtool:
        ecmtool_json = inpath_ecmtool + '/' + core_name + '_ecmtool_result_times.json'
    
    if memory:
        lookup_time = 'max_rss'
        if log_scale:
            ax_label = 'memory [kB]'
        else:
            ax_label = 'memory [GB]'
    else:
        lookup_time = time
        ax_label = 'wall clock time [s]'
    
    ### loading in results
    # ECMproject
    dictionary_own = load_time_results(own_tool_json)
    cores_own, means_own, stdevs_own = get_plotting_input(dictionary_own, lookup_time)
    print('n_processes ECMproject')
    print(cores_own)
    
    if ecmtool:
    # ecmtool
        dictionary_ecmtool = load_time_results(ecmtool_json)
        cores_ecmtool, means_ecmtool, stdevs_ecmtool = get_plotting_input(dictionary_ecmtool, lookup_time)
        print('n_processes ecmtool')
        print(cores_ecmtool)
    
    # plotting results
    fig, ax = plt.subplots()
    #plt.title(core_name)

    if stdev:
        plt.errorbar(cores_own, means_own, stdevs_own, label='ECMproject', linestyle='-', markersize=4, marker='o', capsize=3)
    else:
        plt.plot(cores_own, means_own, label='ECMproject', linestyle='-', markersize=4, marker='o')
    if ecmtool:
        if stdev:
            plt.errorbar(cores_ecmtool, means_ecmtool, stdevs_ecmtool, label='ecmtool mplrs', linestyle='-', markersize=4, marker='s', capsize=3)
        else:
            plt.plot(cores_ecmtool, means_ecmtool, label='ecmtool mplrs', linestyle='-', markersize=4, marker='s')
    
    # legend
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles] # take out errorbars of legend objects if object stems from plt.errorbar
    if ecmtool:
        ncols = 2
    else:
        ncols = 1
    plt.legend(handles, labels, loc="center", bbox_to_anchor=(0.5, -0.15), ncols=ncols)
    
    plt.margins(x=0.03, tight=True)
    ax.tick_params(axis='both', direction='in', which='both')
    
    # x-axis
    plt.xlabel('number of processes')
    plt.xticks(n_processes, n_processes)
    
    # y-axis
    #plt.ylabel(ax_label)
    #if log_scale:
    #    plt.yscale("log")
    
    # y axis
    plt.ylabel(ax_label)
    
    if log_scale:
        plt.yscale('log')
        if not memory:
            plt.text(lowest_n_processes, 60 + 5, '1 minute')
            plt.text(lowest_n_processes, 3600 + 200, '1 hour')
            #plt.text(lowest_n_processes, 86400 + 2000, '1 day')
        else:
            plt.text(lowest_n_processes, 1000 + 100, '1 MB')
            plt.text(lowest_n_processes, 1000000 + 100000, '1 GB')
            #plt.text(lowest_n_processes, 1000000000 + 100000000, '1 TB')
    else:
        plt.yscale('linear')
        if not memory:
            y_max = 45000
            plt.ylim(0,y_max)
            interval = y_max/100
            #plt.text(lowest_n_processes, 60 + interval, '1 minute')
            plt.text(lowest_n_processes, 3600 + interval, '1 hour')
            plt.text(lowest_n_processes, 43200 + interval, '1/2 day')
        else:
            plt.yticks(range(0, 1000000 + 1, 100000))
            y_max = 1050000
            plt.ylim(0,y_max)
            interval = y_max/100
            #plt.text(lowest_n_processes, 1000 + interval, '1 MB')
            #plt.text(lowest_n_processes, 1000000 + interval, '1 GB')
            #plt.text(lowest_n_processes, 1000000000 + interval, '1 TB')
    
    # grid
    if log_scale:
        if not memory:
            ax.axhline(60, linestyle='--', color='grey', linewidth=0.5)
            ax.axhline(3600, linestyle='--', color='grey', linewidth=0.5)
            #ax.axhline(86400, linestyle='--', color='grey', linewidth=0.5)
        else:
            ax.axhline(1000, linestyle='--', color='grey', linewidth=0.5)
            ax.axhline(1000000, linestyle='--', color='grey', linewidth=0.5)
            #ax.axhline(1000000000, linestyle='--', color='grey', linewidth=0.5)
    else:
        if not memory:
            #ax.axhline(60, linestyle='--', color='grey', linewidth=0.5)
            ax.axhline(3600, linestyle='--', color='grey', linewidth=0.5)
            ax.axhline(43200, linestyle='--', color='grey', linewidth=0.5)
        else:
            plt.grid(visible=True, which='major', axis='y', color='grey', linestyle='--', linewidth=0.5)
    
    # save and close file
    if log_scale:
        plt.savefig(outpath + core_name + '_log_' + lookup_time + '_processes' + '.png', bbox_inches='tight')
    else:
        plt.savefig(outpath + core_name + '_' + lookup_time + '_processes' + '.png', bbox_inches='tight')
    plt.close()
    