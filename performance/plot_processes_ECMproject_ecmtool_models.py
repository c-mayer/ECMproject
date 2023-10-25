#!/usr/bin/env python

"""Python script to plot mean times of model list.

@author: Christian Mayer
"""

# https://matplotlib.org/stable/gallery/lines_bars_and_markers/masked_demo.html
# https://stackoverflow.com/questions/36455083/how-to-plot-and-work-with-nan-values-in-matplotlib


### import statements
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import container
import statistics as stat
import json
from time import time
import os


### functions

def load_time_results(filename):
    '''Loads in dictionary of time result file.'''
    return json.load(open(filename))
    

def get_plotting_input(dictionary, lookup_time, json_file):
    core_list = [int(cores) for cores in dictionary.keys()]
    mean_list = [dictionary[str(cores)][lookup_time]['mean'] for cores in core_list]
    try:
        stdev_list = [dictionary[str(cores)][lookup_time]['stdev'] for cores in core_list]
    except KeyError:
        stdev_list = len(core_list)*[0]
        print(f'No stdevs for {json_file}.')
    return (core_list, mean_list, stdev_list)
    

#result = create_plotting_input('time_total', n_processes, own_dictionary=own_dictionary)
#print(result)


if __name__ == '__main__':
    
    # set up architecture for further performance measurement
    dir_path = os.path.dirname(os.path.realpath(__file__))
    outpath = dir_path + '/figures/processes/'
    # create figures directory
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    
    inpath_ECMproject = dir_path + '/times/ECMproject/process/'
    inpath_ecmtool = dir_path + '/times/ecmtool/'
    
    n_processes = [3, 5, 10, 20, 30, 40, 50, 60]
    lowest_n_processes = min(n_processes)
    lookup_time = 'time_total'
    log_scale = True
    memory = True
    ecmtool = True
    stdev = False
    
    if memory:
        if log_scale:
            ax_label = 'memory [kB]'
        else:
            ax_label = 'memory [GB]'
        lookup_time = 'max_rss'
    else:
        ax_label = 'wall clock time [s]'
    
    ### plot total times of all mmsyn models
    x = [0, 5, 10, 15]
    models = ['mmsyn_sm' + str(num).zfill(2) for num in x]
    
    fig, ax = plt.subplots()
    
    for model in models:
        if ecmtool:
            json_file = inpath_ecmtool + model + '_ecmtool_result_times.json'
        else:
            json_file = inpath_ECMproject + model + '_mplrs_project_result_times.json'
        dictionary = load_time_results(json_file)
        cores, means, stdevs = get_plotting_input(dictionary, lookup_time, json_file)
        
        # plotting
        if stdev:
            plt.errorbar(cores, means, stdevs, label=model, marker='o', markersize=4, capsize=3)
        else:
            plt.plot(cores, means, label=model, marker='o', markersize=4)

    plt.margins(x=0.02, tight=True)
    ax.tick_params(axis='both', direction='in', which='both')
    
    # x axis
    plt.xlabel('number processes')
    plt.xticks(n_processes, n_processes)
    
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
        if not memory:
            y_max = 11000
            plt.ylim(0,y_max)
            interval = y_max/100
            plt.text(lowest_n_processes, 60 + interval, '1 minute')
            plt.text(lowest_n_processes, 3600 + interval, '1 hour')
            #plt.text(lowest_n_processes, 86400 + interval, '1 day')
        else:
            y_max = 1000000
            plt.ylim(0,y_max)
            interval = y_max/100
            #plt.text(lowest_n_processes, 1000 + interval, '1 MB')
            #plt.text(lowest_n_processes, 1000000 + interval, '1 GB')
            #plt.text(lowest_n_processes, 1000000000 + interval, '1 TB')
    
    # grid
    if not memory:
        ax.axhline(60, linestyle='--', color='grey', linewidth=0.5)
        ax.axhline(3600, linestyle='--', color='grey', linewidth=0.5)
        #ax.axhline(86400, linestyle='--', color='grey', linewidth=0.5)
    else:
        if log_scale:
            ax.axhline(1000, linestyle='--', color='grey', linewidth=0.5)
            ax.axhline(1000000, linestyle='--', color='grey', linewidth=0.5)
            #ax.axhline(1000000000, linestyle='--', color='grey', linewidth=0.5)
        else:
            plt.grid(visible=True, which='major', axis='y', color='grey', linestyle='--', linewidth=0.5)
    
    # legend
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles] # take out errorbars of legend objects if object stems from plt.errorbar
    if  (len(models) % 2) == 0:
        ncols = 2
    elif len(models) <= 1:
        ncols = 1
    else:
        ncols = 3
    plt.legend(handles, labels, loc="center", bbox_to_anchor=(0.5, -0.18), ncols=ncols)
    
    # save file
    if ecmtool:
        pre = 'ecmtool'
    else:
        pre = 'ECMproject'
    
    if log_scale:
        plt.savefig(outpath + pre + '_models_log_' + lookup_time + '.png', bbox_inches='tight') # called function does not close plt object --> we can manipulate it outside the function until we close it
    else:
        plt.savefig(outpath + pre + '_models_' + lookup_time + '.png', bbox_inches='tight') # called function does not close plt object --> we can manipulate it outside the function until we close it
    plt.close() # its important to close it before we start the next plot_time call
