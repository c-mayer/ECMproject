#!/usr/bin/env python

"""Python script to plot times of each section of ECMproject with mfel against ECMproject without mfel.

@author: Christian Mayer
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import container
import statistics as stat
import json
from time import time


def plot_time(model_name, n_processes, core_number_text, label_list, fmts, *args, stdev=True, log_scale=True, legend=True):
    
    fig, ax = plt.subplots()
    
    for arg, label, fmt in zip(args, label_list, fmts):
        
        n_cores = arg[0]
        mean = arg[1]
        if stdev:
            stdevs = arg[2]
        
        if stdev:
            plt.errorbar(n_cores, mean, yerr=stdevs, fmt=fmt, label=label, capsize=3, markersize=4)
        else:
            plt.plot(n_cores, mean, fmt, label=label, markersize=4)


    # axes
    plt.margins(x=0.02, tight=True)
    ax.tick_params(axis='both', direction='in', which='both')
    
    # x axis
    plt.xlabel('number processes')
    plt.xticks(n_processes, n_processes)
    
    
    # y axis
    plt.ylabel('wall clock time [s]')
    if log_scale:
        plt.yscale('log')
    plt.margins(x=0.02, tight=True)
    
    if log_scale:
        plt.text(core_number_text, 60 + 5, '1 minute')
        plt.text(core_number_text, 600 + 40, '10 minutes')
        #plt.text(core_number_text, 1800 + 80, '30 minutes')
        plt.text(core_number_text, 3600 + 200, '1 hour')
        #plt.text(core_number_text, 86400 + 2000, '1 day')
    else:
        plt.yscale('linear')
        y_max = 1100
        plt.ylim(0,y_max)
        interval = y_max/100
        plt.text(core_number_text, 60 + interval, '1 minute')
        plt.text(core_number_text, 600 + interval, '10 minutes')
        #plt.text(core_number_text, 3600 + interval, '1 hour')
        #plt.text(core_number_text, 43200 + interval, '1/2 day')
    
    # grid
    if log_scale:
        ax.axhline(60, linestyle='--', color='grey', linewidth=0.5)
        ax.axhline(600, linestyle='--', color='grey', linewidth=0.5)
        #ax.axhline(1800, linestyle='--', color='grey', linewidth=0.5)
        ax.axhline(3600, linestyle='--', color='grey', linewidth=0.5)
        #ax.axhline(86400, linestyle='--', color='grey', linewidth=0.5)
    else:
        ax.axhline(60, linestyle='--', color='grey', linewidth=0.5)
        ax.axhline(600, linestyle='--', color='grey', linewidth=0.5)
        #ax.axhline(3600, linestyle='--', color='grey', linewidth=0.5)
        #ax.axhline(43200, linestyle='--', color='grey', linewidth=0.5)
        #plt.grid(visible=True, which='major', axis='y', color='grey', linestyle='--', linewidth=0.5)
    
    #plt.title('Time ' + model_name)
    
    # legend
    if legend:
        handles, labels = ax.get_legend_handles_labels()
        handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles] # take out errorbars of legend objects if object stems from plt.errorbar
        plt.legend(handles, labels, loc="center", bbox_to_anchor=(0.5, -0.2))


def get_plotting_input(dictionary, lookup_time, json_file):
    '''Creates list with cores, mean and stdev as input for plotting.'''
    core_list = [int(cores) for cores in dictionary.keys()]
    mean_list = [dictionary[str(cores)][lookup_time]['mean'] for cores in core_list]
    try:
        stdev_list = [dictionary[str(cores)][lookup_time]['stdev'] for cores in core_list]
    except KeyError:
        stdev_list = len(core_list)*[0]
        print(f'No stdevs for {json_file}.')
    return (core_list, mean_list, stdev_list)


if __name__ == '__main__':
    
    start = time()
    
    ### set variables ###
    n_processes = [3, 5, 10, 20, 30, 40, 50, 60]
    core_number_text = 51.5
    model_name = 'mmsyn_sm10'
    outpath = './figures/mfel/'
    log_scale = True
    legend = True
    stdev = False
    ECMproject = True
    mfel = True

    
    ### model names
    ECMproject_json = './times/ECMproject/process/' + model_name + '_mplrs_project_result_times.json'
    mfel_json = './times/ECMproject/mfel/' + model_name + '_mfel_result_times.json'
    

    ### read in json time results ###
    if ECMproject:
        ECMproject_dictionary = json.load(open(ECMproject_json))
    
    if mfel:
        mfel_dictionary = json.load(open(mfel_json))


    ### plot times ###

    if ECMproject and mfel:
        time_list = ['time_total', 'time_initial_setup', 'time_projection', 'time_enumeration', 'time_postprocessing']
        for lookup_time in time_list:
            label_list = [lookup_time + ' ECMproject mfel=False', lookup_time + ' ECMproject mfel=True']
            fmts = ['--bo', '-g^']
            results_ECMproject = get_plotting_input(ECMproject_dictionary, lookup_time, ECMproject_json)
            results_mfel = get_plotting_input(mfel_dictionary, lookup_time, mfel_json)
            plot_time(model_name, n_processes, core_number_text, label_list, fmts, results_ECMproject, results_mfel, stdev=stdev, log_scale=log_scale, legend=legend)
            plt.savefig(outpath + model_name + '_' + lookup_time + '_mfel_false_vs_mfel_true.png', bbox_inches='tight') # called function does not close plt object --> we can manipulate it outside the function until we close it
            plt.close() # its important to close it before we start the next plot_time call

    end = time()
    print(f'Ran in {round(end - start, 5)} seconds.')
