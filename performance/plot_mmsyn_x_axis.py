#!/usr/bin/env python

"""Python script to plot mean times on 20 and 60 cores of ECMproject and ecmtool with all mmsyn models as x-axis.

@author: Christian Mayer
"""

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


if __name__ == '__main__':
    
    ### set up architecture for further performance measurement
    dir_path = os.path.dirname(os.path.realpath(__file__))
    outpath = dir_path + '/figures/complexity/'
    # create figures directory
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    
    ECMproject_inpath = dir_path + '/times/ECMproject/process/'
    ECMproject_pool_inpath = dir_path + '/times/ECMproject/pool/'
    ECMproject_single_inpath = dir_path + '/times/ECMproject/single/'
    ecmtool_inpath = dir_path + '/times/ecmtool/'
    
    ### options
    log_scale = True
    day = True # 1 day mark in y axis
    y_borders = False # if not False, give the border values of the y axis (for log scale use power of.... e.g. (10**-1, 10**5))
    memory = False
    stdev = True
    ecmtool = True
    pool = False
    single = False
    max_complexity = 22 # if 22 --> up to mmsyn_sm22
    legend = True
    lookup_time = 'time_total'
    if ecmtool:
        lookup_time_ecmtool_1 = 'time_total'
        lookup_time_ecmtool_2 = False
            
    ECMproject_label = 'ECMproject'
    ecmtool_label = 'ecmtool'
    
    # range of x-axis for model complexity
    x = list(range(0,max_complexity+1))
    
    if memory:
        ax_label = 'memory [kB]'
        lookup_time = 'max_rss'
        lookup_time_ecmtool_1 = 'max_rss'
        lookup_time_ecmtool_2 = False
    else:
        ax_label = 'wall clock time [s]'
    
    
    ### ECMproject
    # set up model list
    max_model = 22
    x_ECMproject = list(range(0,max_model+1))
    models = ['mmsyn_sm' + str(num).zfill(2) for num in x_ECMproject]
    # process list (which process number shall be investigated on which complexity)
    process_list = ['20']*18 + ['60']*5
    if len(models) == len(process_list):
        ECMproject_means = []
        ECMproject_stdevs = []
        for processes, model in zip(process_list, models):
            # reading in data
            json_file = ECMproject_inpath + model + '_mplrs_project_result_times.json'
            dictionary = load_time_results(json_file)
            # get all means
            means = dictionary[processes][lookup_time]['mean']
            ECMproject_means.append(means)
            if stdev:
                # get all sandard deviations
                try:
                    stdevs = dictionary[processes][lookup_time]['stdev']
                except KeyError:
                    stdevs = 0
                ECMproject_stdevs.append(stdevs)
    else:
        print(f'process_list: {process_list}')
        print(f'models: {models}')
        raise Exception("Process list does not fit to list of models on single.")
    
    
    ### ecmtool
    if ecmtool:
        # set up model list
        max_model = 18
        x_ecmtool = list(range(0,max_model+1))
        models = ['mmsyn_sm' + str(num).zfill(2) for num in x_ecmtool]
        # process list
        process_list = ['20']*18 + ['60']*1
        if len(models) == len(process_list):
            ecmtool_means = []
            ecmtool_stdevs = []
            for processes, model in zip(process_list, models):
                # reading in data
                json_file = ecmtool_inpath + model + '_ecmtool_result_times.json'
                dictionary = load_time_results(json_file)
                # get all means
                means = dictionary[processes][lookup_time_ecmtool_1]['mean']
                if lookup_time_ecmtool_2:
                    means_2 = dictionary[processes][lookup_time_ecmtool_2]['mean']
                    means = means + means_2
                ecmtool_means.append(means)
                if stdev:
                    # get all sandard deviations
                    try:
                        stdevs = dictionary[processes][lookup_time_ecmtool_1]['stdev']
                    except KeyError:
                        stdevs = 0
                    if lookup_time_ecmtool_2:
                        stdevs = 0
                    ecmtool_stdevs.append(stdevs)
        else:
            print(f'process_list: {process_list}')
            print(f'models: {models}')
            raise Exception("Process list does not fit to list of models on single.")
    

    ### ECMproject with pool option
    if pool:
        # set up model list
        max_model = 18
        x_pool = list(range(0,max_model+1))
        models = ['mmsyn_sm' + str(num).zfill(2) for num in x_pool]
        # process list
        process_list = ['20']*18 + ['60']*1
        if len(models) == len(process_list):
            pool_means = []
            pool_stdevs = []
            for processes, model in zip(process_list, models):
                # reading in data
                json_file = ECMproject_pool_inpath + model + '_mplrs_project_result_times.json'
                dictionary = load_time_results(json_file)
                # get all means
                means = dictionary[processes][lookup_time]['mean']
                pool_means.append(means)
                if stdev:
                    # get all sandard deviations
                    try:
                        stdevs = dictionary[processes][lookup_time]['stdev']
                    except KeyError:
                        stdevs = 0
                    pool_stdevs.append(stdevs)
        else:
            print(f'process_list: {process_list}')
            print(f'models: {models}')
            raise Exception("Process list does not fit to list of models on single.")
        
        
    ### ECMproject with single
    if single:
        # set up model list
        max_model = 17
        x_single = list(range(0,max_model+1))
        models = ['mmsyn_sm' + str(num).zfill(2) for num in x_single]
        # process list
        process_list = ['20']*18
        if len(models) == len(process_list):
            single_means = []
            single_stdevs = []
            for processes, model in zip(process_list, models):
                # reading in data
                json_file = ECMproject_single_inpath + model + '_mplrs_project_result_times.json'
                dictionary = load_time_results(json_file)
                # get all means
                means = dictionary[processes][lookup_time]['mean']
                single_means.append(means)
                if stdev:
                    # get all sandard deviations
                    try:
                        stdevs = dictionary[processes][lookup_time]['stdev']
                    except KeyError:
                        stdevs = 0
                    single_stdevs.append(stdevs)
        else:
            print(f'process_list: {process_list}')
            print(f'models: {models}')
            raise Exception("Process list does not fit to list of models on single.")
    
    
    ### plotting
    fig, ax = plt.subplots()
    
    if stdev:
        plt.errorbar(x_ECMproject, ECMproject_means, ECMproject_stdevs, label=ECMproject_label, linestyle='-', color='b', markersize=4, marker='o', capsize=3)
    else:
        plt.plot(x_ECMproject, ECMproject_means, label=ECMproject_label, linestyle='-', color='b', markersize=4, marker='o')
    
    if ecmtool:
        if stdev:
            plt.errorbar(x_ecmtool, ecmtool_means, ecmtool_stdevs, label=ecmtool_label, linestyle='--', color='orange', markersize=4, marker='s', capsize=3)
        else:
            plt.plot(x_ecmtool, ecmtool_means, label=ecmtool_label, linestyle='--', color='orange', markersize=4, marker='s')
    
    if single:
        if stdev:
            plt.errorbar(x_single, single_means, single_stdevs, label='ECMproject not parallel', linestyle='--', color='g', markersize=4, marker='x', capsize=3)
        else:
            plt.plot(x_single, single_means, label='ECMproject not parallel', linestyle='--', color='g', markersize=4, marker='x')
    
    if pool:
        if stdev:
            plt.errorbar(x_pool, pool_means, pool_stdevs, label='ECMproject parallel with -pool', linestyle='-.', color='r', markersize=4, marker='^', capsize=3)
        else:
            plt.plot(x_pool, pool_means, label='ECMproject parallel with pool', linestyle='-.', color='r', markersize=4, marker='^')

    
    ax.set_xlim(-0.5, max_complexity+0.5)
    #plt.margins(x=0.02, tight=True)
    ax.tick_params(axis='both', direction='in', which='both')
    
    # legend
    if legend:
        handles, labels = ax.get_legend_handles_labels()
        handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles] # take out errorbars of legend objects if object stems from plt.errorbar
        # set number of columns for legend
        if ecmtool:
            ncols = 2
        elif pool and single:
            ncols = 2
        else:
            ncols = 1

        if (len(handles) == 3) and single and pool:
            # set order of handle items
            order = [1,2,0] # 0 = parallel not pool; 1 = single, 2 = parallel pool
            plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc="center", bbox_to_anchor=(0.45, -0.2), ncols=ncols)
        else:
            plt.legend(handles, labels, loc="center", bbox_to_anchor=(0.5, -0.16), ncols=ncols)
    
    # x-axis
    plt.xlabel('complexity of the medium')
    plt.xticks(x, x)
    # grey zone
    plt.axvspan(18-0.5, max_complexity+0.5, facecolor='grey', alpha=0.3)
    
    # y-axis
    plt.ylabel(ax_label)
    if log_scale:
        plt.yscale("log")
        if not memory:
            if y_borders:
                y_min, y_max = y_borders
                plt.ylim(y_min, y_max)
            plt.text(0, 60 + 5, '1 minute')
            plt.text(0, 3600 + 200, '1 hour')
            if day:
                plt.text(0, 86400 + 2000, '1 day')
        else:
            if y_borders:
                y_min, y_max = y_borders
                plt.ylim(y_min, y_max)
            plt.text(0, 1000 + 100, '1 MB')
            plt.text(0, 1000000 + 100000, '1 GB')
            plt.text(0, 1000000000 + 100000000, '1 TB')
    else:
        if not memory:
            if y_borders:
                y_min, y_max = y_borders
                plt.ylim(y_min, y_max)
            else:
                plt.ylim(0, 100000)
            interval = y_max/100
            plt.text(0, 60 + interval, '1 minute')
            plt.text(0, 3600 + interval, '1 hour')
            if day:
                plt.text(0, 86400 + interval, '1 day')
        else:
            if y_borders:
                y_min, y_max = y_borders
                plt.ylim(y_min, y_max)
            else:
                plt.ylim(0, 1000000)
            interval = y_max/100
            plt.text(0, 1000 + interval, '1 MB')
            plt.text(0, 1000000 + interval, '1 GB')
            #plt.text(0, 1000000000 + interval, '1 TB')
            
    # grid
    if not memory:
        ax.axhline(60, linestyle='--', color='grey', linewidth=0.5)
        ax.axhline(3600, linestyle='--', color='grey', linewidth=0.5)
        if day:
            ax.axhline(86400, linestyle='--', color='grey', linewidth=0.5)
    else:
        ax.axhline(1000, linestyle='--', color='grey', linewidth=0.5)
        ax.axhline(1000000, linestyle='--', color='grey', linewidth=0.5)
        ax.axhline(1000000000, linestyle='--', color='grey', linewidth=0.5)
        plt.grid(visible=True, which='major', axis='y', color='grey', linestyle='--', linewidth=0.5)
    
    # save and close file
    pre = ''
    if single and pool:
        pre = 'ECMproject_single_pool'
    elif ecmtool:
        pre = 'ECMproject_ecmtool'
    elif pool and not single:
        pre = 'ECMproject_pool'
    elif not pool and single:
        pre = 'ECMproject_single'
    else:
        pre = 'ECMproject'
    
    if log_scale:
        plt.savefig(outpath + pre + '_log_' + lookup_time + '.png', bbox_inches='tight')
    else:
        plt.savefig(outpath + pre + '_' + lookup_time + '.png', bbox_inches='tight')
    plt.close()
