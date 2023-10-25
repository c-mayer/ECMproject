#!/usr/bin/env python

"""Python script to plot times of each section of ECMproject against ecmtool and mfel against mplrs.

@author: Christian Mayer
"""

### plotting scheme ###

## e_coli_core: (repeats = 5)
# plot own, ecmtool and mfel on total time --> show, that mfel is quite slow
# plot only own and mfel comparison of each time --> show, that mfel is slower here, postprocessing

## mmsyn: (repeats = 3)

## e_coli_core combo with mmsyn??:
# only ecmtool and own
# plot total time
# plot preprocessing vs initial setup
# plot first vertex enumeration + intermediate processing vs projection #(first vertex enumeration + intermediate processing) > projection
# plot second vertex enumeration + postprocessing vs enumeration + postprocessing #(second vertex enumeration + postprocessing) == (enumeration + postprocessing)
# plot second vertex enumeration vs enumeration
# plot postprocessing vs postprocessing

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import container
import statistics as stat
import json
from time import time


def plot_time(core_name, n_cores, lowest_core_n, label_list, fmts, *args, stdev=True, log_scale=True):
    
    fig, ax = plt.subplots()
    
    for arg, label, fmt in zip(args, label_list, fmts):

        mean = [x[0] for x in arg]
        if stdev:
            stdevs = [x[1] for x in arg]
        
        if stdev:
            plt.errorbar(n_cores, mean, yerr=stdevs, fmt=fmt, label=label, capsize=3, markersize=4)
        else:
            plt.plot(n_cores, mean, fmt, label=label, markersize=4)
    
    plt.ylabel('wall clock time [s]')
    if log_scale:
        plt.yscale('log')
    plt.margins(x=0.02, tight=True)
    plt.xlabel('number processes')
    
    
    if log_scale:
        plt.text(lowest_core_n, 60 + 5, '1 minute')
        plt.text(lowest_core_n, 600 + 40, '10 minutes')
        plt.text(lowest_core_n, 3600 + 200, '1 hour')
        #plt.text(lowest_core_n, 86400 + 2000, '1 day')
    else:
        plt.yscale('linear')
        y_max = 1100
        plt.ylim(0,y_max)
        interval = y_max/100
        plt.text(lowest_core_n, 60 + interval, '1 minute')
        plt.text(lowest_core_n, 600 + interval, '10 minutes')
        #plt.text(lowest_core_n, 3600 + interval, '1 hour')
        #plt.text(lowest_core_n, 43200 + interval, '1/2 day')
    
    # grid
    if log_scale:
        ax.axhline(60, linestyle='--', color='grey', linewidth=0.5)
        ax.axhline(600, linestyle='--', color='grey', linewidth=0.5)
        ax.axhline(3600, linestyle='--', color='grey', linewidth=0.5)
        #ax.axhline(86400, linestyle='--', color='grey', linewidth=0.5)
    else:
        ax.axhline(60, linestyle='--', color='grey', linewidth=0.5)
        ax.axhline(600, linestyle='--', color='grey', linewidth=0.5)
        #ax.axhline(3600, linestyle='--', color='grey', linewidth=0.5)
        #ax.axhline(43200, linestyle='--', color='grey', linewidth=0.5)
        #plt.grid(visible=True, which='major', axis='y', color='grey', linestyle='--', linewidth=0.5)
    
    #plt.title('Time ' + core_name)
    
    # legend
    #handles, labels = ax.get_legend_handles_labels()
    #handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles] # take out errorbars of legend objects if object stems from plt.errorbar
    #plt.legend(handles, labels, loc="center", bbox_to_anchor=(0.5, -0.2))


def create_plotting_input(lookup_time, n_cores, **kwargs):  
    '''Creates list with mean and stdev as input for plotting. **kwargs can be: own_dictionary, ecmtool_dictionary and mfel_dictionary.'''
    result = []
    if 'own_dictionary' in kwargs.keys():
        mean_stdev_list_own = [(own_dictionary[str(cores)][lookup_time]['mean'], own_dictionary[str(cores)][lookup_time]['stdev']) for cores in n_cores]
        result.append(mean_stdev_list_own)
    if 'ecmtool_dictionary' in kwargs.keys():
        mean_stdev_list_ecmtool = [(ecmtool_dictionary[str(cores)][lookup_time]['mean'], ecmtool_dictionary[str(cores)][lookup_time]['stdev']) for cores in n_cores]
        result.append(mean_stdev_list_ecmtool)
    if 'mfel_dictionary' in kwargs.keys():
        mean_stdev_list_mfel = [(mfel_dictionary[str(cores)][lookup_time]['mean'], mfel_dictionary[str(cores)][lookup_time]['stdev']) for cores in n_cores]
        result.append(mean_stdev_list_mfel)
    if len(result) <= 1:
        result = result[0]
    return result


if __name__ == '__main__':
    
    start = time()
    
    ### set variables ###
    n_cores = ['3', '5', '10', '20', '30', '40', '50', '60']
    lowest_core_n = '3'
    #own_tool_json = './e_coli_core_corrected_mplrs_project_result_times.json'
    own_tool_json = './times/ECMproject/pool/mmsyn_sm10_mplrs_project_result_times.json'
    #ecmtool_json = './e_coli_core_corrected_ecmtool_result_times.json'
    #ecmtool_json = './mmsyn_sm10_mplrs_project_result_times.json'
    mfel_json = './times/ECMproject/mfel/mmsyn_sm10_mfel_result_times.json'
    core_name = 'mmsyn_sm10'
    outpath = './figures/mfel/'
    log_scale = True
    stdev = False
    mfel = True
    ecmtool = False
    
    
    ### read in json time results ###
    own_dictionary = json.load(open(own_tool_json))
    
    ## times ##
    #total
    #initial_setup
    #projection
    #enumeration
    #postprocessing
    
    if ecmtool:
        ecmtool_dictionary = json.load(open(ecmtool_json))
    
    ## times ##
    #total
    #preprocessing
    #first vertex enumeration
    #intermediate processing
    #second vertex enumeration
    #postprocessing
    
    if mfel:
        mfel_dictionary = json.load(open(mfel_json))
    
    ## times ##
    #total
    #initial_setup
    #projection
    #enumeration
    #postprocessing

    
    ### create input datasets for plotting ###
    # total
    if mfel and ecmtool:
        mean_stdev_own_total, mean_stdev_ecmtool_total, mean_stdev_mfel_total = create_plotting_input('time_total', n_cores, own_dictionary=own_dictionary, ecmtool_dictionary=ecmtool_dictionary, mfel_dictionary=mfel_dictionary)
    elif not mfel and ecmtool:
        mean_stdev_own_total, mean_stdev_ecmtool_total = create_plotting_input('time_total', n_cores, own_dictionary=own_dictionary, ecmtool_dictionary=ecmtool_dictionary)
    elif mfel and not ecmtool:
        mean_stdev_own_total, mean_stdev_mfel_total = create_plotting_input('time_total', n_cores, own_dictionary=own_dictionary, mfel_dictionary=mfel_dictionary)
    else:
        mean_stdev_own_total = create_plotting_input('time_total', n_cores, own_dictionary=own_dictionary)
    
    # initial_setup, preprocessing
    mean_stdev_own_initial_setup = create_plotting_input('time_initial_setup', n_cores, own_dictionary=own_dictionary)
    if ecmtool:
        mean_stdev_ecmtool_preprocessing = create_plotting_input('time_preprocessing', n_cores, ecmtool_dictionary=ecmtool_dictionary)
    if mfel:
        mean_stdev_mfel_initial_setup = create_plotting_input('time_initial_setup', n_cores, mfel_dictionary=mfel_dictionary)
    
    # projection, first vertex enumeration
    mean_stdev_own_projection = create_plotting_input('time_projection', n_cores, own_dictionary=own_dictionary)
    if ecmtool:
        mean_stdev_ecmtool_first_vertex_enumeration = create_plotting_input('time_first_vertex_enumeration', n_cores, ecmtool_dictionary=ecmtool_dictionary)
    if mfel:
        mean_stdev_mfel_projection = create_plotting_input('time_projection', n_cores, mfel_dictionary=mfel_dictionary)
    
    # intermediate processing
    if ecmtool:
        mean_stdev_ecmtool_intermediate_processing = create_plotting_input('time_intermediate_processing', n_cores, ecmtool_dictionary=ecmtool_dictionary)
    
    # enumeration, second vertex enumeration
    mean_stdev_own_enumeration = create_plotting_input('time_enumeration', n_cores, own_dictionary=own_dictionary)
    if ecmtool:
        mean_stdev_ecmtool_second_vertex_enumeration = create_plotting_input('time_second_vertex_enumeration', n_cores, ecmtool_dictionary=ecmtool_dictionary)
    if mfel:
        mean_stdev_mfel_enumeration = create_plotting_input('time_enumeration', n_cores, mfel_dictionary=mfel_dictionary)
    
    # postprocessing
    mean_stdev_own_postprocessing = create_plotting_input('time_postprocessing', n_cores, own_dictionary=own_dictionary)
    if ecmtool:
        mean_stdev_ecmtool_postprocessing = create_plotting_input('time_postprocessing', n_cores, ecmtool_dictionary=ecmtool_dictionary)
    if mfel:
        mean_stdev_mfel_postprocessing = create_plotting_input('time_postprocessing', n_cores, mfel_dictionary=mfel_dictionary)
    
    ## calculate wanted times
    
    # first vertex enumeration + intermediate processing vs projection
    if ecmtool:
        all_values_ecmtool_first_vert_inter = [np.array(ecmtool_dictionary[str(cores)]['time_first_vertex_enumeration']['all_values']) + np.array(ecmtool_dictionary[str(cores)]['time_intermediate_processing']['all_values']) for cores in n_cores]
        # first_vert_inter
        mean_stdev_ecmtool_first_vert_inter = [(stat.mean(all_values), stat.stdev(all_values)) for all_values in all_values_ecmtool_first_vert_inter]
    
    # second vertex enumeration + postprocessing vs enumeration + postprocessing
    all_values_own_enum_post = [np.array(own_dictionary[str(cores)]['time_enumeration']['all_values']) + np.array(own_dictionary[str(cores)]['time_postprocessing']['all_values']) for cores in n_cores]
    # enum_post
    mean_stdev_own_enum_post = [(stat.mean(all_values), stat.stdev(all_values)) for all_values in all_values_own_enum_post]
    if ecmtool:
        all_values_ecmtool_sec_vert_post = [np.array(ecmtool_dictionary[str(cores)]['time_second_vertex_enumeration']['all_values']) + np.array(ecmtool_dictionary[str(cores)]['time_postprocessing']['all_values']) for cores in n_cores]
        # sec_vert_post
        mean_stdev_ecmtool_sec_vert_post = [(stat.mean(all_values), stat.stdev(all_values)) for all_values in all_values_ecmtool_sec_vert_post]


    ### plot times ###
    
    ## total time
    if mfel and ecmtool:
        label_list = ['total_time ECMproject mfel=False', 'total_time ecmtool', 'total_time ECMproject mfel=True']
        fmts = ['-bo', '-rs', '-g^']
        plot_time(core_name, n_cores, lowest_core_n, label_list, fmts, mean_stdev_own_total, mean_stdev_ecmtool_total, mean_stdev_mfel_total, stdev=stdev, log_scale=log_scale)
    if mfel and not ecmtool:
        label_list = ['total_time ECMproject mfel=False', 'total_time ECMproject mfel=True']
        fmts = ['-bo', '-g^']
        plot_time(core_name, n_cores, lowest_core_n, label_list, fmts, mean_stdev_own_total, mean_stdev_mfel_total, stdev=stdev, log_scale=log_scale)
    if ecmtool and not mfel:
        label_list = ['total_time ECMproject', 'total_time ecmtool']
        fmts = ['-bo', '-rs']
        plot_time(core_name, n_cores, lowest_core_n, label_list, fmts, mean_stdev_own_total, mean_stdev_ecmtool_total, stdev=stdev, log_scale=log_scale)
    if not ecmtool and not mfel:
        label_list = ['total_time ECMproject']
        fmts = ['-bo']
        plot_time(core_name, n_cores, lowest_core_n, label_list, fmts, mean_stdev_own_total, stdev=stdev, log_scale=log_scale)
    plt.savefig(outpath + core_name + '_time_total_mfel_false_vs_mfel_true.png', bbox_inches='tight') # called function does not close plt object --> we can manipulate it outside the function until we close it
    plt.close() # its important to close it before we start the next plot_time call
    
    ## comparison mfel=False vs mfel=True on each different time (show where mfel is slow)
    if mfel:
        time_list = ['time_initial_setup', 'time_projection', 'time_enumeration', 'time_postprocessing']
        for new_time in time_list:
            label_list = [new_time + ' ECMproject mfel=False', new_time + ' ECMproject mfel=True']
            fmts = ['-bo', '-g^']
            mean_stdev_own, mean_stdev_mfel = create_plotting_input(new_time, n_cores, own_dictionary=own_dictionary, mfel_dictionary=mfel_dictionary)
            plot_time(core_name, n_cores, lowest_core_n, label_list, fmts, mean_stdev_own, mean_stdev_mfel, stdev=stdev, log_scale=log_scale)
            plt.savefig(outpath + core_name + '_' + new_time + '_mfel_false_vs_mfel_true.png', bbox_inches='tight') # called function does not close plt object --> we can manipulate it outside the function until we close it
            plt.close() # its important to close it before we start the next plot_time call

    
    ## compare each different time on ECMproject and ecmtool
    if ecmtool:
        label_list = ['total_time ECMproject', 'total_time ecmtool', 'preprocessing ECMproject', 'preprocessing ecmtool', 'projection ECMproject', 'first_vertex_enumeration ecmtool', 'intermediate_processing ecmtool', 'enumeration ECMproject', 'second_vertex_enumeration ecmtool', 'postprocessing ECMproject', 'postrpocessing ecmtool']
        fmts = ['-ro', '-rs', '-go', '-gs', '-bo', '-bs', '-ks', '-mo', '-ms', '-co', '-cs']
        plot_time(core_name, n_cores, lowest_core_n, label_list, fmts, mean_stdev_own_total, mean_stdev_ecmtool_total, mean_stdev_own_initial_setup, mean_stdev_ecmtool_preprocessing, mean_stdev_own_projection, mean_stdev_ecmtool_first_vertex_enumeration, mean_stdev_ecmtool_intermediate_processing, mean_stdev_own_enumeration, mean_stdev_ecmtool_second_vertex_enumeration, mean_stdev_own_postprocessing, mean_stdev_ecmtool_postprocessing, stdev=stdev, log_scale=log_scale)
        plt.savefig(outpath + core_name + '_all_values.png', bbox_inches='tight') # called function does not close plt object --> we can manipulate it outside the function until we close it
        plt.close() # its important to close it before we start the next plot_time call

        ## first vertex enumeration + intermediate processing vs projection
        plot_time(core_name, n_cores, lowest_core_n, ['projection ECMproject', 'first vertex enumeration + intermediate processing ecmtool'], ['-bo', '-rs'], mean_stdev_own_projection, mean_stdev_ecmtool_first_vert_inter, stdev=stdev, log_scale=log_scale)
        plt.savefig(outpath + core_name + '_sum_proj.png', bbox_inches='tight') # called function does not close plt object --> we can manipulate it outside the function until we close it
        plt.close() # its important to close it before we start the next plot_time call

        ## second vertex enumeration + postprocessing vs enumeration + postprocessing
        plot_time(core_name, n_cores, lowest_core_n, ['enumeration + postprocessing ECMproject', 'second vertex enumeration + postprocessing ecmtool'], ['-bo', '-rs'], mean_stdev_own_enum_post, mean_stdev_ecmtool_sec_vert_post, stdev=stdev, log_scale=log_scale)
        plt.savefig(outpath + core_name + '_sum_enum_post.png', bbox_inches='tight') # called function does not close plt object --> we can manipulate it outside the function until we close it
        plt.close() # its important to close it before we start the next plot_time call
    
    end = time()
    print(f'Ran in {round(end - start, 5)} seconds.')
