"""Analyzing the hotnet2 results for the MI_Spatial project for Controls
   env: BasePy
"""
import pandas as pd
import os
import sys

# add path so we can source the helper functions
cur_dir = current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(parent_dir)
os.chdir("/ix/djishnu/Hanxi/MI_Spatial/network_analysis/hotnet2/053124_Control/")
import analyze_hotent2_helper as helper

def main(parent_directory, p_thresh, m_thresh):
    # find the folders that have significant results
    significant_dirs = helper.find_significant_folders(parent_directory, p_thresh, m_thresh)
    # find which modules have significant p-values
    significant_dirs = helper.find_significant_module_sizes(significant_dirs, p_thresh, m_thresh)
    # find the significant modules
    sig_comps = helper.find_significant_modules(significant_dirs)
    # find the edges that are in the significant modules (and logic)
    filtered_edges = helper.find_hotnet2_edges(sig_comps, edges)
    return filtered_edges

#%% ################################################
#set parameters and read in the PPI
h5_path = '../HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4.h5'
edges = helper.read_in_ppi(h5_path)
p_thresh = 0.1
m_thresh = 3
###################################################
#%% ################################################
#run the functions for CD68
parent_directory = 'CD68_hotnet_out/homosapiens_binary_co_complex_feb2023-cd_68_control'
filtered_edges = main(parent_directory, p_thresh, m_thresh)
filtered_edges.to_csv('CD68_hotnet_out/filtered_edges.txt', sep = '\t', index = False)
###################################################

#%% ################################################
#run the functions for AF
parent_directory = 'AF_hotnet_out/homosapiens_binary_co_complex_feb2023-active_fibroblast_control'
filtered_edges = main(parent_directory, p_thresh, m_thresh)
filtered_edges.to_csv('AF_hotnet_out/filtered_edges.txt', sep = '\t', index = False)
# no significant modules
###################################################

#%% ################################################
#run the functions for RF
parent_directory = 'RF_hotnet_out/homosapiens_binary_co_complex_feb2023-resting_fibroblast_control'
filtered_edges = main(parent_directory, p_thresh, m_thresh)
filtered_edges.to_csv('RF_hotnet_out/filtered_edges.txt', sep = '\t', index = False)
###################################################

#%% ################################################
#run the functions for CD68 close to AF
parent_directory = 'CD68_close_hotnet_out/homosapiens_binary_co_complex_feb2023-cd_68_close'
filtered_edges = main(parent_directory, p_thresh, m_thresh)
filtered_edges.to_csv('CD68_hotnet_out/filtered_edges.txt', sep = '\t', index = False)
###################################################
