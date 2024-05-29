"""Analyzing the hotnet2 results for the MI_Spatial project
"""
import pandas as pd
import numpy as np
import os
import h5py
import networkx as nx
import matplotlib.pyplot as plt
os.chdir("/ix/djishnu/Hanxi/MI_Spatial/network_analysis/hotnet2/051424/")

def find_significant_folders(parent_directory, p_thresh, m_thresh = 3):
    """Find all directories in the parent directory that contain a significance.txt file with significant results
    Args:
        parent_directory (_str_): the directory above delta folders
        p_thresh (_float_): p value threshold
        m_thresh (int, optional): threshold for module size. Defaults to 3.

    Returns:
        dict: key as the directory path and value as the dataframe of the significance.txt file.
    """    
    significant_paths = {}
    # Walk through all directories in the parent directory
    for dirpath, dirnames, filenames in os.walk(parent_directory):
        if 'significance.txt' in filenames:
            # Construct the full file path
            file_path = os.path.join(dirpath, 'significance.txt')
            # Read the file
            try:
                data = pd.read_csv(file_path, sep = "\t")
                if any((data['p-value'] < p_thresh) & (data['Size'] > m_thresh)):
                    significant_paths[dirpath] = data
            except Exception as e:
                print(f"Failed to read or process {file_path}: {e}")
    return significant_paths

def find_significant_module_sizes(significant_paths, p_thresh, m_thresh = 3):
    """Etract the row in significant.txt that has the p-value < p_thresh and smalles size (but > m_thresh).
       Each delta should only have one row.
    Args:
        significant_paths (_dict_): key as the directory path and value as the dataframe of the significance.txt file.
        p_thresh (_float_): p value threshold
        m_thresh (int, optional): module size threshold. Defaults to 3.
    Returns:
        dict: key as the directory path and value as one row of the significant.txt that is "most significant".
    """    
    # range through the significant_dirs
    for dirpath, data in significant_paths.items():
        sig_module = data[(data['p-value']< p_thresh) & (data['Size'] > m_thresh)]
        # if there are 2 module sizes that are significant
        if len(sig_module) > 1:
            # only get the 1st row of sig_module
            sig_module = sig_module.iloc[0]
        assert len(sig_module) == 1, "There should be only one significant module"
        significant_paths[dirpath] = sig_module
    return significant_paths
            

def find_significant_modules(significant_dirs):
    """Find the significant modules in the significant directories.

    Args:
        significant_dirs (_dict_): dictionary with key as the directory path and value as one row of the dataframe of the significance.txt file.

    Returns:
        _dataframe_: a pandas dataframe with all the significant modules from all selected delta folders. 
    """    
    sig_comps = []
    for dirpath, data in significant_dirs.items():
        comps = pd.read_csv(os.path.join(dirpath, 'components.txt'), sep = "\t", header = None)
        # get the rows that has more than a certain number of non-nan values
        comps = comps[comps.count(axis = 1) > data['Size'].values[0]]
        sig_comps.append(comps)
    sig_comps = pd.concat(sig_comps)
    return sig_comps

def read_in_ppi(h5_path):
    """Read in the PPI from the h5 file

    Args:
        h5_path (_str_): path to the h5 file   

    Returns:
        _dataframe_: a pandas dataframe with 2 columns: the source and target of the edges
    """
    with h5py.File(h5_path, 'r') as file:
        print("Keys: %s" % file.keys())
        edges = pd.DataFrame(file['edges'][:])
    for column in edges.columns:
        edges[column] = edges[column].str.decode('utf-8')
    edges.columns = ['source', 'target']
    return edges

def find_hotnet2_edges(sig_comps, edges):
    """extract the edges that are in the significant components from the reference PPI

    Args:
        sig_comps (_df_): a dataframe with the significant components across delta
        edges (_df_): edges from the h5 file

    Returns:
        _dataframe_: dataframe contain the final edges
    """    
    gene_names = set()
    for column in sig_comps.columns:
        gene_names.update(sig_comps[column].unique())
    filtered_edges = edges[edges.apply(lambda x: x['source'] in gene_names and x['target'] in gene_names, axis=1)]
    return filtered_edges


#%% run the functions
h5_path = '../HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4.h5'
edges = read_in_ppi(h5_path)

parent_directory = 'CD68_hotnet_out/homosapiens_binary_co_complex_feb2023-cd_68'
significant_dirs = find_significant_folders(parent_directory, p_thresh = 0.1, m_thresh = 3)
significant_dirs = find_significant_module_sizes(significant_dirs, p_thresh = 0.1, m_thresh = 3)
sig_comps = find_significant_modules(significant_dirs)
filtered_edges = find_hotnet2_edges(sig_comps, edges)
filtered_edges.to_csv('filtered_edges.csv', index = False)

G = nx.from_pandas_edgelist(filtered_edges, 'source' , 'target', create_using=nx.Graph())
# Draw the network
pos = nx.spectral_layout(G)
nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=200, edge_color='k', linewidths=1, font_size=5)
plt.show()

