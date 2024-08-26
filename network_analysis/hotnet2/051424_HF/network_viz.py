import pandas as pd 
import numpy as np
import os

os.chdir("/ix/djishnu/Hanxi/MI_Spatial/network_analysis/hotnet2/051424_HF/")

# read all txt files in a foler
def read_SLIDE_LFs(folder):
    """read all SLIDE latent factors from results folder

    Args:
        folder (_str_): the string of path containing SLIDE results

    Returns:
        dfs: a pandas dataframe containing all latent factors
    """
    files = os.listdir(folder)
    txt_files = [f for f in files if f.endswith('.txt')]
    dfs = []
    for f in txt_files:
        dfs.append(pd.read_csv(folder + f, sep = '\t'))
    df = pd.concat(dfs)
    return df

def get_overlapping_nodes(df, filtered_edges, linked_TFs, human_TFs):
    """get the overlapping nodes between SLIDE and HotNet2

    Args:
        df (_pandas dataframe_): the dataframe containing SLIDE latent factors
        filtered_edges (_pandas dataframe_): the dataframe containing HotNet2 edges

    Returns:
        node_attributes: a pandas dataframe containing node attributes. \
                         overlap==1 if the node is in SLIDE, 2 if the node is in linked_TFs.
                         TF==1 if the node is a TF.
    """
    unique_nodes = pd.Series(list(set(filtered_edges['source']).union(filtered_edges['target'])))
    node_attributes = pd.DataFrame({'node': unique_nodes})
    node_attributes['overlap'] = node_attributes['node'].isin(df['names']).astype(int)
    node_attributes['overlap'] = node_attributes.apply(
        lambda x: 2 if x['node'] in linked_TFs else x['overlap'], axis=1
    )
    node_attributes['TF'] = node_attributes['node'].isin(human_TFs[0]).astype(int)
    return node_attributes


#%% ############################# 
# AF
human_TFs = pd.read_csv("../../../Cell_Oracle/allTFs_hg38_Scenic.txt", header = None)
filtered_edges = pd.read_csv("AF_hotnet_out/filtered_edges.txt", sep = '\t')
df = read_SLIDE_LFs("../../../ER_SLIDE/ActiveFibro/121223/results/SLIDE_Results/")

# made this list from the figures in /ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/Active_Fibro/021324/linked_figures
linked_TFs = ['ATF3', 'EP300', 'KLF2', 'KLF4', 'KLF6', 'KLF12', 'MYC', 'REST', 
              'SP3', 'TAL1', 'THAP1', 'YY1']

AF_node_attribute = get_overlapping_nodes(df, filtered_edges, linked_TFs, human_TFs)
AF_node_attribute.to_csv("AF_hotnet_out/AF_node_attribute.txt", sep = '\t', index = False)
#%% #############################
# RF
filtered_edges = pd.read_csv("RF_hotnet_out/filtered_edges.txt", sep = '\t')
df = read_SLIDE_LFs("../../../ER_SLIDE/RestingFibro/121423/results/SLIDE_Results/")

linked_TFs = ["BCLAF1", "E2F4", "JUND", "KLF2", "KLF6", "KLF13", "MAZ", "SP3", "STAT2", "ZBTZ7A"]

RF_node_attribute = get_overlapping_nodes(df, filtered_edges, linked_TFs, human_TFs)
RF_node_attribute.to_csv("RF_hotnet_out/RF_node_attribute.txt", sep = '\t', index = False)
#%% #############################
# CD68
filtered_edges = pd.read_csv("CD68_hotnet_out/filtered_edges.txt", sep = '\t')
df = read_SLIDE_LFs("../../../ER_SLIDE/CD68/121423/results/SLIDE_Results/")

linked_TFs = ['ATF3', 'EGR1', "EP300", "IRF1", "KLF3", "MAFK", "MAZ"]

CD68_node_attribute = get_overlapping_nodes(df, filtered_edges, linked_TFs, human_TFs)
CD68_node_attribute.to_csv("CD68_hotnet_out/CD68_node_attribute.txt", sep = '\t', index = False)

#%% #############################
# CD68
filtered_edges = pd.read_csv("CD68_close_hotnet_out/filtered_edges.txt", sep = '\t')
df = read_SLIDE_LFs("../../../ER_SLIDE/Within_Region/121223_hf_mac/results/SLIDE_Results/")

# gene names from /ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/CD68/021424/linked_figures
linked_TFs = ['ATF3', 'EGR1', "EP300", "IRF1", "KLF3", "MAFK", "MAZ"]

CD68_close_node_attribute = get_overlapping_nodes(df, filtered_edges, linked_TFs, human_TFs)
CD68_close_node_attribute.to_csv("CD68_close_hotnet_out/CD68_close_node_attribute.txt", sep = '\t', index = False)
