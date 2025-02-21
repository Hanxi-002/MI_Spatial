import pandas as pd
import numpy as np
import scanpy as sc
import anndata 
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
import commot as ct
import pickle
import itertools
#%% load into the data
raw_counts = pd.read_csv("CCC/Commot/raw_x.txt", sep="\t", index_col=0)
raw_counts.columns = raw_counts.columns.str.slice(0, -4)

norm_counts = pd.read_csv("CCC/Commot/norm_x.txt", sep="\t", index_col=0)
norm_counts.columns = norm_counts.columns.str.slice(0, -4)

meta = pd.read_excel("RawData/Final_Annotation_Dutta_sample_info.xlsx", index_col=0)
meta.index = meta['DCCnames']

# proximity of CD68 to AF. 0 encodes furthur away and 1 encodes closer to AF.
proximity = pd.read_csv("ER_SLIDE/Within_Region/121223_hf_mac/Data/y.csv", index_col=0)
proximity.index = proximity.index.str.slice(0, -4)
proximity.columns = ['mac_proximity']

meta_merged = pd.merge(meta, proximity, how='left', left_index=True, right_index=True)
#  0 is further away, 1 is closer to AF and 2 is not available as in not mac cells
meta_merged['mac_proximity'].fillna(2, inplace=True)
meta = meta_merged

meta_hf = meta[meta['Status'] == 'HF']
meta_control = meta[meta['Status'] == 'Control']

raw_counts_hf = raw_counts.loc[:, meta[meta['Status'] == 'HF'].index]
raw_counts_control = raw_counts.loc[:, meta[meta['Status'] == 'Control'].index]
norm_counts_hf = norm_counts.loc[:, meta[meta['Status'] == 'HF'].index]
norm_counts_control = norm_counts.loc[:, meta[meta['Status'] == 'Control'].index]
assert sum(raw_counts_hf.columns != norm_counts_hf.columns) == 0
#%% Initiate anndata object
def initiate_adata(raw_counts, norm_counts):
    """initiate anndata object with raw and normalized counts

    Args:
        raw_counts (df): raw counts, in the shape of genes x cells
        norm_counts (df): normalized counts, in the shape of genes x cells

    Returns:
        anndata: anndata object
    """
    # initiate anndata object with raw counts
    adata = anndata.AnnData(X=raw_counts.T, 
                            obs=pd.DataFrame(index=raw_counts.columns), 
                            var=pd.DataFrame(index=raw_counts.index))

    # add raw data to adata
    adata.raw = adata
    adata.layers["raw_counts"] = adata.raw.X.copy()
    # now add the normalized counts
    adata.X = norm_counts.T
    return adata

def add_metadata(adata, meta):
    """add metadata to anndata object

    Args:
        adata (anndata): anndata object
        meta (df): metadata, in the shape of cells x metadata

    Returns:
        anndata: anndata object with metadata
    """
    adata.obs = meta
    adata.obsm['spatial'] = meta[['ROICoordinateX', 'ROICoordinateY']]
    adata.obsm['spatial'] = np.asarray(adata.obsm['spatial'])
    return adata

adata_hf = initiate_adata(raw_counts_hf, norm_counts_hf)
adata_control = initiate_adata(raw_counts_control, norm_counts_control)
adata_hf = add_metadata(adata_hf, meta_hf)
adata_control = add_metadata(adata_control, meta_control)

#%% Preprocess data
# try log
# try mean centering

def preprocess(adata):
    """preprocess the data

    Args:
        adata (anndata): anndata object

    Returns:
        anndata: preprocessed anndata object
    """
    # log transform
    sc.pp.log1p(adata)
    # mean centering
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.4)

    sc.pl.umap(adata, color=['leiden', 'Status'])

preprocess(adata_hf)
preprocess(adata_control)


#%% CCC

def filter_LR_pairs(adata):
    df_cellchat = ct.pp.ligand_receptor_database(species='human', signaling_type='Secreted Signaling', database='CellChat')
    print(df_cellchat.shape)

    # comeback to do my own filter
    df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata, min_cell_pct=0.001)
    print(df_cellchat_filtered.shape)
    df_cellchat_filtered.columns = ['Ligand', 'Receptor', 'PathwayNames', 'Method']
    return df_cellchat_filtered

df_cellchat_filtered_hf = filter_LR_pairs(adata_hf)
df_cellchat_filtered_control = filter_LR_pairs(adata_control)
# check if the filtered pairs are the same
df_cellchat_filtered_hf == df_cellchat_filtered_control
#df_cellchat_filtered.to_csv('CCC/Commot/df_cellchat_filtered.csv')

#distance_matrix = distance_matrix(adata_hf.obsm["spatial"], adata_hf.obsm["spatial"])

ct.tl.spatial_communication(adata_hf, database_name='cellchat', 
                            df_ligrec=df_cellchat_filtered_hf, dis_thr=30000, 
                            heteromeric=True, pathway_sum=True)

ct.tl.spatial_communication(adata_control, database_name='cellchat', 
                            df_ligrec=df_cellchat_filtered_control, dis_thr=30000, 
                            heteromeric=True, pathway_sum=True)

#pickle.dump(adata_commot_hf, open("CCC/Commot/adata_commot_hf.pkl", "wb"))
#pickle.dump(adata_commot_control, open("CCC/Commot/adata_commot_control.pkl", "wb"))
#ct.tl.communication_direction(adata, database_name='cellchat', pathway_name='IL13-IL13RA2', k=5)

# ct.pl.plot_cell_communication(adata, database_name='cellchat', pathway_name='IL13-IL13RA2', plot_method='grid', background_legend=True,
#     scale=0.00003, ndsize=8, grid_density=0.4, summary='sender', background='cluster', clustering='leiden', cmap='Alphabet',
#     normalize_v = True, normalize_v_quantile=0.995)

#%% plot sum of communication
# first assign the big cell types we want to see in the network
# assigning mac_hf_close, mac_hf_far, mac_control, AF_hf, AF_control, RF_hf, RF_control
#adata = pickle.load(open("CCC/Commot/adata_commot.pkl", "rb"))
cell_types_hf = pd.DataFrame(adata_hf.obs['mac_proximity'])
cell_types_hf[cell_types_hf == 0] = 'mac_hf_far'
cell_types_hf[cell_types_hf == 1] = 'mac_hf_close'
mask = (adata_hf.obs['Status'] == 'HF') & (adata_hf.obs['SegmentLabel'] == 'active fibroblast')
cell_types_hf[mask] = 'AF_hf'
mask = (adata_hf.obs['Status'] == 'HF') & (adata_hf.obs['SegmentLabel'] == 'resting fibroblast')
cell_types_hf[mask] = 'RF_hf'
cell_types_hf.columns = ['cell_types']
# these are the regions that were not included in the ER because they didn't have full 3 AOIs
cell_types_hf.loc['DSP-1001660011816-D-A11'] = 'mac_hf_far'
cell_types_hf.loc['DSP-1001660012268-C-G05'] = 'mac_hf_close'
cell_types_hf.loc['DSP-1001660012268-C-D03'] = 'mac_hf_far'
cell_types_hf.loc['DSP-1001660012268-C-D05'] = 'mac_hf_far'
cell_types_hf.loc['DSP-1001660011816-D-E04'] = 'mac_hf_far'
cell_types_hf.loc['DSP-1001660011816-D-E06'] = 'mac_hf_far'
cell_types_hf.loc['DSP-1001660011816-D-E08'] = 'mac_hf_far'
cell_types_hf['cell_types'].unique()

cell_types_control = pd.DataFrame(adata_control.obs['mac_proximity'])
mask = (adata_control.obs['Status'] == 'Control') & (adata_control.obs['SegmentLabel'] == 'CD 68')
cell_types_control[mask] = 'mac_control'
mask = (adata_control.obs['Status'] == 'Control') & (adata_control.obs['SegmentLabel'] == 'active fibroblast')
cell_types_control[mask] = 'AF_control'
mask = (adata_control.obs['Status'] == 'Control') & (adata_control.obs['SegmentLabel'] == 'resting fibroblast')
cell_types_control[mask] = 'RF_control'
cell_types_control.columns = ['cell_types']
cell_types_control['cell_types'].unique()


def sum_cell_type_comm(adata, cell_types, df_cellchat_filtered):
    """After running commot, sum the communication between each pair of cell types in HF and Control.
       For example, the sending from AF to RF or the sending from RF to AF in HF or Control.
    Args:
        adata (anndata): adata with commot results
        cell_types (df): sample names as index and 1 column with cell type names
        df_cellchat_filtered (df): filtered ligand receptor pairs from commot

    Returns:
        dict: dictionary with cell type pairs as keys and sum of communication as values
    """
    cell_type_pairs = list(itertools.product(cell_types['cell_types'].unique(), cell_types['cell_types'].unique()))
    dict_comm = {pair: 0 for pair in cell_type_pairs}
    # add each pair in cell_type_pairs is a key in dict_comm
    for cell_type_pair in cell_type_pairs:
        print(cell_type_pair)
        for _, row in df_cellchat_filtered.iterrows():
            pair = row['Ligand'] + '-' + row['Receptor']
            pair_comm = 'commot-cellchat-' + pair
            pair_comm = pd.DataFrame(adata.obsp[pair_comm].toarray())
            pair_comm.columns, pair_comm.index = adata.obs.index, adata.obs.index
            r_idx = cell_types[cell_types['cell_types'] == cell_type_pair[0]].index
            c_idx = cell_types[cell_types['cell_types'] == cell_type_pair[1]].index
            pair_comm = pair_comm[pair_comm.index.isin(r_idx)]
            pair_comm = pair_comm.loc[:, pair_comm.columns.isin(c_idx)]
            assert pair_comm.shape[0] == len(r_idx), print(f"{pair_comm.shape[0]} doesn't match with {len(r_idx)}") 
            assert pair_comm.shape[1] == len(c_idx), print(f"{pair_comm.shape[1]} doesn't match with {len(c_idx)}") 
            # get the sum of communication
            sum_comm = pair_comm.sum().sum()
            dict_comm[cell_type_pair] += sum_comm
    return dict_comm


dict_comm_hf = sum_cell_type_comm(adata_hf, cell_types_hf, df_cellchat_filtered_hf)
dict_comm_control = sum_cell_type_comm(adata_control, cell_types_control, df_cellchat_filtered_control)

dict_comm_hf_df = pd.DataFrame(
    [(source, target, weight) for (source, target), weight in dict_comm_hf.items()],
    columns=['source', 'target', 'weight']
)

dict_comm_control_df = pd.DataFrame(
    [(source, target, weight) for (source, target), weight in dict_comm_control.items()],
    columns=['source', 'target', 'weight']
)

dict_comm_hf_df.to_csv('CCC/Commot/CCC_hf_graph.csv')
dict_comm_control_df.to_csv('CCC/Commot/CCC_control_graph.csv')
#%%
# import networkx as nx
# import matplotlib.pyplot as plt
# import pandas as pd
# from matplotlib.patches import FancyArrowPatch
# from matplotlib import patheffects

# # Function to draw curved edges (even without arrows)
# def draw_curved_edges(G, pos, ax, rad=0.2):
#     for edge in G.edges():
#         source, target = edge
#         if source == target:  # self-loop case
#             loop = FancyArrowPatch(posA=pos[source], posB=pos[source], connectionstyle=f"arc3,rad={rad}",
#                                    color='gray', linewidth=2, arrowstyle='-', mutation_scale=15)
#             ax.add_patch(loop)
#         else:
#             curve = FancyArrowPatch(posA=pos[source], posB=pos[target], connectionstyle=f"arc3,rad={rad}",
#                                     color='gray', linewidth=2, arrowstyle='-', mutation_scale=15)
#             ax.add_patch(curve)

# # Sample dictionary with (source, target) as keys and weight as values
# pair_dict = {
#     ('apple', 'banana'): 3, 
#     ('banana', 'apple'): 5, 
#     ('apple', 'orange'): 1, 
#     ('banana', 'orange'): 6,
#     ('apple', 'apple'): 2  # Self-loop
# }

# # Convert the dictionary to a pandas DataFrame
# df = pd.DataFrame(
#     [(source, target, weight) for (source, target), weight in pair_dict.items()],
#     columns=['source', 'target', 'weight']
# )

# # Create a directed graph using NetworkX
# G = nx.DiGraph()

# # Add edges with weights from the DataFrame
# for _, row in df.iterrows():
#     G.add_edge(row['source'], row['target'], weight=row['weight'])

# # Set up the positions of nodes using a circular layout
# pos = nx.circular_layout(G)

# # Create a matplotlib figure and axis
# fig, ax = plt.subplots()

# # Draw the nodes with smaller size
# nx.draw_networkx_nodes(G, pos, node_size=1000, node_color='lightblue', ax=ax)

# # Draw the curved edges manually using FancyArrowPatch
# draw_curved_edges(G, pos, ax, rad=0.2)

# # Move the node labels slightly next to the nodes by adjusting their positions
# label_pos = {key: (x, y + 0.05) for key, (x, y) in pos.items()}
# nx.draw_networkx_labels(G, label_pos, font_size=12, font_color="black", ax=ax)

# # Turn off the axis
# ax.set_axis_off()
# plt.show()



 #%% plot vector field

# plt.figure(figsize=(8, 8))

# df1 = adata.obs[['ROICoordinateX', 'ROICoordinateY', 'SegmentLabel']]
# df1.columns = ['x', 'y', 'label']
# df2 = pd.DataFrame(IL13)
# df2.columns = ['dx', 'dy']
# df2.index = adata.obs.index
# df2  = df2 * 200

# df = pd.concat([df1, df2], axis=1)

# labels = df['label'].unique()
# for l in labels:
#     df_plot = df[df['label'] == l]
#     plt.scatter(df_plot['x'], df_plot['y'], color = "red", s = 10)
#     # Plot vectors using quiver
#     plt.quiver(df_plot['x'], df_plot['y'], df_plot['dx'], df_plot['dy'], angles='xy', scale_units='xy', scale=1, color='blue')

#     # Set limits and labels
#     plt.xlabel('X Coordinate')
#     plt.ylabel('Y Coordinate')
#     plt.title(l+'-IL13-IL13RA1 Sender Signal')
#     plt.grid(False)

#     # Show the plot
#     plt.show()
