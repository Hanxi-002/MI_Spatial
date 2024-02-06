#%%
import sys
sys.path.append("/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/COAnalyses")
from adata_oracle import *
from oracle_links import *
import numpy as np
import pandas as pd
import dill
import os
import matplotlib.pyplot as plt
import celloracle as co
import scanpy as sc
%matplotlib inline
# %%
print(os.getcwd())
os.chdir("/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/CD68/020124")


adata_oracle = adata_oracle('CD68_norm.csv')
adata_oracle.norm_data.shape
adata_oracle.adata.X.shape
adata_oracle.qc_filter()
adata_oracle.adata.X.shape
adata_oracle.add_raw_count('CD68.csv')


meta = pd.read_excel('/ix/djishnu/Hanxi/MI_Spatial/RawData/Final_Annotation_Dutta_sample_info.xlsx', index_col=0)
meta = meta[meta['CD 68'] == True]
celltype_df = meta[['DCCnames', 'Status']]
celltype_df.index = meta['DCCnames'] + ".dcc"

assert celltype_df.index.isin(adata_oracle.adata.obs.index).sum() == adata_oracle.norm_data.shape[1], "length of meta data not match with data"
celltype_df = celltype_df.sort_index()
# make sure that celltype_df is in the same order as data
assert sum(celltype_df.index == adata_oracle.adata.obs.index) == adata_oracle.norm_data.shape[1], "index of meta data not match with data"

adata_oracle.adata.obs['cell_type'] = 'CD 68'
adata_oracle.adata.obs['Status'] = celltype_df['Status']

adata_oracle.mean_centering()
adata_oracle.pca()

sc.pp.neighbors(adata_oracle.adata, n_neighbors=3, n_pcs=50)
sc.tl.leiden(adata_oracle.adata)
sc.tl.louvain(adata_oracle.adata)
sc.pl.umap(adata_oracle.adata, color=['leiden', 'louvain', 'Status'])

adata_oracle.plot_dim_reduction(method = "umap", color = "louvain")
adata_oracle.plot_dim_reduction(method = "umap", color = "Status")

dill.dump(adata_oracle, open('adata_oracle_CD68_CO_020124.pkl', 'wb'))




#load the pickle file through dill
with open('adata_oracle_CD68_CO_020124.pkl', 'rb') as f:
    adata_oracle = dill.load(f)

adata_oracle.adata

#%%

oracle = co.Oracle()
print("Metadata columns :", list(adata_oracle.adata.obs.columns))
print("Dimensional reduction: ", list(adata_oracle.adata.obsm.keys()))
adata_oracle.adata.X = adata_oracle.adata.layers["raw_count"].copy()

oracle.import_anndata_as_raw_count(adata=adata_oracle.adata,
                                   cluster_column_name='louvain',
                                   embedding_name='X_umap')

base_GRN = co.data.load_human_promoter_base_GRN()
oracle.import_TF_data(TF_info_matrix=base_GRN)

oracle.perform_PCA()

# Select important PCs
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.show()
print(n_comps)
n_comps = min(n_comps, 50)

n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")
k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)

#adata = co.load_hdf5('allcell_allcondition_adata')
oracle.to_hdf5("CD68.celloracle.oracle")


############## calculate GRN
links = oracle.get_links(cluster_name_for_GRN_unit="louvain", alpha=10,
                         verbose_level=10)

links.to_hdf5(file_path="CD68.celloracle.links")
#%%

oracle_links = oracle_links("CD68.celloracle.links")
oracle_links.links.links_dict['0'].head()
oracle_links.filter_pval()
oracle_links.calc_weighted_logp()
oracle_links.links.links_dict['0'].head()

# load the human TF list
# use pandas to read a txt file where each row contains one gene nam
base_GRN = co.data.load_human_promoter_base_GRN()
human_TF = pd.read_csv('/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/allTFs_hg38_Scenic.txt', header=None)
oracle_links.filter_human_TF(human_TF=human_TF)

#this might take a few minutes
oracle_links.get_top_links(percentile = 0.1)
oracle_links.TF_degree_dict # the degree of each NODE in each cluster
oracle_links.filtered_links_dict # the top links for each NODE in each cluster

#%%
def load_SLIDE_res(folder_path):
    """Load SLIDE restulst as a dictionary.
        Args:
            folder_path (str): the path to the folder containing all SLIDE results.
    """
    latent_factors = {}  
    for filename in os.listdir(folder_path):
        if filename.endswith('.txt'): 
            file_path = os.path.join(folder_path, filename)
            df = pd.read_csv(file_path, delimiter='\t')
            latent_factors[filename[0:-4]] = df
    return latent_factors
#%%

latent_factors = load_SLIDE_res('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/CD68/121423/results/SLIDE_Results')
latent_factors.keys()

oracle_links.find_TF_overlap_SLIDE(latent_factors = latent_factors, human_TF = human_TF)
for c in range(len(oracle_links.degree_1_overlap.keys())):
    print(oracle_links.degree_1_overlap[str(c)])

oracle_links.find_gene_TF_SLIDE(latent_factors, human_TF)
oracle_links.gene_TF_SLIDE
oracle_links.linked_TF_SLIDE

dill.dump(oracle_links, open('CD68_oracle_links.pkl', 'wb'))


# %% in silico perturbation'
print(co.__version__)
#plt.rcParams["font.family"] = "arial"
plt.rcParams["figure.figsize"] = [6,6]
%config InlineBackend.figure_format = 'retina'
plt.rcParams["savefig.dpi"] = 600

save_folder = "figures"
os.makedirs(save_folder, exist_ok=True)

oracle = co.load_hdf5("CD68.celloracle.oracle")

oracle_links.links.filtered_links = oracle_links.filtered_links_dict
oracle.get_cluster_specific_TFdict_from_Links(links_object = oracle_links.links)
oracle.fit_GRN_for_simulation(alpha=10,
                              use_cluster_specific_TFdict=True)



#%%
overlap_TF = []
for df in oracle_links.degree_1_overlap.values():
    overlap_TF.extend(df['overlap'].tolist())
overlap_TF = list(set(overlap_TF))
len(overlap_TF)
#%%
overlap_TF_toy = ['SIRT6',
 'JUNB',
 'JUND',
 'ZNF680',
 'MBD2',
 'MAZ',
 'CEBPD']

for goi in overlap_TF_toy:
    goi  = "JUNB"
    sc.get.obs_df(oracle.adata, keys=goi, layer="imputed_count").hist()
    plt.show()
    oracle.simulate_shift(perturb_condition={goi: 0.0},
                        n_propagation=3)

    # Get transition probability
    oracle.estimate_transition_prob(n_neighbors=50,
                                    knn_random=True,
                                    sampled_fraction=1)

    # Calculate embedding
    oracle.calculate_embedding_shift(sigma_corr=0.05)

    fig, ax = plt.subplots(1, 2,  figsize=[13, 6])

    scale = 50
    # Show quiver plot
    oracle.plot_quiver(scale=scale, ax=ax[0])
    ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")

    # Show quiver plot that was calculated with randomized graph.
    oracle.plot_quiver_random(scale=scale, ax=ax[1])
    ax[1].set_title(f"Randomized simulation vector")

    plt.show()

# %%
n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=50)
oracle.suggest_mass_thresholds(n_suggestion=12)

min_mass = 0.14
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

fig, ax = plt.subplots(1, 2,  figsize=[13, 6])
scale_simulation = 10
# Show quiver plot
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")

# Show quiver plot that was calculated with randomized graph.
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")

plt.show()

fig, ax = plt.subplots(figsize=[8, 8])

oracle.plot_cluster_whole(ax=ax, s=10)
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)
# %%

oracle.delta_embedding.shape
oracle.delta_embedding_random
adata_oracle.adata.obs