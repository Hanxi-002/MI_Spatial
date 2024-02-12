#%%
import sys
sys.path.append("/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/COAnalyses")
from adata_oracle import *
from oracle_links import *
from helper_funcs import *
import numpy as np
import pandas as pd
import dill
import os
import matplotlib.pyplot as plt
import celloracle as co
import scanpy as sc
%matplotlib inline
# %% creating the andata object
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


# #load the pickle file through dill
# with open('adata_oracle_CD68_CO_020124.pkl', 'rb') as f:
#     adata_oracle = dill.load(f)

# adata_oracle.adata

#%% creating the oracle object and get context specific GRN
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


##########################save both the adata and adata_oracle object###########
#adata = co.load_hdf5('allcell_allcondition_adata')
oracle.to_hdf5("CD68.celloracle.oracle")

#adata_oracle = dill.load(open('adata_oracle_CD68_CO_020124.pkl', 'rb'))
adata_oracle.oracle = oracle
dill.dump(adata_oracle, open('adata_oracle_CD68_CO_020124.pkl', 'wb'))

#%%
############################# calculate GRN ##############################
links = oracle.get_links(cluster_name_for_GRN_unit="louvain", alpha=10,
                         verbose_level=10)

# just like oracle, we save link object but also save it as part of the oracle_links
links.to_hdf5(file_path="CD68.celloracle.links")

#perform link analysis
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

######################find the overlap between SLIDE and oracle######################
# load the oracle_link object from pkl file
#oracle_links = dill.load(open('CD68_oracle_links.pkl', 'rb'))
#oracle = co.load_hdf5("CD68.celloracle.oracle")

latent_factors = load_SLIDE_res('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/CD68/121423/results/SLIDE_Results')
latent_factors.keys()

# the TFs that appear in the latent factors
oracle_links.find_TF_overlap_SLIDE(latent_factors = latent_factors, oracle = adata_oracle.oracle)
for c in range(len(oracle_links.overlap_TF_SLIDE.keys())):
    print(oracle_links.overlap_TF_SLIDE[str(c)])

# the TFs that regulates the latent factors
oracle_links.find_gene_TF_SLIDE(latent_factors, adata_oracle.oracle)
oracle_links.gene_TF_SLIDE
oracle_links.linked_TF_SLIDE

######################### perform ridge regression to calculate GRN ####################

oracle_links.fit_ridge_regression(adata_oracle.oracle)
#dill.dump(oracle_links, open('CD68_oracle_links.pkl', 'wb'))
#dill.dump(adata_oracle, open('adata_oracle_CD68_CO_020124.pkl', 'wb'))


############################get the TFs as a list for perturbation##########################
# function in helper_funcs.py
overlap_TF_list = gene_dict_to_list(oracle_links.overlap_TF_SLIDE)
print(overlap_TF_list)
linked_TF_list = threshold_linked_TFs(oracle_links, thresh = 10)
print(linked_TF_list)
len(linked_TF_list)
#linked_TF_list = TF_list


#%% in silicco perturbation
######################### single TF perturbation ####################
#plt.rcParams["font.family"] = "arial"
plt.rcParams["figure.figsize"] = [6,6]
%config InlineBackend.figure_format = 'retina'
plt.rcParams["savefig.dpi"] = 600

save_folder = "figures"
os.makedirs(save_folder, exist_ok=True)

#goi = 'JUND'

for goi in linked_TF_list:
    print(goi)
    plot_gene_hist(adata_oracle.oracle.adata, goi, save_folder)

    calc_cell_identity_shifts(oracle=adata_oracle.oracle, goi=goi, save_folder = save_folder, \
                              scale = 50, perturb_value = 0.0)

    get_optimal_min_mass(n_grid=40, oracle=adata_oracle.oracle)

    calc_vector_fields(min_mass = 0.14, goi = goi, oracle=adata_oracle.oracle, \
                       save_folder = save_folder, scale_simulation = 10)

    plot_vector_filed_on_cluster(oracle=adata_oracle.oracle, goi = goi, save_folder = save_folder, scale_simulation=10)    

    calc_cluster_vector_diff(oracle=adata_oracle.oracle, adata=adata_oracle.oracle.adata, goi=goi, save_folder=save_folder)

