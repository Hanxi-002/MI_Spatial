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
#%%

print(os.getcwd())
os.chdir("/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/Active_Fibro/021323/")

adata_oracle = adata_oracle('active_fibro_norm.csv')
adata_oracle.norm_data.shape
adata_oracle.adata.X.shape
adata_oracle.qc_filter()
adata_oracle.adata.X.shape
adata_oracle.add_raw_count('active_fibro.csv')


meta = pd.read_excel('/ix/djishnu/Hanxi/MI_Spatial/RawData/Final_Annotation_Dutta_sample_info.xlsx', index_col=0)
meta = meta[meta['active fibroblast'] == True]
celltype_df = meta[['DCCnames', 'Status']]
celltype_df.index = meta['DCCnames'] + ".dcc"

assert celltype_df.index.isin(adata_oracle.adata.obs.index).sum() == adata_oracle.norm_data.shape[1], "length of meta data not match with data"
celltype_df = celltype_df.sort_index()
# make sure that celltype_df is in the same order as data
assert sum(celltype_df.index == adata_oracle.adata.obs.index) == adata_oracle.norm_data.shape[1], "index of meta data not match with data"

adata_oracle.adata.obs['cell_type'] = 'active fibroblast'
adata_oracle.adata.obs['Status'] = celltype_df['Status']

adata_oracle.mean_centering()
adata_oracle.pca()

sc.pp.neighbors(adata_oracle.adata, n_neighbors=3, n_pcs=50)
adata_oracle.umap()
sc.tl.leiden(adata_oracle.adata)
sc.tl.louvain(adata_oracle.adata)

sc.pl.umap(adata_oracle.adata, color=['leiden', 'louvain', 'Status'], save = 'active_fibroblast_umap.pdf')
#%%
#%% creating the oracle object and get context specific GRN
##############################################
         #creating the oracle object
##############################################
oracle = co.Oracle()
print("Metadata columns :", list(adata_oracle.adata.obs.columns))
print("Dimensional reduction: ", list(adata_oracle.adata.obsm.keys()))
adata_oracle.adata.X = adata_oracle.adata.layers["raw_count"].copy()

oracle.import_anndata_as_raw_count(adata=adata_oracle.adata,
                                   cluster_column_name='leiden',
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
oracle.to_hdf5("active_fibroblast.celloracle.oracle")
adata_oracle.oracle = oracle
dill.dump(adata_oracle, open('adata_oracle_active_fibroblast_CO_021324.pkl', 'wb'))
#%%
######################################
         # calculate GRN
##############################################
links = oracle.get_links(cluster_name_for_GRN_unit="leiden", alpha=10,
                         verbose_level=10)

# just like oracle, we save link object but also save it as part of the oracle_links
links.to_hdf5(file_path="active_fibroblast.celloracle.links")
#%%
links = co.load_hdf5("active_fibroblast.celloracle.links")
oracle_links = oracle_links("active_fibroblast.celloracle.links")
oracle_links.links.links_dict['0'].head()
oracle_links.filter_pval()
oracle_links.calc_weighted_logp()
oracle_links.links.links_dict['0'].head()

base_GRN = co.data.load_human_promoter_base_GRN()

oracle_links.get_top_links(percentile = 0.1)
oracle_links.TF_degree_dict # the degree of each NODE in each cluster
oracle_links.filtered_links_dict # the top links for each NODE in each cluster
