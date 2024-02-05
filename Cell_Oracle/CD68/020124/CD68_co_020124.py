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
