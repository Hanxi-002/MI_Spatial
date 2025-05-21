import sys
sys.path.append("Cell_Oracle/COAnalyses")
import os
# set working directory

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

# read in the adata object
adata = co.load_hdf5("Cell_Oracle/all_cell/012924/allcell.celloracle.oracle")
adata = adata.adata


legend_elements = None
n_cell_types = len(adata.obs['cell_type'].unique())
fig, axes = plt.subplots(1, n_cell_types, figsize=(12, 4)) 
for i, cell in enumerate(adata.obs['cell_type'].unique()):
    cell_adata = adata[adata.obs['cell_type'] == cell]
    sc.pl.umap(cell_adata, color='Status', palette = ['#0000FF', '#FF0000'], \
               title=cell, show=False, size = 100, ax=axes[i])


plt.tight_layout()
plt.savefig("Cell_Oracle/all_cell/012924/umap_cell_types_v2.pdf", dpi=300)
plt.show()