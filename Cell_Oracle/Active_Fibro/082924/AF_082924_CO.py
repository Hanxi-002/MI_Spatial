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
'''
Load the objects from the previous cell oracle analysis at CD68/020124.
'''

# Load the data
print(os.getcwd())
os.chdir("/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/Active_Fibro/082924")

# loading the privious adata_oracle object and oracle_links object (both custom)
adata_oracle = dill.load(open('../021324/active_fibroblast_adata_oracle_co_020124.pkl', 'rb'))
oracle_links = dill.load(open('../021324/active_fibroblast_oracle_links_co_020124.pkl', 'rb'))

#%% ###################################################################
'''
Repeat 021324 analysis to save the oracle object post perturbation. 
This was not done in the prvious analysis. 
'''
overlap_TF_list = gene_dict_to_list(oracle_links.overlap_TF_SLIDE)
print(overlap_TF_list)
linked_TF_list = oracle_links.threshold_linked_TFs(thresh = 8)
print(linked_TF_list)
len(linked_TF_list)


#plt.rcParams["font.family"] = "arial"
plt.rcParams["figure.figsize"] = [6,6]
%config InlineBackend.figure_format = 'retina'
plt.rcParams["savefig.dpi"] = 600

save_folder = "toy"
os.makedirs(save_folder, exist_ok=True)

for goi in overlap_TF_list:
    print(goi)
    plot_gene_hist(adata_oracle.oracle.adata, goi, save_folder)

    calc_cell_identity_shifts(oracle=adata_oracle.oracle, goi=goi, save_folder = save_folder, \
                              scale = 40, perturb_value = 0.0)

    get_optimal_min_mass(n_grid=40, oracle=adata_oracle.oracle)

    calc_vector_fields(min_mass = 0.14, goi = goi, oracle=adata_oracle.oracle, \
                       save_folder = save_folder, scale_simulation = 10)

    plot_vector_filed_on_cluster(oracle=adata_oracle.oracle, goi = goi, save_folder = save_folder, scale_simulation=10)    

    calc_cluster_vector_diff(oracle=adata_oracle.oracle, adata=adata_oracle.oracle.adata, goi=goi,\
                              save_folder=save_folder, method = 'leiden')

    dill.dump(adata_oracle, open('perturbed_oracle/021324_CO/overlap_objects/' + goi + '_adata_oracle.pkl', 'wb'))
#%%
#%% ###################################################################
'''
Relax the linked_TF threshold and perform in sillico perturbation. 
'''
print(os.getcwd())
os.chdir("/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/Active_Fibro/082924")

adata_oracle = dill.load(open('../021324/active_fibroblast_adata_oracle_co_020124.pkl', 'rb'))
oracle_links = dill.load(open('../021324/active_fibroblast_oracle_links_co_020124.pkl', 'rb'))

# produce 46 linked TFs
linked_TF_list = oracle_links.threshold_linked_TFs(thresh = 6)
print(linked_TF_list)
len(linked_TF_list)

#plt.rcParams["font.family"] = "arial"
plt.rcParams["figure.figsize"] = [6,6]
%config InlineBackend.figure_format = 'retina'
plt.rcParams["savefig.dpi"] = 600

save_folder = "figures/new_linked_TF"
os.makedirs(save_folder, exist_ok=True)


for goi in linked_TF_list:
    print(goi)
    plot_gene_hist(adata_oracle.oracle.adata, goi, save_folder)

    calc_cell_identity_shifts(oracle=adata_oracle.oracle, goi=goi, save_folder = save_folder, \
                              scale = 40, perturb_value = 0.0)

    get_optimal_min_mass(n_grid=40, oracle=adata_oracle.oracle)

    calc_vector_fields(min_mass = 0.13, goi = goi, oracle=adata_oracle.oracle, \
                       save_folder = save_folder, scale_simulation = 15)

    plot_vector_filed_on_cluster(oracle=adata_oracle.oracle, goi = goi, save_folder = save_folder, scale_simulation=10)    

    calc_cluster_vector_diff(oracle=adata_oracle.oracle, adata=adata_oracle.oracle.adata, goi=goi,\
                              save_folder=save_folder, method = 'louvain')
    
    dill.dump(adata_oracle, open('perturbed_oracle/082924_CO/linked_objects/' + goi + '_adata_oracle.pkl', 'wb'))