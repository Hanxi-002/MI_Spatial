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
%matplotlib inline

#%%
print(os.getcwd())
os.chdir("/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/CD68/021424")

adata_oracle = dill.load(open('adata_oracle_CD68_CO_020124.pkl', 'rb'))
oracle_links = dill.load(open('CD68_oracle_links.pkl', 'rb'))

latent_factors = load_SLIDE_res('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/results/SLIDE_Results')
latent_factors.keys()

# delete two attributes of latent_factors
# only keep Z99 as it's the stand alone
del latent_factors['gene_list_Z131']
del latent_factors['gene_list_Z45']

oracle_links.find_TF_overlap_SLIDE(latent_factors = latent_factors, oracle = adata_oracle.oracle)

#%%
for c in range(len(oracle_links.overlap_TF_SLIDE.keys())):
    print(oracle_links.overlap_TF_SLIDE[str(c)])

# the TFs that regulates the latent factors
# 
oracle_links.find_gene_TF_SLIDE(latent_factors, adata_oracle.oracle)
oracle_links.gene_TF_SLIDE
oracle_links.linked_TF_SLIDE

overlap_TF_list = gene_dict_to_list(oracle_links.overlap_TF_SLIDE)
print(overlap_TF_list)
linked_TF_list = oracle_links.threshold_linked_TFs(thresh = 7)
print(linked_TF_list)
len(linked_TF_list)


#%% in silicco perturbation
######################### single TF perturbation ####################
#plt.rcParams["font.family"] = "arial"
plt.rcParams["figure.figsize"] = [6,6]
%config InlineBackend.figure_format = 'retina'
plt.rcParams["savefig.dpi"] = 600

save_folder = "toy_figures"
os.makedirs(save_folder, exist_ok=True)

#%%
# run one goi to check the parameter settings
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



#%%
prob_mat = adata_oracle.oracle.transition_prob
prob_mat.shape
prob_mat.sum(axis = 0)