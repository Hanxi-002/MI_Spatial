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
import scanpy as sc
%matplotlib inline
#%%
############## load data
print(os.getcwd())
os.chdir("/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/all_cell/012923/")

adata_oracle = ao.adata_oracle('all_cell_norm.csv')
adata_oracle.qc_filter()
adata_oracle.add_raw_count('all_cell.csv')


############## add in the meta data cell type into the annData object
meta = pd.read_excel('../RawData/Final_Annotation_Dutta_sample_info.xlsx', index_col=0)
celltype_df = meta[['SegmentLabel', 'Status']]
celltype_df.index = meta['DCCnames'] + ".dcc"
assert celltype_df.index.isin(adata_oracle.norm_data.columns).sum() == adata_oracle.norm_data.shape[1], "length of meta data not match with data"

celltype_df = celltype_df.sort_index()
# make sure that celltype_df is in the same order as data
assert sum(celltype_df.index == adata_oracle.norm_data.columns) == adata_oracle.norm_data.shape[1], "index of meta data not match with data"

adata_oracle.adata.obs['cell_type'] = celltype_df['SegmentLabel']
adata_oracle.adata.obs['Status'] = celltype_df['Status']

############## dimensionality reduction and clustering
adata_oracle.mean_centering()
adata_oracle.pca()
adata_oracle.louvain()
adata_oracle.umap()
adata_oracle.plot_dim_reduction(method = "umap", color = "louvain")
adata_oracle.plot_dim_reduction(method = "umap", color = "Status")


############## save the adata object
#adata_oracle.save_adata('all_cell_adata.h5ad')
#adata = co.load_hdf5('all_cell_adata.h5ad')
#adata = sc.read_h5ad('all_cell_adata.h5ad')
#sc.pl.umap(adata, color="louvain")
#sc.pl.umap(adata, color="Status")

#print(os.getcwd())

#%% check if the UMAP clustering driven by ROI/AOI location?
# # get a object with hf macrophage only
# hf_mac = adata_oracle.adata[adata_oracle.adata.obs['Status'] == 'HF']
# hf_mac  = hf_mac[hf_mac.obs['cell_type'] == 'CD 68']

# sc.pl.umap(hf_mac, color='Status')

# # get the meta data for hf macrophage
# meta['DCCnames'] = meta['DCCnames'] + ".dcc"
# location = meta[meta['DCCnames'].isin(list(hf_mac.obs.index))]
# print(location.shape)
# location.index = location['DCCnames'] 

# # add the x and y coordinates to the hf macrophage object
# hf_mac.obs['ROICoordinateX'] = location['ROICoordinateX']
# hf_mac.obs['ROICoordinateY'] = location['ROICoordinateY']

# # plot the UMAP again with axis labels
# plt.plot(hf_mac.obsm['X_umap'][:,0], hf_mac.obsm['X_umap'][:,1], 'o')

# # according to the plot above, we can seperate UMAP into 3 clusters.
# # for each cluster, subset into a new anndata and plot the locations of each AOI in that cluster
# c1 = hf_mac[(hf_mac.obsm['X_umap'][:, 0] <= 0 ) & (hf_mac.obsm['X_umap'][:, 1] <= 0 )]
# plt.plot(c1.obsm['X_umap'][:,0], c1.obsm['X_umap'][:,1], 'o')
# plt.xlim(-10, 20)
# plt.ylim(-10, 12)
# c1.shape

# c2 = hf_mac[(hf_mac.obsm['X_umap'][:, 0] >= 0 ) & (hf_mac.obsm['X_umap'][:, 1] <= 2.5 )]
# plt.plot(c2.obsm['X_umap'][:,0], c2.obsm['X_umap'][:,1], 'o')
# plt.xlim(-10, 20)
# plt.ylim(-10, 12)
# c2.shape

# c3 = hf_mac[(hf_mac.obsm['X_umap'][:, 0] >= 0 ) & (hf_mac.obsm['X_umap'][:, 1] >= 2.5 )]
# plt.plot(c3.obsm['X_umap'][:,0], c3.obsm['X_umap'][:,1], 'o')
# plt.xlim(-10, 20)
# plt.ylim(-10, 12)
# c3.shape

# # plot the AOIs in each cluster with their location
# plt.plot(c1.obs['ROICoordinateX'], c1.obs['ROICoordinateY'], 'o')
# plt.plot(c2.obs['ROICoordinateX'], c2.obs['ROICoordinateY'], 'o')
# plt.plot(c3.obs['ROICoordinateX'], c3.obs['ROICoordinateY'], 'o')
# plt.xlim(7000, 35000)
# plt.ylim(7000, 40000)
# plt.xlabel('X coordinate')
# plt.ylabel('Y coordinate')
# plt.show()

# #%% to make this more quantified, calculat the inter and intra cluster distance for each aoi
# c1_dist_array =pd.DataFrame({'x' : c1.obs['ROICoordinateX'], 'y' : c1.obs['ROICoordinateY']})
# c1_intra = dist.squareform(dist.pdist(c1_dist_array.to_numpy()))

# c2_dist_array =pd.DataFrame({'x' : c2.obs['ROICoordinateX'], 'y' : c2.obs['ROICoordinateY']})
# c2_intra = dist.squareform(dist.pdist(c2_dist_array.to_numpy()))

# c3_dist_array =pd.DataFrame({'x' : c3.obs['ROICoordinateX'], 'y' : c3.obs['ROICoordinateY']})
# c3_intra = dist.squareform(dist.pdist(c3_dist_array.to_numpy()))

# plt.imshow(c1_intra, cmap = 'Blues')
# plt.colorbar()
# plt.xlabel('Cluster 1')
# plt.ylabel('Cluster 1')
# plt.show()

# plt.imshow(c2_intra, cmap = 'Blues')
# plt.colorbar()
# plt.xlabel('Cluster 2')
# plt.ylabel('Cluster 2')
# plt.show()

# plt.imshow(c3_intra, cmap = 'Blues')
# plt.colorbar()
# plt.xlabel('Cluster 3')
# plt.ylabel('Cluster 3')
# plt.show()

# # calculate the inter cluster distance
# c1_c2_inter = dist.cdist(c1_dist_array.to_numpy(), c2_dist_array.to_numpy())
# c1_c3_inter = dist.cdist(c1_dist_array.to_numpy(), c3_dist_array.to_numpy())
# c2_c3_inter = dist.cdist(c2_dist_array.to_numpy(), c3_dist_array.to_numpy())

# plt.imshow(c1_c2_inter, cmap = 'Blues')
# plt.colorbar()
# plt.xlabel('Cluster 2')
# plt.ylabel('Cluster 1')
# plt.show()

# plt.imshow(c1_c3_inter, cmap = 'Blues')
# plt.colorbar()
# plt.xlabel('Cluster 3')
# plt.ylabel('Cluster 1')
# plt.show()

# plt.imshow(c2_c3_inter, cmap = 'Blues')
# plt.colorbar()
# plt.xlabel('Cluster 3')
# plt.ylabel('Cluster 2')
# plt.show()

#%%
############## moving to cell oracle
oracle = co.Oracle()
print("Metadata columns :", list(adata_oracle.adata.obs.columns))
print("Dimensional reduction: ", list(adata_oracle.adata.obsm.keys()))

#adata_raw = adata.copy()
adata_oracle.adata.X = adata_oracle.adata.layers["raw_count"].copy()

# Instantiate Oracle object.
oracle.import_anndata_as_raw_count(adata=adata_oracle.adata,
                                   cluster_column_name='louvain',
                                   embedding_name='X_umap')

base_GRN = co.data.load_human_promoter_base_GRN()
oracle.import_TF_data(TF_info_matrix=base_GRN)

# Perform PCA
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
oracle.to_hdf5("allcell.celloracle.oracle")


############## calculate GRN
links = oracle.get_links(cluster_name_for_GRN_unit="louvain", alpha=10,
                         verbose_level=10)

links.to_hdf5(file_path="allcell.celloracle.links")

#%%  load the links object and filter.
'''
Load link object.
Filter solely based on p value.
Calculate weighted logp (coef_abs * -logp).

(??)For each cluster, filter out links(edges) that neither source or target nodes are \
    in the human TF list.(???)

    ***** reversed *****

For each TF, rank links(edges) by weighted logp and take top percentile.\
    (save degree and new edges to object.)

'''

oracle_links = oracle_links("allcell.celloracle.links")
oracle_links.filter_pval()
oracle_links.calc_weighted_logp()
oracle_links.links.links_dict['0'].head()

# load the human TF list
# use pandas to read a txt file where each row contains one gene nam
human_TF = pd.read_csv('/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/allTFs_hg38_Scenic.txt', header=None)
oracle_links.filter_human_TF(human_TF=human_TF)

#this might take a few minutes
oracle_links.get_top_links(percentile = 0.1)
oracle_links.TF_degree_dict # the degree of each NODE in each cluster
oracle_links.filtered_links_dict # the top links for each NODE in each cluster

#using dill to save the object
#dill.dump(oracle_links, open('allcell_oracle_links.pkl', 'wb'))


#%% Pull out TFs discovered in SLIDE
# load the SLIDE results
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

folder_path = '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423/SLIDE_Results'  # Replace with your folder path
latent_factors = load_SLIDE_res(folder_path)
print(os.getcwd())
#oracle_links = dill.load(open('allcell_oracle_links.pkl', 'rb'))
#oracle_links.degree_1_SLIDE_overlap(latent_factors = latent_factors, human_TF = human_TF)
#dill.dump(oracle_links, open('allcell_oracle_links.pkl', 'wb'))


#%%
# use dill to load saved link object
oracle_links = dill.load(open('allcell_oracle_links.pkl', 'rb'))
oracle_links.degree_1_overlap
for c in range(len(oracle_links.degree_1_overlap.keys())):
    print(oracle_links.degree_1_overlap[str(c)])

oracle_links.find_gene_TF_SLIDE(latent_factors, human_TF)
oracle_links.gene_TF_SLIDE

for c in range(len(oracle_links.gene_TF_SLIDE.keys())):
    #print(f"Cluster {str(c)} has {len(oracle_links.gene_TF_SLIDE[str(c)])} edges discovered by SLIDE.")
    print(oracle_links.gene_TF_SLIDE[str(c)])



# %%
    

gene_TF = dict()
linked_TF = dict()
# for each louvain cluster
for c in range(len(oracle_links.filtered_links_dict.keys())):
    df = oracle_links.filtered_links_dict[str(c)]
    # for each latent factor
    gene_TF_df = pd.DataFrame()
    linked_TF_df = pd.DataFrame()
    for lf in latent_factors.keys():
        lf_df = latent_factors[lf]
        lf_genes = lf_df['names'].tolist()
        lf_genes = [x for x in lf_genes if x not in human_TF[0].tolist()]
        # find the edges of the latent factor that are also in the GRN
        overlap_df = df[df['source'].isin(lf_genes) | df['target'].isin(lf_genes)]
        overlap_df['latent_factor'] = lf

        # pull out the TF names that link to genes in the LF
        linked_TF_list = list(set(overlap_df['source']).union(set(overlap_df['target'])))
        linked_TF_list = [x for x in linked_TF_list if x in human_TF[0].tolist()]
        temp_df = pd.DataFrame({'cluster':str(c), 'latent_factor': lf, 'linked_TF': linked_TF})
        gene_TF_df = gene_TF_df.append(overlap_df)
        linked_TF_df = linked_TF_df.append(temp_df)
    gene_TF[str(c)] = gene_TF_df
    linked_TF[str(c)] = linked_TF_df
oracle_links.gene_TF_SLIDE = gene_TF
oracle_links.linked_TF_SLIDE = linked_TF
