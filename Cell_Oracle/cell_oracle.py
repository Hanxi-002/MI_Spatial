import os
import sys
import matplotlib.pyplot as plt
import pickle as pkl
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import anndata as ad
import celloracle as co
import scipy.spatial.distance as dist
%matplotlib inline

#%% load data in annData format

# input x should be a dataframe with genes as columns and cells as rows

# creating an object of class adata_oracle that can contain all information
class adata_oracle:

    def __init__(self, file_name):
        # the normalized data csv
        self.norm_data = pd.read_csv(file_name, index_col=0)
        # the adata
        self.adata = ad.AnnData(X=self.norm_data.T)

    def qc_filter(self, min_genes = 200, min_cells = 3):
        sc.pp.filter_cells(self.adata, min_genes=min_genes)
        sc.pp.filter_genes(self.adata, min_cells=min_cells)
        self.adata.var_names = self.norm_data.index
        self.adata.obs_names = self.norm_data.columns
    
    def add_raw_count(self, raw_file_name):
        # read raw count
        raw_data = pd.read_csv(raw_file_name, index_col=0)
        # add raw count to adata
        adata_raw = ad.AnnData(X=raw_data.T)
        self.adata.raw = adata_raw
        self.adata.layers['raw_count'] = self.adata.raw.X.copy()
    
    def mean_centering(self):
        gene_means = np.mean(self.adata.X, axis=0)
        self.adata.X = self.adata.X - gene_means

    def pca(self, n_comps = 50):
        sc.tl.pca(self.adata, n_comps = n_comps, svd_solver='arpack')
        sc.pl.pca(self.adata)
    
    def umap(self):
        sc.tl.umap(self.adata)
    
    def tsne(self):
        sc.tl.tsne(self.adata)
    
    def louvain(self, n_neighbors = 8, n_pcs = 50):
        sc.pp.neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.louvain(self.adata)
    
    def plot_dim_reduction(self, method = "umap", color = "louvain"):
        if color not in self.adata.obs.columns:
            print("color not in adata.obs.columns")
        else:
            if method == "umap":
                sc.pl.umap(self.adata, color=color)
            if method == "tsne":
                sc.pl.tsne(self.adata, color=color)
            else: 
                print("method not supported")
                return None

    def save_adata(self, file_name):
        self.adata.write_h5ad(file_name)


#%% main 
############## load data
print(os.getcwd())

adata_oracle = adata_oracle('all_cell_norm.csv')
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
adata_oracle.save_adata('all_cell_adata.h5ad')





#%% check if the UMAP clustering driven by ROI/AOI location?
# get a object with hf macrophage only
hf_mac = adata_oracle.adata[adata_oracle.adata.obs['Status'] == 'HF']
hf_mac  = hf_mac[hf_mac.obs['cell_type'] == 'CD 68']

sc.pl.umap(hf_mac, color='Status')

# get the meta data for hf macrophage
meta['DCCnames'] = meta['DCCnames'] + ".dcc"
location = meta[meta['DCCnames'].isin(list(hf_mac.obs.index))]
print(location.shape)
location.index = location['DCCnames'] 

# add the x and y coordinates to the hf macrophage object
hf_mac.obs['ROICoordinateX'] = location['ROICoordinateX']
hf_mac.obs['ROICoordinateY'] = location['ROICoordinateY']

# plot the UMAP again with axis labels
plt.plot(hf_mac.obsm['X_umap'][:,0], hf_mac.obsm['X_umap'][:,1], 'o')

# according to the plot above, we can seperate UMAP into 3 clusters.
# for each cluster, subset into a new anndata and plot the locations of each AOI in that cluster
c1 = hf_mac[(hf_mac.obsm['X_umap'][:, 0] <= 0 ) & (hf_mac.obsm['X_umap'][:, 1] <= 0 )]
plt.plot(c1.obsm['X_umap'][:,0], c1.obsm['X_umap'][:,1], 'o')
plt.xlim(-10, 20)
plt.ylim(-10, 12)
c1.shape

c2 = hf_mac[(hf_mac.obsm['X_umap'][:, 0] >= 0 ) & (hf_mac.obsm['X_umap'][:, 1] <= 2.5 )]
plt.plot(c2.obsm['X_umap'][:,0], c2.obsm['X_umap'][:,1], 'o')
plt.xlim(-10, 20)
plt.ylim(-10, 12)
c2.shape

c3 = hf_mac[(hf_mac.obsm['X_umap'][:, 0] >= 0 ) & (hf_mac.obsm['X_umap'][:, 1] >= 2.5 )]
plt.plot(c3.obsm['X_umap'][:,0], c3.obsm['X_umap'][:,1], 'o')
plt.xlim(-10, 20)
plt.ylim(-10, 12)
c3.shape

# plot the AOIs in each cluster with their location
plt.plot(c1.obs['ROICoordinateX'], c1.obs['ROICoordinateY'], 'o')
plt.plot(c2.obs['ROICoordinateX'], c2.obs['ROICoordinateY'], 'o')
plt.plot(c3.obs['ROICoordinateX'], c3.obs['ROICoordinateY'], 'o')
plt.xlim(7000, 35000)
plt.ylim(7000, 40000)
plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')
plt.show()

#%% to make this more quantified, calculat the inter and intra cluster distance for each aoi
c1_dist_array =pd.DataFrame({'x' : c1.obs['ROICoordinateX'], 'y' : c1.obs['ROICoordinateY']})
c1_intra = dist.squareform(dist.pdist(c1_dist_array.to_numpy()))

c2_dist_array =pd.DataFrame({'x' : c2.obs['ROICoordinateX'], 'y' : c2.obs['ROICoordinateY']})
c2_intra = dist.squareform(dist.pdist(c2_dist_array.to_numpy()))

c3_dist_array =pd.DataFrame({'x' : c3.obs['ROICoordinateX'], 'y' : c3.obs['ROICoordinateY']})
c3_intra = dist.squareform(dist.pdist(c3_dist_array.to_numpy()))

plt.imshow(c1_intra, cmap = 'Blues')
plt.colorbar()
plt.xlabel('Cluster 1')
plt.ylabel('Cluster 1')
plt.show()

plt.imshow(c2_intra, cmap = 'Blues')
plt.colorbar()
plt.xlabel('Cluster 2')
plt.ylabel('Cluster 2')
plt.show()

plt.imshow(c3_intra, cmap = 'Blues')
plt.colorbar()
plt.xlabel('Cluster 3')
plt.ylabel('Cluster 3')
plt.show()

# calculate the inter cluster distance
c1_c2_inter = dist.cdist(c1_dist_array.to_numpy(), c2_dist_array.to_numpy())
c1_c3_inter = dist.cdist(c1_dist_array.to_numpy(), c3_dist_array.to_numpy())
c2_c3_inter = dist.cdist(c2_dist_array.to_numpy(), c3_dist_array.to_numpy())

plt.imshow(c1_c2_inter, cmap = 'Blues')
plt.colorbar()
plt.xlabel('Cluster 2')
plt.ylabel('Cluster 1')
plt.show()

plt.imshow(c1_c3_inter, cmap = 'Blues')
plt.colorbar()
plt.xlabel('Cluster 3')
plt.ylabel('Cluster 1')
plt.show()

plt.imshow(c2_c3_inter, cmap = 'Blues')
plt.colorbar()
plt.xlabel('Cluster 3')
plt.ylabel('Cluster 2')
plt.show()

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
#%%

links = co.load_hdf5('links.celloracle.links')

# %%
