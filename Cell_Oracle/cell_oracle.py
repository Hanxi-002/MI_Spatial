import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import anndata as ad
import celloracle as co
import scipy.spatial.distance as dist
%matplotlib inline

#%% load data in annData format
print(os.getcwd())
data = pd.read_csv('all_cell_norm.csv', index_col=0)
adata = ad.AnnData(X=data.T)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var_names = data.index
adata.obs_names = data.columns

#%% add in the meta data cell type into the annData object
meta = pd.read_excel('../RawData/Final_Annotation_Dutta_sample_info.xlsx', index_col=0)
celltype_df = meta[['SegmentLabel', 'Status']]
celltype_df.index = meta['DCCnames'] + ".dcc"
celltype_df.index.isin(data.columns).sum()
celltype_df = celltype_df.sort_index()
celltype_df.index == data.columns

adata.obs['cell_type'] = celltype_df['SegmentLabel']
adata.obs['status'] = celltype_df['Status']

#%%
#Z score the data first
#sc.pp.scale(adata, max_value=None)  # Scale the data

# try to just mean centering
gene_means = np.mean(adata.X, axis=0)
adata.X = adata.X - gene_means
sc.tl.pca(adata, svd_solver='arpack')  # Run PCA
sc.pl.pca(adata)

sc.pp.neighbors(adata, n_neighbors=8, n_pcs=40)  # Compute the neighborhood graph
sc.tl.louvain(adata)

sc.tl.umap(adata)
sc.tl.tsne(adata)
sc.pl.umap(adata, color='louvain')
#sc.pl.umap(adata, color='cell_type')
sc.pl.umap(adata, color='status')
sc.pl.tsne(adata, color='status')

# umap_coords = pd.DataFrame(adata.obsm["X_umap"])
# umap_coords.index = adata.obs_names

#%% check if the UMAP clustering driven by ROI/AOI location?
# get a object with hf macrophage only
hf_mac = adata[adata.obs['status'] == 'HF']
hf_mac  = hf_mac[hf_mac.obs['cell_type'] == 'CD 68']
sc.pl.umap(hf_mac, color='status')

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
plt.xlim(-10, 15)
plt.ylim(-5, 15)
c1.shape
plt.plot(c1.obs['ROICoordinateX'], c1.obs['ROICoordinateY'], 'o')
plt.xlim(7000, 35000)
plt.ylim(7000, 40000)

c2 = hf_mac[(hf_mac.obsm['X_umap'][:, 0] >= 0 ) & (hf_mac.obsm['X_umap'][:, 1] <= 0 )]
plt.plot(c2.obsm['X_umap'][:,0], c2.obsm['X_umap'][:,1], 'o')
plt.xlim(-10, 15)
plt.ylim(-5, 15)
c2.shape
plt.plot(c2.obs['ROICoordinateX'], c2.obs['ROICoordinateY'], 'o')
plt.xlim(7000, 35000)
plt.ylim(7000, 40000)

c3 = hf_mac[(hf_mac.obsm['X_umap'][:, 0] >= 0 ) & (hf_mac.obsm['X_umap'][:, 1] >= 0 )]
plt.plot(c3.obsm['X_umap'][:,0], c3.obsm['X_umap'][:,1], 'o')
plt.xlim(-10, 15)
plt.ylim(-5, 15)
c3.shape
plt.plot(c3.obs['ROICoordinateX'], c3.obs['ROICoordinateY'], 'o')
plt.xlim(7000, 35000)
plt.ylim(7000, 40000)

# to make this more quantified, calculat the inter and intra cluster distance for each aoi
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

#%% check the clustering quality by comparing with the true cell type in meta data

# get the cluster assignment from louvain clustering
cluster_assignment = pd.DataFrame(adata.obs['louvain'])

# merging
merged_df = pd.merge(celltype_df, cluster_assignment, left_index=True, right_index=True)

for i in range(0, 8):
    print("Cluster " + str(i) + ":")
    subset_df = merged_df[merged_df['louvain'] == str(i)]
    print(subset_df['SegmentLabel'].value_counts())
    print("")


# %%
