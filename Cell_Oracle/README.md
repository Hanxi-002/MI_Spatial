# Cell Oracle Analysis 

### Background


### Links

### Notes
1. Use all cell type in both control and HF to figure out the process. 
2. Input data is the quantile normalized reads from the RData object (GeoMX V3).
3. Louvain clustering. 

### Discussion
1. The cell types doesn't seem like they can be seperated clearly from UMAP or tsne and louvain clustering. 
2. UMAP can clearly seperate HF v Control but not tsne. 
3. Checking if the clusters on UMAP is confounded by the X and Y coordinates of the AOI (specifically for Macrophages in HF samples.) We see some trend of this happening but the location does not fully determine the clustering. 
