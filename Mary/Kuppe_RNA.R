{\rtf1\ansi\ansicpg1252\cocoartf2821
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 Kuppe_RNA.R\
\
# MDC 2025 #\
\
# =====================\
## LOAD LIBRARIES ##\
# =====================\
library(Seurat)\
library(future)\
library(ggplot2)\
library(SeuratDisk)\
library(dplyr)\
library(stringr)\
\
# combined_seurat = seurat object\
\
# =====================\
## PRE-PROCESSING ##\
# =====================\
Kramann.Visium.obj <- SCTransform(combined_seurat, assay = "Spatial")\
Kramann.Visium.obj <- RunPCA(combined_seurat, verbose = TRUE)\
Kramann.Visium.obj <- RunPCA(combined_seurat, npcs = 30, features = rownames(combined_seurat))\
Kramann.Visium.obj <- RunUMAP(Kramann.Visium.obj, dims = 1:30)\
Kramann.Visium.obj <- FindNeighbors(Kramann.Visium.obj, reduction = "pca", dims = 1:30)\
Kramann.Visium.obj <- FindClusters(Kramann.Visium.obj, resolution = 0.3)\
\
#############  SUBSETS #############\
# =====================\
## VARIABLES ##\
# =====================\
df_sample_expr <- data.frame(\
  Gene = rownames(Kramann.Visium.obj@assays$Spatial$data)\
)\
zones <- unique(Kramann.Visium.obj@meta.data$major_labl)\
for (i in seq_along(zones)) \{\
  print(paste("Zone:", zones[i]))\
  col_name <- paste0('Expression_', zones[i])\
  df_sample_expr[[col_name]] <- rowMeans(GetAssayData(object = subset(Kramann.Visium.obj, subset = major_labl == zones[i]),\
                                                      assay = "Spatial", \
                                                      layer = "data"))\
\}\
\
df_ctype_expr <- data.frame(\
  Gene = rownames(Kramann.Visium.obj@assays$Spatial$data)\
)\
ctypes <- unique(Kramann.Visium.obj@meta.data$celltype_niche)\
for (i in seq_along(ctypes)) \{\
  print(paste("Cell Type:", ctypes[i]))\
  col_name <- paste0('Expression_', ctypes[i])\
  df_ctype_expr[[col_name]] <- rowMeans(GetAssayData(object = subset(Kramann.Visium.obj, subset = celltype_niche == ctypes[i]),\
                                                     assay = "Spatial", \
                                                     layer = "data"))\
\}\
\
\
# ------- SUBSET by sample GROUP & CELL TYPE ------- #\
# make empty data.frame with just gene names\
df_cell_expr <- data.frame(\
  Gene = rownames(Kramann.Visium.obj@assays$Spatial$data)\
)\
cell_type <- unique(Kramann.Visium.obj@meta.data$celltype_niche)\
for (i in seq_along(zones)) \{\
  for (x in seq_along(cell_type)) \{\
    print(paste("Zone:", zones[i], "Cell type:", cell_type[x]))\
    # Create the column name based on the combination of zone and cell type\
    col_name <- paste0('Expression_', zones[i], '_', cell_type[x])\
    # Use tryCatch to handle errors when subset is empty\
    subset_obj <- tryCatch(\{\
      subset(Kramann.Visium.obj, subset = major_labl == zones[i] & \
               celltype_niche == cell_type[x])\
    \}, error = function(e) \{\
      # Print error message and return NULL if no cells are found\
      print(paste("Error with Zone:", zones[i], "Cell type:", cell_type[x], "- no cells found"))\
      return(NULL)  # Return NULL to skip the rest of the iteration\
    \})\
    # If subset_obj is NULL, skip to the next iteration\
    if (is.null(subset_obj) || ncol(subset_obj) == 0) \{\
      next  # Skip to the next iteration\
    \}\
    # Calculate the rowMeans for the subsetted data and assign to the new column\
    df_cell_expr[[col_name]] <- rowMeans(GetAssayData(object = subset_obj, \
                                                      assay = "Spatial", \
                                                      layer = "data"))\
  \}\
\}\
\
\
## FIND MARKERS ##\
# ------- COMBINED SUBSETS by sample LOCATION ------- #\
# initialize a list to store the seurat object subsets\
seurat_subsets<- list()\
# loop through each zone and create subsets\
for (i in seq_along(zones)) \{\
  print(paste("Subsetting for zone:", zones[i]))\
  # subset the Seurat onject for the current zone\
  seurat_subsets[[zones[i]]] <- subset(Kramann.Visium.obj, subset = major_labl == zones[i])\
\}\
# access the specific subsets\
seurat_subsets$CTRL\
# Initialize a list to store the marker results\
markers_list <- list()\
# Exclude "CTRL" from zones to avoid merging with itself\
zones_to_process <- zones[zones != "CTRL"]\
# Loop through each zone and perform the analysis\
for (i in seq_along(zones_to_process)) \{ \
  print(paste("Processing zone:", zones_to_process[i]))\
  # Merge CTRL with current zone\
  combined_obj <- merge(seurat_subsets$CTRL, y = seurat_subsets[[zones_to_process[i]]], \
                        add.cell.ids = c("CTRL", zones_to_process[i]))\
  # Set the identities based on the 'major_labl' field\
  Idents(combined_obj) <- "major_labl"\
  # Prepare SCT data for FindMarkers\
  combined_obj <- PrepSCTFindMarkers(combined_obj)\
  # Perform differential expression analysis (CTRL vs the current zone)\
  markers <- FindMarkers(combined_obj, ident.1 = "CTRL", ident.2 = zones_to_process[i])\
  # Store the markers in the list with the named entry\
  markers_list[[paste0("markers_CTRL_", zones_to_process[i])]] <- markers\
\}\
# access the results for each zone, e.g., markers for CTRL vs IZ\
markers_list$markers_CTRL_IZ\
\
##_______________________________________________________________##\
# ------- COMBINED SUBSETS by sample LOCATION & CELL TYPE ------- #\
##_______________________________________________________________##\
# Initialize a list to store Seurat object subsets for cell type comparisons\
seurat_subsets_cell <- list()\
# Loop through each zone and cell type combination to create subsets\
for (zone in zones) \{\
  for (cell in cell_type) \{\
    print(paste("Subsetting for zone:", zone, "and cell type:", cell))\
    # Use tryCatch to handle any errors when subsetting\
    tryCatch(\{\
      # Create a subset for each zone and cell type combination\
      subset_obj <- subset(Kramann.Visium.obj, subset = major_labl == zone & \
                             celltype_niche == cell)\
      # Check if the subset has any cells\
      if (ncol(subset_obj) > 0) \{\
        seurat_subsets_cell[[paste0(zone, "_", cell)]] <- subset_obj\
      \} else \{\
        # Print message if no cells found for the combination\
        print(paste("No cells found for Zone:", zone, "and Cell type:", cell))\
      \}\
    \}, error = function(e) \{\
      # Skip and print a message if an error occurs (e.g., empty subset)\
      print(paste("Error subsetting for Zone:", zone, "and Cell type:", cell, "- skipping"))\
    \})\
  \}\
\}\
# Now perform comparison between the same cell type in different zones\
markers_list_celltype <- list()\
# Loop through each cell type to perform comparison across zones\
for (cell in cell_type) \{\
  for (i in seq_along(zones_to_process)) \{\
    zone <- zones_to_process[i]\
    print(paste("Comparing cell type:", cell, "between CTRL and", zone))\
    # Check if both subsets for CTRL and the current zone exist\
    if (paste0("CTRL_", cell) %in% names(seurat_subsets_cell) && \
        paste0(zone, "_", cell) %in% names(seurat_subsets_cell)) \{\
      # Merge CTRL with current zone for the specific cell type\
      combined_obj_cell <- merge(seurat_subsets_cell[[paste0("CTRL_", cell)]], \
                                 y = seurat_subsets_cell[[paste0(zone, "_", cell)]],\
                                 add.cell.ids = c(paste0("CTRL_", cell), paste0(zone, "_", cell)))\
      # Set identities based on 'major_labl'\
      Idents(combined_obj_cell) <- "major_labl"\
      # Prepare the SCT assay for differential expression analysis (on the merged object)\
      combined_obj_cell <- PrepSCTFindMarkers(combined_obj_cell)\
      # Perform differential expression analysis (CTRL vs current zone for the same cell type)\
      markers_cell <- FindMarkers(combined_obj_cell, ident.1 = "CTRL", ident.2 = zone)\
      # Store results\
      markers_list_celltype[[paste0("markers_", cell, "_CTRL_vs_", zone)]] <- markers_cell\
    \} else \{\
      print(paste("Skipping comparison for Cell type:", cell, "between CTRL and", zone, "- one of the subsets is empty"))\
    \}\
  \}\
\}\
\
\
\
# =====================\
## DEGs ##\
# =====================\
# FIBROBLASTS #\
fibroblasts <- subset(Kramann.Visium.obj, subset = celltype_niche == "ctniche_4")\
fibroblasts <- PrepSCTFindMarkers(fibroblasts)\
markers_fibro <- FindMarkers(fibroblasts, ident.1 = "IZ", ident.2 = "CTRL")\
write.csv(markers_fibro, "DEGs_fibroblasts.csv")\
\
# MACROPHAGES #\
macrophages <- subset(Kramann.Visium.obj, subset = celltype_niche == "ctniche_5")\
fibroblasts <- PrepSCTFindMarkers(fibroblasts)\
markers_macro <- FindMarkers(macrophages, ident.1 = "IZ", ident.2 = "CTRL")\
write.csv(markers_macro, "DEGs_macrophages.csv")\
}