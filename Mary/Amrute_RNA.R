{\rtf1\ansi\ansicpg1252\cocoartf2821
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 Amrute_RNA.R\
\
# MDC 2025 #\
\
# =====================\
## LOAD LIBRARIES ##\
# =====================\
library(Seurat)\
library(Matrix)\
library(ggplot2)\
library(ggrepel)\
library(gridExtra)\
library(sctransform)\
library(tidyverse)\
library(patchwork)\
\
# Lavine = seurat object\
\
#############  SUBSETS ############# \
control <- subset(x = Lavine, subset = HF.etiology == 'Donor')\
dim(control)\
ICM <- subset(x = Lavine, subset = HF.etiology == 'ICM')\
dim(ICM)\
NICM <- subset(x = Lavine, subset = HF.etiology == 'NICM')\
dim(NICM)\
AMI <- subset(x = Lavine, subset = HF.etiology == 'AMI')\
dim(AMI)\
HF <- subset(Lavine, subset = HF.etiology %in% c("ICM", "NICM", "AMI"))\
dim(HF) # should equal the sum of all three experimental groups\
\
# =====================\
## PRE-PROCESSING ##\
# =====================\
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats\
Lavine[["percent.mt"]] <- PercentageFeatureSet(Lavine, pattern = "^MT-")\
# Visualize QC metrics as a violin plot\
VlnPlot(Lavine, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)\
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used\
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.\
\
plot1 <- FeatureScatter(Lavine, feature1 = "nCount_RNA", feature2 = "percent.mt")\
plot2 <- FeatureScatter(Lavine, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")\
plot1 + plot2\
# Lavine <- subset(Lavine, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 15)\
Lavine <- subset(Lavine, subset = nCount_RNA > 500 & nCount_RNA < 25000 & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 15)\
Lavine <- NormalizeData(Lavine)\
Lavine <- FindVariableFeatures(Lavine, selection.method = "vst", nfeatures = 2000)\
# Identify the 10 most highly variable genes\
top10 <- head(VariableFeatures(Lavine), 10)\
all.genes <- rownames(Lavine)\
Lavine <- ScaleData(Lavine, features = all.genes)\
Lavine <- RunPCA(Lavine, features = VariableFeatures(object = Lavine))\
# Examine and visualize PCA results a few different ways\
print(Lavine[["pca"]], dims = 1:10, nfeatures = 5)\
VizDimLoadings(Lavine, dims = 1:2, reduction = "pca")\
DimPlot(Lavine, reduction = "pca") + NoLegend()\
# DimHeatmap(Lavine, dims = 1, cells = 500, balanced = TRUE)\
DimHeatmap(Lavine, dims = 1:15, cells = 500, balanced = TRUE)\
ElbowPlot(Lavine)\
Lavine <- FindNeighbors(Lavine, dims = 1:2)\
Lavine <- FindClusters(Lavine, resolution = 0.1)\
# Look at cluster IDs of the first 5 cells\
head(Idents(Lavine), 5)\
Lavine <- RunUMAP(Lavine, dims = 1:10)\
# note that you can set `label = TRUE` or use the LabelClusters function to help label\
# individual clusters\
DimPlot(Lavine, reduction = "umap")\
\
# =====================\
## DEGs ##\
# =====================\
fibroblasts <- readRDS(file = "Fibroblast_Seurat_obj.RDS")\
macrophages <- readRDS(file = "Macrophage_Seurat_obj.RDS")\
\
# ----- Get inputs for SLIDE ----- #\
# Fibroblasts\
# ICM + CTRL genes\
ICM_CTRL_fib <- subset(fibroblasts, subset = HF.etiology == c("Donor", "ICM"))\
ICM_CTRL_fib_x <- ICM_CTRL_fib@assays$RNA$data\
ICM_CTRL_fib_x <- t(ICM_CTRL_fib_x)\
write.csv(ICM_CTRL_fib_x, "ICM_CTRL_fib_x.csv", row.names = TRUE)\
cell_names <- colnames(ICM_CTRL_fib)\
cell_labels <- ifelse(ICM_CTRL_fib$HF.etiology == "Donor", 0, 1)\
ICM_CTRL_fib_y <- data.frame(Cell = cell_names, Label = cell_labels)\
write.csv(ICM_CTRL_fib_y, "ICM_CTRL_fib_y.csv", row.names = TRUE)\
\
# Fibroblasts \
# AMI + CTRL genes\
AMI_CTRL_fib <- subset(fibroblasts, subset = HF.etiology == c("Donor", "AMI"))\
AMI_CTRL_fib_x <- AMI_CTRL_fib@assays$RNA$data\
AMI_CTRL_fib_x <- t(AMI_CTRL_fib_x)\
write.csv(AMI_CTRL_fib_x, "AMI_CTRL_fib_x.csv", row.names = TRUE)\
cell_names <- colnames(AMI_CTRL_fib)\
cell_labels <- ifelse(AMI_CTRL_fib$HF.etiology == "Donor", 0, 1)\
AMI_CTRL_fib_y <- data.frame(Cell = cell_names, Label = cell_labels)\
write.csv(AMI_CTRL_fib_y, "AMI_CTRL_fib_y.csv", row.names = TRUE)\
\
\
# Macrophages\
# ICM + CTRL genes\
ICM_CTRL_mac <- subset(macrophages, subset = HF.etiology == c("Donor", "ICM"))\
ICM_CTRL_mac_x <- ICM_CTRL_mac@assays$RNA$data\
ICM_CTRL_mac_x <- t(ICM_CTRL_mac_x)\
write.csv(ICM_CTRL_mac_x, "ICM_CTRL_mac_x.csv", row.names = TRUE)\
cell_names <- colnames(ICM_CTRL_mac)\
cell_labels <- ifelse(ICM_CTRL_mac$HF.etiology == "Donor", 0, 1)\
ICM_CTRL_mac_y <- data.frame(Cell = cell_names, Label = cell_labels)\
write.csv(ICM_CTRL_mac_y, "ICM_CTRL_mac_y.csv", row.names = TRUE)\
\
# Macrophages \
# AMI + CTRL genes\
AMI_CTRL_mac <- subset(macrophages, subset = HF.etiology == c("Donor", "AMI"))\
AMI_CTRL_mac_x <- AMI_CTRL_mac@assays$RNA$data\
AMI_CTRL_mac_x <- t(AMI_CTRL_mac_x)\
write.csv(AMI_CTRL_mac_x, "AMI_CTRL_mac_x.csv", row.names = TRUE)\
cell_names <- colnames(AMI_CTRL_mac)\
cell_labels <- ifelse(AMI_CTRL_mac$HF.etiology == "Donor", 0, 1)\
AMI_CTRL_mac_y <- data.frame(Cell = cell_names, Label = cell_labels)\
write.csv(AMI_CTRL_mac_y, "AMI_CTRL_mac_y.csv", row.names = TRUE)\
\
# =====================\
## VARIABLES ##\
# =====================\
##________________________________________________________##\
# ------- SUBSETS by sample HF TYPE or CELL TYPE ------- #\
##________________________________________________________##\
Lavine <- JoinLayers(Lavine)\
# Remove "NA" entries from zones and ctypes\
zones <- unique(Lavine@meta.data$HF.etiology)\
zones <- zones[!is.na(zones) & zones != "NA"]\
ctypes <- unique(Lavine@meta.data$celltype)\
ctypes <- ctypes[!is.na(ctypes) & ctypes != "NA"]\
\
# HF type #\
df_sample_expr <- data.frame(\
  Gene = rownames(Lavine@assays$RNA$counts)\
)\
# zones <- unique(Lavine@meta.data$HF.etiology)\
for (i in seq_along(zones)) \{\
  # Skip if zone is NA\
  if (is.na(zones[i]) || zones[i] == "NA") \{\
    next\
  \}\
  print(paste("Zone:", zones[i]))\
  col_name <- paste0('Expression_', zones[i])\
  df_sample_expr[[col_name]] <- rowMeans(GetAssayData(object = subset(Lavine, subset = HF.etiology == zones[i]),\
                                                      assay = "RNA", \
                                                      layer = "counts"))\
\}\
# cell type #\
df_ctype_expr <- data.frame(\
  Gene = rownames(Lavine@assays$RNA$counts)\
)\
# ctypes <- unique(Lavine@meta.data$celltype)\
for (i in seq_along(ctypes)) \{\
  # Skip if cell type is NA\
  if (is.na(ctypes[i]) || ctypes[i] == "NA") \{\
    next\
  \}\
  print(paste("Cell Type:", ctypes[i]))\
  col_name <- paste0('Expression_', ctypes[i])\
  df_ctype_expr[[col_name]] <- rowMeans(GetAssayData(object = subset(Lavine, subset = celltype == ctypes[i]),\
                                                     assay = "RNA", \
                                                     layer = "counts"))\
\}\
##________________________________________________________##\
# ------- SUBSETS by sample LOCATION & CELL TYPE ------- #\
##________________________________________________________##\
# Make empty data.frame with just gene names\
df_cell_expr <- data.frame(\
  Gene = rownames(Lavine@assays$RNA$counts)\
)\
# Get unique cell types\
cell_type <- unique(Lavine@meta.data$celltype)\
# Loop through zones and cell types\
for (i in seq_along(zones)) \{\
  for (x in seq_along(cell_type)) \{\
    print(paste("Zone:", zones[i], "Cell type:", cell_type[x]))\
    # Create the column name based on the combination of zone and cell type\
    col_name <- paste0('Expression_', zones[i], '_', cell_type[x])\
    # Attempt to subset the Seurat object\
    try(\{\
      subset_obj <- subset(Lavine, subset = HF.etiology == zones[i] & \
                             celltype == cell_type[x])\
      # Check if the subset object contains any cells\
      if (ncol(subset_obj) > 0) \{  # ncol gives the number of cells\
        # Calculate the rowMeans for the subsetted data and assign to the new column\
        df_cell_expr[[col_name]] <- rowMeans(GetAssayData(object = subset_obj, \
                                                          assay = "RNA", \
                                                          layer = "counts"))\
      \} else \{\
        # If no cells found, create an empty column or fill with NA\
        df_cell_expr[[col_name]] <- NA\
        message(paste("No cells found for Zone:", zones[i], "and Cell type:", cell_type[x], "- skipping..."))\
      \}\
    \}, silent = TRUE)  # Silent will suppress the error and continue the loop\
  \}\
\}\
##________________________________________________________##\
# ------- COMBINED SUBSETS by sample LOCATION ------- #\
##________________________________________________________##\
# Initialize a list to store the Seurat object subsets\
seurat_subsets <- list()\
# Loop through each zone and create subsets\
for (i in seq_along(zones)) \{\
  # Skip if the zone is NA\
  if (is.na(zones[i]) || zones[i] == "NA") \{\
    next\
  \}\
  print(paste("Subsetting for zone:", zones[i]))\
  # Subset the Seurat object for the current zone\
  seurat_subsets[[zones[i]]] <- subset(Lavine, subset = HF.etiology == zones[i])\
\}\
# Access the specific subsets\
seurat_subsets$Donor\
# Initialize a list to store the marker results\
markers_list <- list()\
# Exclude "Donor" and "NA" from zones to avoid merging with itself or NA values\
zones_to_process <- zones[!zones %in% c("Donor", "NA")]\
# Loop through each zone and perform the analysis\
for (i in seq_along(zones_to_process)) \{\
  # Skip if the zone is "NA"\
  if (zones_to_process[i] == "NA" || is.na(zones_to_process[i])) \{\
    next\
  \}\
  print(paste("Processing zone:", zones_to_process[i]))\
  # Access the current zone's subset\
  current_zone_subset <- seurat_subsets[[zones_to_process[i]]]\
  # Check if the current zone subset is NULL or empty\
  if (is.null(current_zone_subset) || ncol(current_zone_subset) == 0) \{\
    print(paste("No cells found for zone:", zones_to_process[i], "- skipping."))\
    next  # Skip to the next iteration of the loop\
  \}\
  # Merge Donor with current zone\
  combined_obj <- merge(seurat_subsets$Donor, y = current_zone_subset,\
                        add.cell.ids = c("Donor", zones_to_process[i]))\
  # Join layers to ensure they are ready for analysis\
  combined_obj <- JoinLayers(combined_obj)\
  # Set the identities based on the 'HF.etiology' field\
  Idents(combined_obj) <- "HF.etiology"\
  # Perform differential expression analysis (Donor vs the current zone)\
  markers <- FindMarkers(combined_obj, ident.1 = "Donor", ident.2 = zones_to_process[i])\
  # Store the markers in the list with a named entry\
  markers_list[[paste0("markers_normal_", zones_to_process[i])]] <- markers\
\}\
# access the results for each zone, e.g., markers for normal vs cardiomyopathy\
markers_list$markers_normal_ICM\
##________________________________________________________##\
# ------- COMBINED SUBSETS by sample CELL TYPE ------- #\
##________________________________________________________##\
# initialize a list to store the seurat object subsets\
seurat_subsets_ctypes <- list()\
# Loop through each cell type and create subsets\
for (i in seq_along(ctypes)) \{\
  # Skip if the cell type is "NA" or NA\
  if (ctypes[i] == "NA" || is.na(ctypes[i])) \{\
    next\
  \}\
  print(paste("Subsetting for cell type:", ctypes[i]))\
  # Subset the Seurat object for the current cell type\
  seurat_subsets_ctypes[[ctypes[i]]] <- subset(Lavine, subset = celltype == ctypes[i])\
\}\
# access the specific subsets\
seurat_subsets_ctypes$Fibroblast\
# initialize a list to store the marker results\
# Properly filter out "Fibroblast" and NA values\
ctypes_to_process <- ctypes[!is.na(ctypes) & ctypes != "Fibroblast" & ctypes != "NA"]\
# Initialize a list to store the marker results\
markers_list_ctypes <- list()\
# Loop through each cell type and perform the analysis\
for (i in seq_along(ctypes_to_process)) \{\
  print(paste("Processing cell type:", ctypes_to_process[i]))\
  # Skip any cell type named "NA"\
  if (is.na(ctypes_to_process[i]) || ctypes_to_process[i] == "NA") \{\
    next\
  \}\
  # Merge Fibroblast with the current cell type\
  combined_obj_ctypes <- merge(seurat_subsets_ctypes$Fibroblast, \
                               y = seurat_subsets_ctypes[[ctypes_to_process[i]]],\
                               add.cell.ids = c("Fibroblast", ctypes_to_process[i]))\
  # Join layers to ensure they are ready for analysis\
  combined_obj_ctypes <- JoinLayers(combined_obj_ctypes)\
  # Set the identities based on the 'celltype' field\
  Idents(combined_obj_ctypes) <- "celltype"\
  # Perform differential expression analysis (Fibroblast vs the current cell type)\
  markers <- FindMarkers(combined_obj_ctypes, ident.1 = "Fibroblast", ident.2 = ctypes_to_process[i])\
  # Store the markers in the list with a named entry\
  markers_list_ctypes[[paste0("markers_Fibroblast_", ctypes_to_process[i])]] <- markers\
\}\
\
# access the results for each zone, e.g., markers for CTRL vs IZ\
# markers_list_ctypes$markers_Fibroblast_Myeloid\
\
##_______________________________________________________________##\
# ------- COMBINED SUBSETS by sample LOCATION & CELL TYPE ------- #\
##_______________________________________________________________##\
# Initialize a list for cell type subsets\
seurat_subsets_cell <- list()\
# Loop through each zone and cell type combination to create subsets\
for (zone in zones) \{\
  for (cell in cell_type) \{\
    # Skip "NA" zones and cell types\
    if (zone == "NA" || is.na(zone) || cell == "NA" || is.na(cell)) \{\
      next  # Skip to the next iteration if "NA" zone or cell type\
    \}\
    print(paste("Subsetting for zone:", zone, "and cell type:", cell))\
    # Check if any cells match the criteria, ensuring has_cells is TRUE/FALSE\
    has_cells <- sum(Lavine@meta.data$HF.etiology == zone & \
                       Lavine@meta.data$celltype == cell, na.rm = TRUE) > 0\
    # Proceed only if has_cells is TRUE (not NA)\
    if (!is.na(has_cells) && has_cells) \{\
      subset_data <- subset(Lavine, \
                            subset = HF.etiology == zone & \
                              celltype == cell)\
      # Check if the subset is valid\
      if (ncol(subset_data) > 0 && nrow(subset_data) > 0) \{\
        seurat_subsets_cell[[paste0(zone, "_", cell)]] <- subset_data\
      \} else \{\
        print(paste("Empty subset created for zone:", zone, "and cell type:", cell, "- skipping."))\
      \}\
    \} else \{\
      print(paste("No cells found for zone:", zone, "and cell type:", cell, "- skipping."))\
    \}\
  \}\
\}\
# Now perform comparison between the same cell type in different zones\
markers_list_celltype <- list()\
# Loop through each cell type to perform comparison across zones\
for (cell in cell_type) \{\
  # Skip "NA" cell types\
  if (cell == "NA" || is.na(cell)) \{\
    next  # Skip to the next iteration if "NA" cell type\
  \}\
  for (i in seq_along(zones_to_process)) \{\
    zone <- zones_to_process[i]\
    # Skip "NA" zones\
    if (zone == "NA" || is.na(zone)) \{\
      next  # Skip to the next iteration if "NA" zone\
    \}\
    print(paste("Comparing cell type:", cell, "between normal and", zone))\
    # Define the current zone subset based on the current loop\
    current_zone_subset <- seurat_subsets_cell[[paste0(zone, "_", cell)]]\
    # Access the normal subset\
    normal_subset <- seurat_subsets_cell[[paste0("normal_", cell)]]\
    # Check if both subsets are valid\
    if (is.null(normal_subset) || !validObject(normal_subset) || \
        is.null(current_zone_subset) || !validObject(current_zone_subset) ||\
        ncol(normal_subset) == 0 || \
        ncol(current_zone_subset) == 0) \{\
      print(paste("Invalid subset found for zone:", zone, "or normal subset for cell type:", cell, "- skipping."))\
      next  # Skip to the next iteration of the loop\
    \}\
    # Merge normal with current zone for the specific cell type\
    combined_obj_cell <- merge(normal_subset, \
                               y = current_zone_subset,\
                               add.cell.ids = c(paste0("normal_", cell), paste0(zone, "_", cell)))\
    # Join layers to ensure they are ready for analysis\
    combined_obj_cell <- JoinLayers(combined_obj_cell)\
    # Check if the combined object is valid\
    if (!validObject(combined_obj_cell)) \{\
      print(paste("Combined object for normal and", zone, "for cell type:", cell, "is invalid. Skipping."))\
      next\
    \}\
    # Set identities based on 'disease__ontology_label'\
    Idents(combined_obj_cell) <- "HF.etiology"\
    # Perform differential expression analysis (normal vs current zone for the same cell type)\
    markers_cell <- FindMarkers(combined_obj_cell, ident.1 = "Donor", ident.2 = zone)\
    # Store results\
    markers_list_celltype[[paste0("markers_", cell, "_normal_vs_", zone)]] <- markers_cell\
  \}\
\}\
\
# =====================\
## GENE ANALYSIS ##\
# ===================== \
# Find the olfactory genes\
# Filter genes that start with "OR" but exclude those starting with "ORC" or "ORA"\
or_strings <- grep("^OR(?!C|A)", row.names(Lavine), value = TRUE, perl = TRUE)\
#or_strings <- grep("^OR", row.names(Lavine), value = TRUE)\
or_strings\
# Find the GPR genes\
gpr_strings <- grep("^GPR", row.names(Lavine), value = TRUE)\
gpr_strings\
\
#CCR2 <- grep("^CCR2", row.names(Lavine), value = TRUE)\
#CCR2\
CCR <- grep("^CCR", row.names(Lavine), value = TRUE)\
CCR\
\
# Calculate average expression for each OR or GPR gene across all cells\
expr_OR <- GetAssayData(Lavine, assay = "RNA", layer = "counts")[or_strings, ]\
average_exp_OR<- rowMeans(expr_OR)\
expr_GPR <- GetAssayData(Lavine, assay = "RNA", layer = "counts")[gpr_strings, ]\
average_exp_GPR<- rowMeans(expr_GPR)\
\
# Fibroblasts\
donor_fibroblast_subset <- subset(seurat_subsets_cell$Donor_Fibroblast, subset = HF.etiology == "Donor" & celltype == "Fibroblast")\
AMI_fibroblast_subset <- subset(seurat_subsets_cell$AMI_Fibroblast, subset = HF.etiology == "AMI" & celltype == "Fibroblast")\
fibroblast_donor_counts_OR <- GetAssayData(donor_fibroblast_subset, assay = "RNA", layer = "counts")[or_strings, ]\
fibroblast_AMI_counts_OR <- GetAssayData(AMI_fibroblast_subset, assay = "RNA", layer = "counts")[or_strings, ]\
fibroblast_donor_counts_GPR <- GetAssayData(donor_fibroblast_subset, assay = "RNA", layer = "counts")[gpr_strings, ]\
fibroblast_AMI_counts_GPR <- GetAssayData(AMI_fibroblast_subset, assay = "RNA", layer = "counts")[gpr_strings, ]\
FIB_healthy_rna_data_OR <- rowMeans(GetAssayData(donor_fibroblast_subset, assay = "RNA", layer = "counts")[or_strings, ])\
FIB_AMI_rna_data_OR <- rowMeans(GetAssayData(AMI_fibroblast_subset, assay = "RNA", layer = "counts")[or_strings, ])\
FIB_healthy_rna_data_GPR <- rowMeans(GetAssayData(donor_fibroblast_subset, assay = "RNA", layer = "counts")[gpr_strings, ])\
FIB_AMI_rna_data_GPR <- rowMeans(GetAssayData(AMI_fibroblast_subset, assay = "RNA", layer = "counts")[gpr_strings, ])\
## CCR ##\
fibroblast_donor_counts_CCR <- GetAssayData(donor_fibroblast_subset, assay = "RNA", layer = "counts")[CCR, ]\
fibroblast_AMI_counts_CCR <- GetAssayData(AMI_fibroblast_subset, assay = "RNA", layer = "counts")[CCR, ]\
FIB_healthy_rna_data_CCR <- rowMeans(GetAssayData(donor_fibroblast_subset, assay = "RNA", layer = "counts")[CCR, ])\
FIB_AMI_rna_data_CCR <- rowMeans(GetAssayData(AMI_fibroblast_subset, assay = "RNA", layer = "counts")[CCR, ])\
\
\
# Myeloid\
donor_myeloid_subset <- subset(seurat_subsets_cell$Donor_Myeloid, subset = HF.etiology == "Donor" & celltype == "Myeloid")\
AMI_myeloid_subset <- subset(seurat_subsets_cell$AMI_Myeloid, subset = HF.etiology == "AMI" & celltype == "Myeloid")\
myeloid_donor_counts_OR <- GetAssayData(donor_myeloid_subset, assay = "RNA", layer = "counts")[or_strings, ]\
myeloid_AMI_counts_OR <- GetAssayData(AMI_myeloid_subset, assay = "RNA", layer = "counts")[or_strings, ]\
myeloid_donor_counts_GPR <- GetAssayData(donor_myeloid_subset, assay = "RNA", layer = "counts")[gpr_strings, ]\
myeloid_AMI_counts_GPR <- GetAssayData(AMI_myeloid_subset, assay = "RNA", layer = "counts")[gpr_strings, ]\
MYE_healthy_rna_data_OR <- rowMeans(GetAssayData(donor_myeloid_subset, assay = "RNA", layer = "counts")[or_strings, ])\
MYE_AMI_rna_data_OR <- rowMeans(GetAssayData(AMI_myeloid_subset, assay = "RNA", layer = "counts")[or_strings, ])\
MYE_healthy_rna_data_GPR <- rowMeans(GetAssayData(donor_myeloid_subset, assay = "RNA", layer = "counts")[gpr_strings, ])\
MYE_AMI_rna_data_GPR <- rowMeans(GetAssayData(AMI_myeloid_subset, assay = "RNA", layer = "counts")[gpr_strings, ])\
\
## CCR ##\
myeloid_donor_counts_CCR <- GetAssayData(donor_myeloid_subset, assay = "RNA", layer = "counts")[CCR, ]\
myeloid_AMI_counts_CCR <- GetAssayData(AMI_myeloid_subset, assay = "RNA", layer = "counts")[CCR, ]\
MYE_healthy_rna_data_CCR <- rowMeans(GetAssayData(donor_myeloid_subset, assay = "RNA", layer = "counts")[CCR, ])\
MYE_AMI_rna_data_CCR <- rowMeans(GetAssayData(AMI_myeloid_subset, assay = "RNA", layer = "counts")[CCR, ])\
\
# =====================\
## INITIALIZE DATA FRAMES ##\
# ===================== \
\
## ----- Healthy vs. Disease ----- ##\
# FIBROBLAST #\
\
# Initialize a data frame to store the results\
comparison_results_healthy_disease <- data.frame(Gene = CCR, p_value = NA, avg_log2FC = NA)\
# Loop through each gene to perform the test\
for (Gene in CCR) \{\
  control_expr_fib <- fibroblast_donor_counts_OR[Gene, ]\
  cardio_expr_fib <- fibroblast_AMI_counts_OR[Gene, ]\
  # Perform Wilcoxon test (can replace with t.test for t-test)\
  test_result <- wilcox.test(control_expr_fib, cardio_expr_fib)\
  # Calculate average log2 fold change\
  avg_log2FC <- log2(mean(cardio_expr_fib + 1) / mean(control_expr_fib + 1))\
  # Store results\
  comparison_results_healthy_disease[comparison_results_healthy_disease$Gene == Gene, "p_value"] <- test_result$p.value\
  comparison_results_healthy_disease[comparison_results_healthy_disease$Gene == Gene, "avg_log2FC"] <- avg_log2FC\
\}\
# Adjust p-values for multiple testing (optional, using Benjamini-Hochberg)\
comparison_results_healthy_disease$p_adj <- p.adjust(comparison_results_healthy_disease$p_value, method = "BH")\
# Create a volcano plot\
volc_plot_healthy_disease <- ggplot(comparison_results_healthy_disease, aes(x = avg_log2FC, y = -log10(p_adj))) +\
  geom_point() +\
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +  # P-value threshold line\
  labs(title = "Olfactory Receptor Genes Differential Expression", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +\
  theme_minimal()\
top_genes <- comparison_results_healthy_disease[order(comparison_results_healthy_disease$p_adj)[1:10], ]  # Get top 10 genes based on p-value\
volc_plot_healthy_disease <- volc_plot_healthy_disease +\
  geom_text(data = top_genes, aes(label = Gene), vjust = -0.5, size = 3)\
print(volc_plot_healthy_disease)\
ggsave("volcano_healthy-disease_FIB_OR.png", volc_plot_healthy_disease)\
\
\
# MYELOID #\
\
# Initialize a data frame to store the results\
comparison_results_healthy_disease <- data.frame(Gene = CCR, p_value = NA, avg_log2FC = NA)\
# Loop through each gene to perform the test\
for (Gene in CCR) \{\
  control_expr_macro <- myeloid_donor_counts_OR[Gene, ]\
  cardio_expr_macro <- myeloid_AMI_counts_OR[Gene, ]\
  # Perform Wilcoxon test (can replace with t.test for t-test)\
  test_result <- wilcox.test(control_expr_macro, cardio_expr_macro)\
  # Calculate average log2 fold change\
  avg_log2FC <- log2(mean(cardio_expr_macro + 1) / mean(control_expr_macro + 1))\
  # Store results\
  comparison_results_healthy_disease[comparison_results_healthy_disease$Gene == Gene, "p_value"] <- test_result$p.value\
  comparison_results_healthy_disease[comparison_results_healthy_disease$Gene == Gene, "avg_log2FC"] <- avg_log2FC\
\}\
# Adjust p-values for multiple testing (optional, using Benjamini-Hochberg)\
comparison_results_healthy_disease$p_adj <- p.adjust(comparison_results_healthy_disease$p_value, method = "BH")\
# Create a volcano plot\
volc_plot_healthy_disease <- ggplot(comparison_results_healthy_disease, aes(x = avg_log2FC, y = -log10(p_adj))) +\
  geom_point() +\
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +  # P-value threshold line\
  labs(title = "Olfactory Receptor Genes Differential Expression", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +\
  theme_minimal()\
top_genes <- comparison_results_healthy_disease[order(comparison_results_healthy_disease$p_adj)[1:10], ]  # Get top 10 genes based on p-value\
volc_plot_healthy_disease <- volc_plot_healthy_disease +\
  geom_text(data = top_genes, aes(label = Gene), vjust = -0.5, size = 3)\
print(volc_plot_healthy_disease)\
ggsave("volcano_healthy-disease_MACRO_OR.png", volc_plot_healthy_disease)\
}