library(Seurat)
library(ggplot2)

################################# Subset Seurat object for Macrophages #################################
Lavine_meyloid <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Myeloid_Seurat_obj.RDS")
dim(Lavine_meyloid)

DefaultAssay(Lavine_meyloid)
FeaturePlot(Lavine_meyloid, features = 'CD68',  cols = c("lightgrey", "blue"), pt.size = 1)
FeaturePlot(Lavine_meyloid, features = 'CD14', cols = c("lightgrey", "blue"), pt.size = 1)


# Extract expression data for GeneX and GeneY
expr_data <- GetAssayData(Lavine_meyloid, slot = "data", assay = "RNA")
geneX_expr <- expr_data["CD68", ]
geneY_expr <- expr_data["CD14", ]

# Add a new metadata column to label cells
Lavine_meyloid$highlight_group <- ifelse(geneX_expr > 2 & geneY_expr > 2, "Highlighted", "Other")

# Plot cells where highlighted cells are distinctly colored
DimPlot(Lavine_meyloid, reduction = "umap", group.by = "highlight_group", 
        cols = c("Other" = "lightgrey", "Highlighted" = "red"), pt.size = 1)


expr_data <- GetAssayData(Lavine_meyloid, slot = "data", assay = "RNA")
cells_to_keep <- colnames(expr_data[, expr_data["CD68", ] > 2 & expr_data["CD14", ] > 2])
Lavine_meyloid_subset <- subset(Lavine_meyloid, cells = cells_to_keep)
print(dim(Lavine_meyloid_subset))
saveRDS(Lavine_meyloid_subset, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Macrophage_Seurat_obj.RDS")

################################# Now Get DEGs #################################
Idents(Lavine_meyloid_subset) <- 'HF.etiology'
AMI_control_DEG <- FindMarkers(Lavine_meyloid_subset, ident.1 = 'AMI', ident.2='Donor', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(AMI_control_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Macrophage_AMI_vs_Control_DEGs.csv')
ICM_control_DEG <- FindMarkers(Lavine_meyloid_subset, ident.1 = 'ICM', ident.2='Donor', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ICM_control_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Macrophage_ICM_vs_Control_DEGs.csv')


###############################################################################
# within Macrophages, get CCR2+ and CCR2- cells
###############################################################################
macrophage_obj <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Macrophage_Seurat_obj.RDS")

DefaultAssay(macrophage_obj)
FeaturePlot(macrophage_obj, features = 'CCR2',  cols = c("lightgrey", "blue"), pt.size = 1)

expr_data <- GetAssayData(macrophage_obj, slot = "data", assay = "RNA")
cells_to_keep <- colnames(expr_data[, expr_data["CCR2", ] > 0]) # only keep macrophages with some CCR2 expression
CCR2_macs = subset(macrophage_obj, cells = cells_to_keep)
dim(CCR2_macs)
# seperate into hi and lo
CCR2_Hi <- colnames(expr_data[, expr_data["CCR2", ] > 1])
CCR2_Lo <- colnames(expr_data[, expr_data["CCR2", ] > 0 & expr_data["CCR2", ] < 1])
print(length(CCR2_Hi) + length(CCR2_Lo))

CCR2_macs@meta.data['CCR2'] <- NA 
CCR2_macs$CCR2[Cells(CCR2_macs) %in% CCR2_Hi] <- "Hi"
CCR2_macs$CCR2[Cells(CCR2_macs) %in% CCR2_Lo] <- "Lo"
saveRDS(CCR2_macs, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Macrophage_CCR2_obj.RDS")

################################# Now Get DEGs #################################
Idents(CCR2_macs) <- 'CCR2'
CCR_Hi_Lo_DEG <- FindMarkers(CCR2_macs, ident.1 = 'Hi', ident.2='Lo', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CCR_Hi_Lo_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Macrophage_CCR2Hi_vs_CCR2Lo_DEGs.csv')


# In CCR2+ macrophages, get AMI vs Control DEGs.
CCR2_macs <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Macrophage_CCR2_obj.RDS")
Idents(CCR2_macs) <- 'CCR2'
CCR2_hi <- subset(CCR2_macs, CCR2 == 'Hi')
Idents(CCR2_hi) <- 'HF.etiology'
AMI_control_DEG <- FindMarkers(CCR2_hi, ident.1 = 'AMI', ident.2='Donor', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(AMI_control_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/CCR2HiMacs_AMI_vs_Control_DEGs.csv')
ICM_control_DEG <- FindMarkers(CCR2_hi, ident.1 = 'ICM', ident.2='Donor', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ICM_control_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/CCR2HiMacs_ICM_vs_Control_DEGs.csv')


