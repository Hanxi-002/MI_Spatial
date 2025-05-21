library(Seurat)

Fibro_obj <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_Seurat_obj.RDS")
dim(Fibro_obj)

DefaultAssay(Fibro_obj)
grep("TGFB", rownames(Fibro_obj), value = TRUE, ignore.case = TRUE)
FeaturePlot(Fibro_obj, features = 'TGFB1',  cols = c("lightgrey", "blue"), pt.size = 1)
FeaturePlot(Fibro_obj, features = 'TGFB2',  cols = c("lightgrey", "blue"), pt.size = 1)
FeaturePlot(Fibro_obj, features = 'TGFB3',  cols = c("lightgrey", "blue"), pt.size = 1)


expr_data <- GetAssayData(Fibro_obj, slot = "data", assay = "RNA")
cells_to_keep <- colnames(expr_data[, expr_data["TGFB1", ] > 1.5 & expr_data["TGFB2", ] > 1.5 & expr_data["TGFB3", ]>1.5])
AF_subset <- subset(Fibro_obj, cells = cells_to_keep)
print(dim(AF_subset))
saveRDS(AF_subset, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_Seurat_obj.RDS")


################################# Now Get DEGs #################################
Idents(AF_subset) <- 'HF.etiology'
AMI_control_DEG <- FindMarkers(AF_subset, ident.1 = 'AMI', ident.2='Donor', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(AMI_control_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_AMI_vs_Control_DEGs.csv')
ICM_control_DEG <- FindMarkers(AF_subset, ident.1 = 'ICM', ident.2='Donor', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ICM_control_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_ICM_vs_Control_DEGs.csv')
