library(Seurat)

Lavine <- readRDS("/ix/djishnu/Mary/MI_data/Lavine_2024/data/expression/Lavine_processed_MDC.rds")
Lavine <- JoinLayers(Lavine)
dim(Lavine)

# let's first see what are all the cell types
print(unique(Lavine@meta.data$celltype))
length(unique(Lavine@meta.data$celltype))

Idents(Lavine) <- 'celltype'
Idents(Lavine)


################################# DEGs for each cell type #################################
celltype_DEGs <- FindAllMarkers(object = Lavine, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(celltype_DEGs, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/celltype_DEGs.csv")

meyloid_DEGs <- celltype_DEGs[celltype_DEGs$cluster == 'Myeloid', ]
fibro_DEGs <- celltype_DEGs[celltype_DEGs$cluster == 'Fibroblast', ]

############################# Myeloid + Fibroblast  DEG calculation #############################
MyeFibro_obj <- subset(x = Lavine, subset = celltype %in% c('Myeloid', 'Fibroblast'))
dim(MyeFibro_obj)
Idents(MyeFibro_obj) <- 'HF.etiology'

# this gives donor(control) against all other HF sub types
all_condition_DEGs <- FindAllMarkers(object = MyeFibro_obj, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_condition_DEGs, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/All_Cell_Type/1_vs_all_DEGs.csv')
# now do AMI with donor
AMI_control_DEG <- FindMarkers(MyeFibro_obj, ident.1 = 'AMI', ident.2='Donor', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(AMI_control_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/All_Cell_Type/AMI_vs_Control_DEGs.csv')

# within meyloid, DEG calculation
myeloid_obj <- subset(x = Lavine, subset = celltype == 'Myeloid')
dim(myeloid_obj)
saveRDS(myeloid_obj, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid_Seurat_obj.RDS")

Idents(myeloid_obj) <- 'HF.etiology'
Idents(myeloid_obj)
# this gives donor(control) against all other HF sub types
all_condition_DEGs <- FindAllMarkers(object = myeloid_obj, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_condition_DEGs, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/1_vs_all_DEGs.csv')
# now do AMI with donor
AMI_control_DEG <- FindMarkers(myeloid_obj, ident.1 = 'AMI', ident.2='Donor', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(AMI_control_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/AMI_vs_Control_DEGs.csv')


# within fibro, DEG calculation
fibro_obj <- subset(x = Lavine, subset = celltype == 'Fibroblast')
dim(fibro_obj)
saveRDS(fibro_obj, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_Seurat_obj.RDS")


#fibro_obj <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_Seurat_obj.RDS")
Idents(fibro_obj) <- 'HF.etiology'
Idents(fibro_obj)
# this gives donor(control) against all other HF sub types
all_condition_DEGs <- FindAllMarkers(object = fibro_obj, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_condition_DEGs, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/1_vs_all_DEGs.csv')
# now do AMI with donor
AMI_control_DEG <- FindMarkers(fibro_obj, ident.1 = 'AMI', ident.2='Donor', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(AMI_control_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AMI_vs_Control_DEGs.csv')
# now do ICM with donor
ICM_control_DEG <- FindMarkers(fibro_obj, ident.1 = 'ICM', ident.2='Donor', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ICM_control_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/ICM_vs_Control_DEGs.csv')
