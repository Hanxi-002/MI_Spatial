library(Seurat)
library(ggplot2)


Fibroblast_Seurat_obj <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_Seurat_obj.RDS")
Fibro_9 <- readRDS('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/fibro_9.RDS')
non_fibro_9 <- subset(Fibroblast_Seurat_obj, cells = setdiff(Cells(Fibroblast_Seurat_obj), Cells(Fibro_9)))

dim(Fibroblast_Seurat_obj)
dim(Fibro_9)
dim(non_fibro_9)

non_fibro_9 <- subset(non_fibro_9, HF.etiology == c('Donor', 'AMI', 'ICM'))

unique(non_fibro_9$HF.etiology)

saveRDS(non_fibro_9, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/non_fibro_9/non_fibro_9.RDS')
