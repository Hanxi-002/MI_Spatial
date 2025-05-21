subset_data <- subset(Geomx_V3, select = Geomx_V3@phenoData@data$`active fibroblast` == TRUE)
dim(subset_data)

expr <- as.matrix(subset_data@assayData$exprs)
norm <- as.matrix(subset_data@assayData$q_norm)

write.csv(expr, "/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/Active_Fibro/021323/active_fibro.csv")
write.csv(norm, "/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/Active_Fibro/021323/active_fibro_norm.csv")
