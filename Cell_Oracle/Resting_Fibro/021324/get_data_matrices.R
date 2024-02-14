subset_data <- subset(Geomx_V3, select = Geomx_V3@phenoData@data$`resting fibroblast` == TRUE)
dim(subset_data)

expr <- as.matrix(subset_data@assayData$exprs)
norm <- as.matrix(subset_data@assayData$q_norm)

write.csv(expr, "/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/Resting_Fibro/021323/resting_fibro.csv")
write.csv(norm, "/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/Resting_Fibro/021323/resting_fibro_norm.csv")
