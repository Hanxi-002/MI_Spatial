subset_data <- subset(Geomx_V3, select = Geomx_V3@phenoData@data$`CD 68` == TRUE)
dim(subset_data)

expr <- as.matrix(subset_data@assayData$exprs)
norm <- as.matrix(subset_data@assayData$q_norm)

write.csv(expr, "/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/CD68/020124/CD68.csv")
write.csv(expr, "/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/CD68/020124/CD68_norm.csv")

