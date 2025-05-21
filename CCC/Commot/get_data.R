
Geomx_V3 <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/Geomx_V3.RDS")
norm_x <- Geomx_V3@assayData$q_norm
raw_x <- Geomx_V3@assayData$exprs

write.table(norm_x, "/ix/djishnu/Hanxi/MI_Spatial/CCC/norm_x.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(raw_x, "/ix/djishnu/Hanxi/MI_Spatial/CCC/raw_x.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
