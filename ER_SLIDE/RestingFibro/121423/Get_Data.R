library(dplyr)

Geomx_V3 <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/Geomx_V3.RDS")
AF_data <- subset(Geomx_V3, select = Geomx_V3@phenoData@data$`resting fibroblast` == TRUE)
AF_mat <- as.matrix(exprs(AF_data))
dim(AF_mat)

# check if the metadata table is in the same order as the expression table
sum(colnames(AF_mat) != row.names(AF_data@protocolData@data))
y <- as.matrix(recode(AF_data@protocolData@data[["Status"]], 'Control' = 0, 'HF' = 1))
row.names(y) <- row.names(AF_data@protocolData@data)
if(dim(AF_mat)[2] != length(y)) {stop("length doesn't match between x and y...")}
AF_mat <- t(AF_mat)

row.names(AF_mat)
old_x <- read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/RestingFibro/022423/x.csv", row.names = 1)
sum(row.names(AF_mat) != row.names(old_x))

write.csv(AF_mat, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/RestingFibro/122423/Data/x.csv')
write.csv(y, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/RestingFibro/122423/Data/y.csv')

check_x <- as.matrix(read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/CD68/121423/Data/x.csv', row.names = 1))