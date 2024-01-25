library(GeomxTools)
library(dplyr)

Geomx_V3 <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/Geomx_V3.RDS")
HF_data <- subset(x = Geomx_V3, select = Geomx_V3@protocolData@data$Status == "HF")
dim(exprs(HF_data))
if (sum(HF_data@protocolData@data$Status != "HF") != 0) {stop("Control samples are also selected...")}

# extract object for regions that are not CD68
# extract object for regions that are AF and RF
HF_data <- subset(x = HF_data, select = HF_data@phenoData@data$`resting fibroblast` == "FALSE")
if (sum(HF_data@phenoData@data$`resing fibroblast`) != 0) {stop("Not all resting fibroblast are all deleted...")}
x <- exprs(HF_data)
dim(x)

if (sum(row.names(HF_data@protocolData@data) != colnames(x)) != 0) {stop('The sample order of x and y are not the same...')}
y <- as.matrix(as.integer(HF_data@phenoData@data$`CD 68`))
row.names(y) <- row.names(HF_data@protocolData@data)

x <- t(x)
dim(x)

write.csv(x, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/HF/CdVActive/121523/Data/x.csv")
write.csv(y, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/HF/CdVActive/121523/Data/y.csv")


