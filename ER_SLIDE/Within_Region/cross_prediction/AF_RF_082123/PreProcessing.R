# region by gene in AF and gene in RF (each gene will show up twice)
#setwd("/ix/djishnu/Hanxi/MI_Spatial")
library(readxl)
library(GeomxTools)

data <- readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Geomx_v3.RDS")
y <- read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/cross_prediction/AF_RF_082123/y.csv",row.names = 1)
#metadata <- read_excel("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/CD68/CD68_Cluster_Label.xlsx")


#hf_data <- data[ , data@protocolData@data[['Status']] == "HF"]

x <- data@assayData$exprs
func <- function(x){x = substr(x, 1, nchar(x)-4)}
toy <- lapply(colnames(x), func)
toy <- unlist(toy)
colnames(x) <- toy

hf_AF <- x[, colnames(x) %in% row.names(y)]
dim(hf_AF)


write.csv(t(hf_AF), "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/cross_prediction/AF_RF_082123/x.csv")
