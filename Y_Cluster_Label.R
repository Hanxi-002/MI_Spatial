# Labels assigned by using GMM.py
# Produce X and Y, Y is spatial cluster label.

library(readxl)

data <- readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Geomx_v2.RDS")
metadata <- read_excel("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/CD68/CD68_Cluster_Label.xlsx")

celltype_data <- subset(data, select = phenoData(data)[["Segment (Name/ Label)"]] == 'CD 68')
celltype_phenodata <- celltype_data@phenoData@data
data_mat <- t(as.matrix(exprs(celltype_data)))

# get Y
y <- c()
for (i in 1:nrow(data_mat)){
  dcc_name <- row.names(data_mat)[i]
  dcc_name <- substr(dcc_name, 1, nchar(dcc_name) - 4)
  dcc_metadata <- metadata[metadata$DCCnames == dcc_name, ]
  if (dim(dcc_metadata)[1] == 0){
    print(dcc_name)
  }
  y <- append(y, dcc_metadata$labels)  
}

Y <- as.data.frame(matrix(nrow = length(y), ncol = 2))
Y['V1'] = metadata$DCCnames
Y['V2'] = y
rownames(Y) <- Y$V1
Y$V1 <- NULL

#output the X and Y
write.csv(data_mat, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/CD68/050123/Data/x.csv')
write.csv(y, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/CD68/050123/Data/y.csv')

