library(readxl)
library(GeomxTools)
data <- readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Geomx_v3.RDS")
annot <- read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/HF_Final_Annotation_RFAF_Ratio.csv")

expr_mat <- data@assayData$exprs
#expr_mat <- data@assayData$q_norm

trim_strings <- function(strings) {
  sapply(strings, function(x) substr(x, 1, nchar(x) - 4))
}

subset_celltype <- function(annot, cell_type, x_mat){
  cell_type_Dcc = annot[annot$SegmentLabel == cell_type, ]$DCCnames
  cell_type_x <- x_mat[, colnames(x_mat) %in% cell_type_Dcc]
}


# match DCC names in the GeoMX object and the annot file
# the annot file here only contains HF samples
colnames(expr_mat) <- trim_strings(colnames(expr_mat))
# subset the expression matrix
expr_mat <- expr_mat[, colnames(expr_mat) %in% annot$DCCnames]
dim(expr_mat)
# get macrophages expression matrix
mac_expr_mat <- subset_celltype(annot, "CD 68", expr_mat)

##################get y matrix#################
name = c()
label = c()
# get the y vector
for (r in colnames(mac_expr_mat)){
  if (annot[annot$DCCnames == r, ]$SegmentLabel != "CD 68"){stop("this is not a CD68 AOI...")}
  name <- append(name, r)
  label <- append(label, annot[annot$DCCnames == r, ]$AFRF_Ratio)
}

y <- cbind(name, label) 
x <- t(mac_expr_mat)
dim(x)
dim(y)

#write.csv(x, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/Data/x.csv')
#write.csv(y, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/Data/y.csv', row.names = FALSE)


