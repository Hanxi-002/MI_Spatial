# region by gene in AF and gene in RF (each gene will show up twice)
library(readxl)
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

top_variance <- function(matrix, n) {
  variances <- apply(matrix, 1, var)
  top_indices <- order(variances, decreasing = TRUE)[1:n]
  return(top_indices)
}



# match DCC names in the GeoMX object and the annot file
# the annot file here only contains HF samples
colnames(expr_mat) <- trim_strings(colnames(expr_mat))
# subset the expression matrix
expr_mat <- expr_mat[, colnames(expr_mat) %in% annot$DCCnames]
dim(expr_mat)

# split expression matrix into active and resting
AF_expr_mat <- subset_celltype(annot, "active fibroblast", expr_mat)
RF_expr_mat <- subset_celltype(annot, "resting fibroblast", expr_mat)
dim(AF_expr_mat)
dim(RF_expr_mat)

# get the highes variant gene in each cell type
AF_var_idx <- top_variance(AF_expr_mat, 2500)
RF_var_idx <- top_variance(RF_expr_mat, 2500)
idx <- intersect(AF_var_idx, RF_var_idx)
length(idx)

var_AF_mat <- AF_expr_mat[idx, ]
var_RF_mat <- RF_expr_mat[idx, ]
dim(var_RF_mat)

# use the intersect of the idx to get the final x matrix
rownames(var_AF_mat) <- paste(rownames(var_AF_mat), "_AF", sep = "")
rownames(var_RF_mat) <- paste(rownames(var_RF_mat), "_RF", sep = "")

combined_matrix <- rbind(var_AF_mat, var_RF_mat)
joint_col_names <- paste(colnames(var_AF_mat), colnames(var_RF_mat), sep = "_")
colnames(combined_matrix) <- joint_col_names
dim(combined_matrix)

name = c()
label = c()
# get the y vector
for (r in colnames(combined_matrix)){
  n1 <- sub("_.*", "", r)
  n2 <- sub("^[^_]*_", "", r)
  
  if (annot[annot$DCCnames == n1, ]$SegmentLabel != "active fibroblast"){stop("AF and RF order is wrong...")}
  if (annot[annot$DCCnames == n2, ]$SegmentLabel != "resting fibroblast"){stop("AF and RF order is wrong...")}
  if (annot[annot$DCCnames == n1, ]$AFRF_Ratio != annot[annot$DCCnames == n2, ]$AFRF_Ratio){stop("the 2 samples have different label...")}
  name <- append(name, r)
  label <- append(label, annot[annot$DCCnames == n1, ]$AFRF_Ratio)
}
 
y <- cbind(name, label) 
x <- t(combined_matrix)
dim(x)
dim(y)


#write.csv(x, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/AF_RF_Concat/111523/Data/x.csv')
#write.csv(y, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/AF_RF_Concat/111523/Data/y.csv', row.names = FALSE)

# colnames(norm_mat) <- trim_strings(colnames(norm_mat))
# norm_mat <- norm_mat[, colnames(norm_mat) %in% annot$DCCnames]
# dim(norm_mat)
# 
# 
# for (r in 1:nrow(y)){
#   row <-y[r, ]
#   subset_data <- hf_data[ , (hf_data@protocolData@data[['ScanLabel']] == row["ScanLabel"]) & (hf_data@protocolData@data[['ROILabel']] == as.numeric(row["ROILabel"]))]
#   cnt = 0
#   for (d in subset_data@protocolData@data[["SampleID"]]){
#     if (row.names(y)[r] == d){
#       cnt = cnt + 1
#     }
#   }
#   if (cnt != 1){stop('cannot find the CD68 dcc name in this subsetted data...')}
#   RF_expr <- subset_data[, subset_data@protocolData@data[['SegmentLabel']] == "resting fibroblast"]
#   RF_expr <- as.matrix(exprs(RF_expr))
#   AF_expr <- subset_data[, subset_data@protocolData@data[['SegmentLabel']] == "active fibroblast"]
#   AF_expr <- as.matrix(exprs(AF_expr))
#   rownames(RF_expr) <- paste0(rownames(RF_expr), "_RF")
#   rownames(AF_expr) <- paste0(rownames(AF_expr), "_AF")
#   tmp <- rbind(RF_expr, AF_expr)
#   x[r, ] <- tmp
#   colnames(x) <- rownames(tmp)
# }
# 
# # get the colnames of x for r = 1, then use this to final check after running the entire for loop 
# #check_name <- colnames(x)
# 
# row.names(x) <- row.names(y)
# y <-as.numeric(y[, 'Labels'])
# #write.csv(x, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/060723_hf/Data/x.csv")
# #write.csv(y, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/060723_hf/Data/y.csv")


