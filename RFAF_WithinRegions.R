# region by gene in AF and gene in RF (each gene will show up twice)
library(readxl)
data <- readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Geomx_v3.RDS")

y <- as.matrix(read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/060723_hf/Data/y_w_ROI.csv",
                        row.names = 1))
#metadata <- read_excel("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/CD68/CD68_Cluster_Label.xlsx")


hf_data <- data[ , data@protocolData@data[['Status']] == "HF"]

x <- matrix(nrow = nrow(y), ncol = (dim(data)[1] * 2))
dim(x)


for (r in 1:nrow(y)){
  row <-y[r, ]
  subset_data <- hf_data[ , (hf_data@protocolData@data[['ScanLabel']] == row["ScanLabel"]) & (hf_data@protocolData@data[['ROILabel']] == as.numeric(row["ROILabel"]))]
  cnt = 0
  for (d in subset_data@protocolData@data[["SampleID"]]){
    if (row.names(y)[r] == d){
      cnt = cnt + 1
    }
  }
  if (cnt != 1){stop('cannot find the CD68 dcc name in this subsetted data...')}
  RF_expr <- subset_data[, subset_data@protocolData@data[['SegmentLabel']] == "resting fibroblast"]
  RF_expr <- as.matrix(exprs(RF_expr))
  AF_expr <- subset_data[, subset_data@protocolData@data[['SegmentLabel']] == "active fibroblast"]
  AF_expr <- as.matrix(exprs(AF_expr))
  rownames(RF_expr) <- paste0(rownames(RF_expr), "_RF")
  rownames(AF_expr) <- paste0(rownames(AF_expr), "_AF")
  tmp <- rbind(RF_expr, AF_expr)
  x[r, ] <- tmp
  colnames(x) <- rownames(tmp)
}

# get the colnames of x for r = 1, then use this to final check after running the entire for loop 
#check_name <- colnames(x)

row.names(x) <- row.names(y)
y <-as.numeric(y[, 'Labels'])
write.csv(x, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/060723_hf/Data/x.csv")
write.csv(y, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/060723_hf/Data/y.csv")


