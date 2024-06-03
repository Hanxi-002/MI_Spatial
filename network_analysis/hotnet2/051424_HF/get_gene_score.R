# getting gene score for hotnet 2 
# cell type specific gene score with non-zero median only in HF samples

Geomx_V3 <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/Geomx_V3.RDS")

main <- function(celltype_obj, geomx_data, ofile) {
  celltype_obj <- subset(celltype_obj, select = celltype_obj@protocolData@data[["Status"]] == "HF")
  median_expr <- celltype_obj@assayData[["q_norm"]]
  row_medians <- apply(median_expr, 1, median)
  write.table(row_medians, sep = "\t", col.names = F, quote = F,
              file = paste0("/ix/djishnu/Hanxi/MI_Spatial/network_analysis/hotnet2/051424/", ofile, ".txt"))
}

#cell_type_name <- "CD 68"
celltype_obj <- subset(Geomx_V3, select = Geomx_V3@phenoData@data$`CD 68` == TRUE)
main(celltype_obj, Geomx_V3, 'CD_68')

celltype_obj <- subset(Geomx_V3, select = Geomx_V3@phenoData@data$`resting fibroblast` == TRUE)
main(celltype_obj, Geomx_V3, "resting_fibroblast")

celltype_obj <- subset(Geomx_V3, select = Geomx_V3@phenoData@data$`active fibroblast` == TRUE)
main(celltype_obj, Geomx_V3, 'active_fibroblast')

##################################################################
# get gene score for CD68 cloe and far away from activated fibroblast

main_dist <- function(celltype_obj, ofile) {
  median_expr <- celltype_obj@assayData[["q_norm"]]
  row_medians <- apply(median_expr, 1, median)
  write.table(row_medians, sep = "\t", col.names = F, quote = F,
              file = paste0("/ix/djishnu/Hanxi/MI_Spatial/network_analysis/hotnet2/051424_HF/", ofile, ".txt"))
}


# read the y label from the SLIDE analyis
y <- as.matrix(read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/Data/y.csv", row.names = 1))
celltype_obj <- subset(Geomx_V3, select = Geomx_V3@phenoData@data$`CD 68` == TRUE)

set_0_y <- y[y==0,]
names(set_0_y)
set_0 <- celltype_obj[, colnames(celltype_obj) %in% names(set_0_y)]
main_dist(set_0, "CD68_far")

set_1_y <- y[y==1,]
names(set_1_y)
set_1 <- celltype_obj[, colnames(celltype_obj) %in% names(set_1_y)]
main_dist(set_1, "CD68_close")

