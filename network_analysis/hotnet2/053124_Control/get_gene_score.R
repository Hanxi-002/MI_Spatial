# getting gene score for hotnet 2 
# cell type specific gene score with non-zero median only in HF samples


Geomx_V3 <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/Geomx_V3.RDS")

main <- function(celltype_obj, geomx_data, ofile) {
  celltype_obj <- subset(celltype_obj, select = celltype_obj@protocolData@data[["Status"]] == "Control")
  median_expr <- celltype_obj@assayData[["q_norm"]]
  row_medians <- apply(median_expr, 1, median)
  write.table(row_medians, sep = "\t", col.names = F, quote = F,
              file = paste0("/ix/djishnu/Hanxi/MI_Spatial/network_analysis/hotnet2/053124/", ofile, ".txt"))
}

#cell_type_name <- "CD 68"
celltype_obj <- subset(Geomx_V3, select = Geomx_V3@phenoData@data$`CD 68` == TRUE)
main(celltype_obj, Geomx_V3, 'CD_68_control')

celltype_obj <- subset(Geomx_V3, select = Geomx_V3@phenoData@data$`resting fibroblast` == TRUE)
main(celltype_obj, Geomx_V3, "resting_fibroblast_control")

celltype_obj <- subset(Geomx_V3, select = Geomx_V3@phenoData@data$`active fibroblast` == TRUE)
main(celltype_obj, Geomx_V3, 'active_fibroblast_control')
