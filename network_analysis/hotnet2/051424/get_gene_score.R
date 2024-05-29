# getting gene score for hotnet 2 
# cell type specific gene score with non-zero median only in HF samples


Geomx_V3 <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/Geomx_V3.RDS")

main <- function(cell_type_name, geomx_data) {
  celltype_obj <- subset(Geomx_V3, select = Geomx_V3@phenoData@data$`CD 68` == TRUE)
  celltype_obj <- subset(celltype_obj, select = celltype_obj@protocolData@data[["Status"]] == "HF")
  median_expr <- celltype_obj@assayData[["q_norm"]]
  row_medians <- apply(median_expr, 1, median)
  write.table(row_medians, sep = "\t", col.names = F, quote = F,
              file = paste0("/ix/djishnu/Hanxi/MI_Spatial/network_analysis/hotnet2/051424/", cell_type_name, ".txt"))
}

#cell_type_name <- "CD 68"
main("CD 68", Geomx_V3)
main("resting fibroblast", Geomx_V3)
main("avtive fibroblast", Geomx_V3)
