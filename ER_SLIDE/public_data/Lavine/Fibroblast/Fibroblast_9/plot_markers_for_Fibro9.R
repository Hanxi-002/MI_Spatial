library(Seurat)
library(ggplot2)
library(dplyr)


Fibroblast_Seurat_obj <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_Seurat_obj.RDS")
AMI_fibro <- subset(Fibroblast_Seurat_obj, HF.etiology == c('AMI', 'ICM'))

plot_gene_expression <- function(seurat_obj, gene, title = NULL) {
  if(is.null(title)) {
    title <- paste(gene, "Expression in Fibroblasts")
  }
  
  FeaturePlot(seurat_obj,
              features = gene,
              reduction = "umap",    
              pt.size = 1,
              order = TRUE,
              cols = c("lightgrey", "purple"),
              max.cutoff = "q95") +
    theme_minimal() +
    ggtitle(title) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
}

# Then you can easily create plots for each gene:
plot_gene_expression(AMI_fibro, "FAP")
plot_gene_expression(AMI_fibro, "THY1")
plot_gene_expression(AMI_fibro, "POSTN")


fibro_9 <- subset(Fibroblast_Seurat_obj, subset = FAP > 1 & POSTN > 2.5)
dim(fibro_9)
unique(fibro_9@meta.data$HF.etiology)

saveRDS(fibro_9, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/fibro_9.RDS')


