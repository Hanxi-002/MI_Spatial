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
############################################ plot z score #############################################

# source('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/calc_Lavine_mac_Z99_score.R')
# 
# seurat_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/fibro_9.RDS"
# er_results_path <- '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/AllLatentFactors.rds'
# z_score_column = 17
# plot_title = "Z17 Score by Condition (Lavine Fibroblast_9)"
# # Run the main script
# results <- main(seurat_path, er_results_path, z_score_column = 17, plot_title)
# write.csv(results$plot$data, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/Z17_Score_by_Condition_(Lavine_Fibroblast_9).csv')
# 
# 
# 
# z_score_column = 37
# plot_title = "Z37 Score by Condition (Lavine Fibroblast_9)"
# # Run the main script
# results <- main(seurat_path, er_results_path, z_score_column = 37, plot_title)
# write.csv(results$plot$data, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/Z37_Score_by_Condition_(Lavine_Fibroblast_9).csv')

