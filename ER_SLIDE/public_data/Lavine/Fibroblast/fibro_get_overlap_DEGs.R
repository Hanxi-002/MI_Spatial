# this script calculates the number of genes that are found in the SLIDE LFs and seurat DEGs
# this script is for Fibroblasts(subset from entire Lavine RDS) and all fibroblasts from GeoMX

source('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/get_overlap_DEGs_helper.R')

#######################################################################################
# Donor vs all HF subtype 
#######################################################################################
############################# read in SLIDE LFs #############################
folder_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllFibro/081224/results/SLIDE_results/0.1_0.5_out_final"
LF_genes <- get_LF_genes(folder_path)
num_LF_genes = length(unique(LF_genes))

#write.csv(LF_genes, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AllFibro_upet_plot_lists/AllFibro_081224_LFs.csv")
############################# get Lavine DEGs #############################
# fibroblast allHF v Control
DEG_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/1_vs_all_DEGs.csv"
Lavine_Fibro_DEG <- read.csv(DEG_path, row.names = 1)
# size match to num_LFs
Lavine_Fibro_DEG <- size_match_LF_DEGs(Lavine_Fibro_DEG, num_LF_genes)
# get the seurat DEGs
seurat_DEGs <- get_seurat_DEGs(Lavine_Fibro_DEG)
############################# get Overlap #############################
overlap_DEGs = intersect(LF_genes, seurat_DEGs)
print(length(overlap_DEGs))
overlap_DEGs = Lavine_Fibro_DEG[Lavine_Fibro_DEG$gene %in% overlap_DEGs, ]


# ggplot(data = overlap_DEGs, aes(x = avg_log2FC, y = p_val_adj)) +
#   geom_point(color = "blue", alpha = 0.7) +  # Scatter points
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  # Horizontal line at y=0
#   geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Vertical line at x=0
#   labs(
#     x = "avg_log2FC",
#     y = "-log10(p_val_adj)",
#     title = "Scatter plot of avg_log2FC vs. p_val_adj"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.title = element_text(size = 14),
#     plot.title = element_text(hjust = 0.5, size = 16)
#   )


#######################################################################################
# Donor vs AMI
#######################################################################################

# no need to re-read in LFs
############################# get Lavine DEGs #############################
# fibroblast allHF v Control
DEG_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AMI_vs_Control_DEGs.csv"
Lavine_Fibro_DEG <- read.csv(DEG_path, row.names = 1)
# size match to num_LFs
Lavine_Fibro_DEG <- size_match_LF_DEGs(Lavine_Fibro_DEG, num_LF_genes)
# get the seurat DEGs
seurat_DEGs <- get_seurat_DEGs(Lavine_Fibro_DEG)
############################# get Overlap #############################
overlap_DEGs = intersect(LF_genes, seurat_DEGs)
print(length(overlap_DEGs))
overlap_DEGs = Lavine_Fibro_DEG[row.names(Lavine_Fibro_DEG) %in% overlap_DEGs, ]


#write.csv(Lavine_Fibro_DEG, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AllFibro_upet_plot_lists/AllFibro_AMI_vs_Control_DEGs_SizeMatched.csv")
#write.csv(overlap_DEGs, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AllFibro_upet_plot_lists/AllFibro_AMI_vs_Control_Overlaps.csv")
#######################################################################################
# Donor vs ICM
#######################################################################################

# no need to re-read in LFs
############################# get Lavine DEGs #############################
# fibroblast allHF v Control
DEG_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/ICM_vs_Control_DEGs.csv"
Lavine_Fibro_DEG <- read.csv(DEG_path, row.names = 1)
# size match to num_LFs
Lavine_Fibro_DEG <- size_match_LF_DEGs(Lavine_Fibro_DEG, num_LF_genes)
# get the seurat DEGs
seurat_DEGs <- get_seurat_DEGs(Lavine_Fibro_DEG)
############################# get Overlap #############################
overlap_DEGs = intersect(LF_genes, seurat_DEGs)
print(length(overlap_DEGs))
overlap_DEGs = Lavine_Fibro_DEG[row.names(Lavine_Fibro_DEG) %in% overlap_DEGs, ]

#write.csv(Lavine_Fibro_DEG, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AllFibro_upet_plot_lists/AllFibro_ICM_vs_Control_DEGs_SizeMatched.csv")
#write.csv(overlap_DEGs, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AllFibro_upet_plot_lists/AllFibro_ICM_vs_Control_Overlaps.csv")
