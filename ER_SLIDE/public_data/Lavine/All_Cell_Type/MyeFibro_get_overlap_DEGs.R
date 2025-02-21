# this script calculates the number of genes that are found in the SLIDE LFs and seurat DEGs
# this script is for MyeFibro (subset from entire Lavine RDS) and all cell type from GeoMX

source('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/get_overlap_DEGs_helper.R')
#library(ggplot2)
#######################################################################################
# Donor vs all HF subtype 
#######################################################################################
############################# read in SLIDE LFs #############################
folder_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423/SLIDE_Run_100824/results/0.1_0.5_out_final"
LF_genes <- get_LF_genes(folder_path)
num_LF_genes = length(unique(LF_genes))

############################# get Lavine DEGs #############################
# fibroblast allHF v Control
DEG_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/All_Cell_Type/1_vs_all_DEGs.csv"
Lavine_MyeFibro_DEG <- read.csv(DEG_path, row.names = 1)
# size match to num_LFs
Lavine_MyeFibro_DEG <- size_match_LF_DEGs(Lavine_MyeFibro_DEG, num_LF_genes)
# get the seurat DEGs
seurat_DEGs <- get_seurat_DEGs(Lavine_MyeFibro_DEG)
############################# get Overlap #############################
overlap_DEGs = intersect(LF_genes, seurat_DEGs)
print(length(overlap_DEGs))
overlap_DEGs = Lavine_MyeFibro_DEG[Lavine_MyeFibro_DEG$gene %in% overlap_DEGs, ]


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
DEG_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/All_Cell_Type/AMI_vs_Control_DEGs.csv"
Lavine_MyeFibro_DEG <- read.csv(DEG_path, row.names = 1)
# size match to num_LFs
Lavine_MyeFibro_DEG <- size_match_LF_DEGs(Lavine_MyeFibro_DEG, num_LF_genes)
# get the seurat DEGs
seurat_DEGs <- get_seurat_DEGs(Lavine_MyeFibro_DEG)
############################# get Overlap #############################
overlap_DEGs = intersect(LF_genes, seurat_DEGs)
print(length(overlap_DEGs))
overlap_DEGs = Lavine_MyeFibro_DEG[row.names(Lavine_MyeFibro_DEG) %in% overlap_DEGs, ]
