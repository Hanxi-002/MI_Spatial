# this script calculates the number of genes that are found in the SLIDE LFs and seurat DEGs
# this script is for Activated Fibroblast (subset from entire Lavine RDS) and CCR2+ v CCR2- SLIDE analyis

source('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/get_overlap_DEGs_helper.R')

#######################################################################################
# Donor vs all AMI subtype 
#######################################################################################
############################# read in SLIDE LFs #############################
folder_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final"
LF_genes <- get_LF_genes(folder_path)
num_LF_genes = length(unique(LF_genes))

#write.csv(unique(LF_genes),"/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_CCR2_upset_plot_lists/AFCCR2_091124_LFs.csv")

############################# get Lavine DEGs #############################
# fibroblast allHF v Control
DEG_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_AMI_vs_Control_DEGs.csv"
Lavine_Fibro_DEG <- read.csv(DEG_path, row.names = 1)
filtered_DEG <- Lavine_Fibro_DEG[Lavine_Fibro_DEG$p_val_adj <= 0.05,]
# size match to num_LFs
Lavine_Fibro_DEG <- size_match_LF_DEGs(Lavine_Fibro_DEG, num_LF_genes)
# get the seurat DEGs
seurat_DEGs <- get_seurat_DEGs(Lavine_Fibro_DEG)
############################# get Overlap #############################
overlap_DEGs = intersect(LF_genes, seurat_DEGs)
print(length(overlap_DEGs))
overlap_DEGs = Lavine_Fibro_DEG[Lavine_Fibro_DEG$gene %in% overlap_DEGs, ]

#write.csv(filtered_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_CCR2_upset_plot_lists/not_size_matched/AFCCR2_AMI_vs_Control_DEGs.csv')
#write.csv(Lavine_Fibro_DEG, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_CCR2_upset_plot_lists/size_matched/AFCCR2_AMI_vs_Control_DEGs_SizeMatched.csv")
#write.csv(overlap_DEGs, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_CCR2_upset_plot_lists/AFCCR2_AMI_vs_Control_Overlaps.csv")
#######################################################################################
# Donor vs all ICM subtype 
#######################################################################################
############################# get Lavine DEGs #############################
# fibroblast allHF v Control
DEG_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_ICM_vs_Control_DEGs.csv"
Lavine_Fibro_DEG <- read.csv(DEG_path, row.names = 1)
filtered_DEG <- Lavine_Fibro_DEG[Lavine_Fibro_DEG$p_val_adj <= 0.05,]
# size match to num_LFs
Lavine_Fibro_DEG <- size_match_LF_DEGs(Lavine_Fibro_DEG, num_LF_genes)
# get the seurat DEGs
seurat_DEGs <- get_seurat_DEGs(Lavine_Fibro_DEG)
############################# get Overlap #############################
overlap_DEGs = intersect(LF_genes, seurat_DEGs)
print(length(overlap_DEGs))
#overlap_DEGs = Lavine_Fibro_DEG[Lavine_Fibro_DEG$gene %in% overlap_DEGs, ]

#write.csv(filtered_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_CCR2_upset_plot_lists/not_size_matched/AFCCR2_ICM_vs_Control_DEGs.csv')
#write.csv(Lavine_Fibro_DEG, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_CCR2_upset_plot_lists/size_matched/AFCCR2_ICM_vs_Control_DEGs_SizeMatched.csv")
#write.csv(overlap_DEGs, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_CCR2_upset_plot_lists/size_matched/AFCCR2_ICM_vs_Control_Overlaps.csv")