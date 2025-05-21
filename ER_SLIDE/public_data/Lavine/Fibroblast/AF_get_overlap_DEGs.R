# this script calculates the number of genes that are found in the SLIDE LFs and seurat DEGs
# this script is for Activated Fibroblast (subset from entire Lavine RDS) and AF HF v AF control SLIDE analyis

source('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/get_overlap_DEGs_helper.R')

#######################################################################################
# Donor vs all AMI subtype 
#######################################################################################
############################# read in SLIDE LFs #############################
folder_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/ActiveFibro/121223/results/SLIDE_Results"
LF_genes <- get_LF_genes(folder_path)
num_LF_genes = length(unique(LF_genes))

#write.csv(unique(LF_genes),"/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_upset_plot_lists/AF_121223_LFs.csv")
############################# get Lavine DEGs #############################
# fibroblast allHF v Control
DEG_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_AMI_vs_Control_DEGs.csv"
Lavine_Fibro_DEG <- read.csv(DEG_path, row.names = 1)
filtered_DEG <- Lavine_Fibro_DEG[Lavine_Fibro_DEG$p_val_adj <= 0.05, ]


# size match to num_LFs
Lavine_Fibro_DEG <- size_match_LF_DEGs(Lavine_Fibro_DEG, num_LF_genes)
# get the seurat DEGs
seurat_DEGs <- get_seurat_DEGs(Lavine_Fibro_DEG)
############################# get Overlap #############################
overlap_DEGs = intersect(LF_genes, seurat_DEGs)
print(overlap_DEGs)
#overlap_DEGs = Lavine_Fibro_DEG[Lavine_Fibro_DEG$gene %in% overlap_DEGs, ]

# write out fies for upset plot
#write.csv(filtered_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_upset_plot_lists/not_size_matched/AF_AMI_vs_Control_DEGs.csv')
#write.csv(Lavine_Fibro_DEG, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_AMI_vs_Control_DEGs_SizeMatched.csv")
#write.csv(overlap_DEGs, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_AMI_vs_Control_Overlaps.csv")
#######################################################################################
# Donor vs ICM
#######################################################################################

# no need to re-read in LFs
############################# get Lavine DEGs #############################
DEG_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_ICM_vs_Control_DEGs.csv"
Lavine_Fibro_DEG <- read.csv(DEG_path, row.names = 1)
filtered_DEG <- Lavine_Fibro_DEG[Lavine_Fibro_DEG$p_val_adj <= 0.05, ]
# size match to num_LFs
Lavine_Fibro_DEG <- size_match_LF_DEGs(Lavine_Fibro_DEG, num_LF_genes)
# get the seurat DEGs
seurat_DEGs <- get_seurat_DEGs(Lavine_Fibro_DEG)
############################# get Overlap #############################
overlap_DEGs = intersect(LF_genes, seurat_DEGs)
print(length(overlap_DEGs))
overlap_DEGs = Lavine_Fibro_DEG[row.names(Lavine_Fibro_DEG) %in% overlap_DEGs, ]

#write.csv(filtered_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_upset_plot_lists/not_size_matched/AF_ICM_vs_Control_DEGs.csv')
#write.csv(Lavine_Fibro_DEG, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_ICM_vs_Control_DEGs_SizeMatched.csv")
#write.csv(overlap_DEGs, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/AF_ICM_vs_Control_Overlaps.csv")

