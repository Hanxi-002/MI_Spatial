source('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/get_overlap_DEGs_helper.R')

############################# read in SLIDE LFs #############################
folder_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/CD68/121423/results/SLIDE_Results"
LF_genes <- get_LF_genes(folder_path)
num_LF_genes = length(unique(LF_genes))

#write.csv(unique(LF_genes), "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Mac_upset_plot_lists/Mac_121423_LFs.csv")


#######################################################################################
# Donor vs AMI
#######################################################################################
############################# get Lavine DEGs #############################
# myeloid allHF v Control
DEG_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Macrophage_AMI_vs_Control_DEGs.csv"
Lavine_Myeloid_DEG <- read.csv(DEG_path, row.names = 1)
filtered_DEG <- Lavine_Myeloid_DEG[Lavine_Myeloid_DEG$p_val_adj <= 0.05, ]

# size match to num_LFs
Lavine_Myeloid_DEG <- size_match_LF_DEGs(Lavine_Myeloid_DEG, num_LF_genes)
# get the seurat DEGs
seurat_DEGs <- get_seurat_DEGs(Lavine_Myeloid_DEG)
############################# get Overlap #############################
overlap_DEGs = intersect(LF_genes, seurat_DEGs)
print(length(overlap_DEGs))
#overlap_DEGs = Lavine_Myeloid_DEG[row.names(Lavine_Myeloid_DEG) %in% overlap_DEGs, ]

#write.csv(filtered_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Mac_upset_plot_lists/not_size_matched/Mac_AMI_vs_Control_DEGs.csv')
#write.csv(Lavine_Myeloid_DEG, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Mac_upset_plot_lists/size_matched/Mac_AMI_vs_Control_DEGs_SizeMatched.csv")
#write.csv(overlap_DEGs, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Mac_upset_plot_lists/Mac_AMI_vs_Control_Overlaps.csv")


#######################################################################################
################################ Donor vs ICM ################################
#######################################################################################
# no need to re-read in LFs
############################# get Lavine DEGs #############################
# myeloid allHF v Control
DEG_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Macrophage_ICM_vs_Control_DEGs.csv"
Lavine_Myeloid_DEG <- read.csv(DEG_path, row.names = 1)
filtered_DEG <- Lavine_Myeloid_DEG[Lavine_Myeloid_DEG$p_val_adj <= 0.05, ]
# size match to num_LFs
Lavine_Myeloid_DEG <- size_match_LF_DEGs(Lavine_Myeloid_DEG, num_LF_genes)
# get the seurat DEGs
seurat_DEGs <- get_seurat_DEGs(Lavine_Myeloid_DEG)
############################# get Overlap #############################
overlap_DEGs = intersect(LF_genes, seurat_DEGs)
print(length(overlap_DEGs))
overlap_DEGs = Lavine_Myeloid_DEG[row.names(Lavine_Myeloid_DEG) %in% overlap_DEGs, ]

#write.csv(filtered_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Mac_upset_plot_lists/not_size_matched/Mac_ICM_vs_Control_DEGs.csv')
#write.csv(Lavine_Myeloid_DEG, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Mac_upset_plot_lists/size_matched/Mac_ICM_vs_Control_DEGs_SizeMatched.csv")
#write.csv(overlap_DEGs, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Mac_upset_plot_lists/Mac_ICM_vs_Control_Overlaps.csv")


#######################################################################################
# CCR2Hi Macrophages Donor vs AMI
#######################################################################################
# no need to re-read in LFs
############################# get Lavine DEGs #############################
DEG_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/CCR2HiMacs_AMI_vs_Control_DEGs.csv"
Lavine_Myeloid_DEG <- read.csv(DEG_path, row.names = 1)
filtered_DEG <- Lavine_Myeloid_DEG[Lavine_Myeloid_DEG$p_val_adj <= 0.5, ]

Lavine_Myeloid_DEG <- size_match_LF_DEGs(Lavine_Myeloid_DEG, num_LF_genes, p_val_adj = 0.5)
seurat_DEGs <- get_seurat_DEGs(Lavine_Myeloid_DEG) #### stuck here, need to ask Jishnu.

############################# get Overlap #############################
overlap_DEGs = intersect(LF_genes, seurat_DEGs)
print(length(overlap_DEGs))
overlap_DEGs = Lavine_Myeloid_DEG[row.names(Lavine_Myeloid_DEG) %in% overlap_DEGs, ]

#write.csv(filtered_DEG, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Mac_upset_plot_lists/not_size_matched/MacCCR2+_AMI_vs_Control_DEGs.csv')
#write.csv(Lavine_Myeloid_DEG, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Mac_upset_plot_lists/size_matched/MacCCR2+_AMI_vs_Control_DEGs_SizeMatched.csv")
#write.csv(overlap_DEGs, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Mac_upset_plot_lists/MacCCR2+_AMI_vs_Control_Overlaps.csv")


#######################################################################################
# CCR2Hi Macrophages Donor vs ICM
#######################################################################################
# no need to re-read in LFs
############################# get Lavine DEGs #############################
DEG_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/CCR2HiMacs_ICM_vs_Control_DEGs.csv"
Lavine_Myeloid_DEG <- read.csv(DEG_path, row.names = 1)
filtered_DEG <- Lavine_Myeloid_DEG[Lavine_Myeloid_DEG$p_val_adj <= 0.5, ]

Lavine_Myeloid_DEG <- size_match_LF_DEGs(Lavine_Myeloid_DEG, num_LF_genes, p_val_adj = 0.5)
seurat_DEGs <- get_seurat_DEGs(Lavine_Myeloid_DEG) #### stuck here, need to ask Jishnu.

#write.csv(filtered_DEG, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Mac_upset_plot_lists/not_size_matched/MacCCR2+_ICM_vs_Control_DEGs.csv")
#write.csv(Lavine_Myeloid_DEG, "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Mac_upset_plot_lists/MacCCR2+_ICM_vs_Control_DEGs_SizeMatched.csv")
############################# get Overlap #############################
overlap_DEGs = intersect(LF_genes, seurat_DEGs)
print(length(overlap_DEGs))
overlap_DEGs = Lavine_Myeloid_DEG[row.names(Lavine_Myeloid_DEG) %in% overlap_DEGs, ]
# no overlaps

