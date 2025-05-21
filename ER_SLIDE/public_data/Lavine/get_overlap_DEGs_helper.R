size_match_LF_DEGs <- function(Lavine_DEGs, num_LF_genes, p_val_adj =0.05){
  # input a dataframe of DEGs (csv of saved DEGs from Seurat) and a integer of number of unique LF genes
  # use abs(avg_log2FC) as standard to get top num_LF_genes number of Lavine DEGs
  Lavine_DEGs <- Lavine_DEGs[Lavine_DEGs$p_val_adj <= p_val_adj, ]
  Lavine_DEGs <- Lavine_DEGs[order(abs(Lavine_DEGs$avg_log2FC), decreasing = TRUE), ]
  Lavine_DEGs <- head(Lavine_DEGs, num_LF_genes)
  # return a datafrmae
  return(Lavine_DEGs)
}


get_seurat_DEGs <- function(Lavine_DEGs){
  if ('gene' %in% colnames(Lavine_DEGs)){
    seurat_DEGs <- Lavine_DEGs$gene
  } else {
    seurat_DEGs <- row.names(Lavine_DEGs)
  }
  return(seurat_DEGs)
}


unpack_LF_genes <- function(LF_list){
  # subroutine of get_LF_genes.
  # this function unpakcs each txt and extract gene names
  # LF list is a list contain all txt file of LFs
  LF_genes <- unlist(lapply(LF_list, function(df) df$names))
  return(LF_genes)
}


get_LF_genes <- function(folder_path){
  # read in SLIDE LFs and return a list of LF genes
  folder_path <- folder_path
  file_list <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)
  LF_list <- lapply(file_list, function(file) {
    read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  })
  names(LF_list) <- basename(file_list)
  LF_genes <- unpack_LF_genes(LF_list)
  return(LF_genes)
}



# ######################################## example case ########################################
# ############################# read in SLIDE LFs #############################
# folder_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423/SLIDE_Run_100824/results/0.1_0.5_out_final"
# LF_genes <- get_LF_genes(folder_path)
# num_LF_genes = length(unique(LF_genes))
# 
# ############################# get Lavine DEGs #############################
# # fibroblast allHF v Control
# DEG_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/All_Cell_Type/1_vs_all_DEGs.csv"
# Lavine_MyeFibro_DEG <- read.csv(DEG_path, row.names = 1)
# # size match to num_LFs
# Lavine_MyeFibro_DEG <- size_match_LF_DEGs(Lavine_MyeFibro_DEG, num_LF_genes)
# # get the seurat DEGs
# seurat_DEGs <- get_seurat_DEGs(Lavine_MyeFibro_DEG)
# ############################# get Overlap #############################
# overlap_DEGs = intersect(LF_genes, seurat_DEGs)
# print(length(overlap_DEGs))
# overlap_DEGs = Lavine_MyeFibro_DEG[Lavine_MyeFibro_DEG$gene %in% overlap_DEGs, ]
