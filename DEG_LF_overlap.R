ExtractDEGs <- function(DEseq_res, FC_thresh, padj_thresh){
  # DESeq_res: the DESeqResults object from running DESeq2
  # The x-axis cutoff on the volcano plot
  # The y_axis cutoff on the volacno plot
  DEseq['log10padj'] <- -log10(DEseq$padj)
  DEGs <- DEseq[(abs(DEseq$log2FoldChange) > FC_thresh) & (DEseq$log10padj > padj_thresh), ]
  cat(length(row.names(DEGs)), "DEGs found.")
  return(DEGs)
}


ReadAllLFs <- function(SLIDE_res_path){
  # read in all the SLIDE results
  folder_path <- SLIDE_res_path
  txt_files <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)
  
  list_of_dataframes <- lapply(txt_files, function(file) {
    read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  })
  names(list_of_dataframes) <- basename(txt_files)  
  return(list_of_dataframes)
}

FindOverlapDEGs <- function(SLIDE_lfs, DEG){
  #SLIDE_lfs: output of ReadAllLFs
  #DEGs: otput of ExtractDEGs
  overlap_gene <- c()
  for (l in SLIDE_lfs){
    overlap_gene <- append(overlap_gene, row.names(DEG)[row.names(DEG) %in% l$names])
    
  }
  overlap_gene <- unique(overlap_gene)
  cat("Out of ", length(row.names(DEG)), "number of DEGs. There are ", length(overlap_gene), "number of overlap DEGs.")
  return(overlap_gene)
}







