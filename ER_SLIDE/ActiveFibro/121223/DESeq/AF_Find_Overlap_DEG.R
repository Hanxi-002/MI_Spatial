source("/ix/djishnu/Hanxi/MI_Spatial/DEG_LF_overlap.R")
setwd("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/ActiveFibro/121223")

DEseq <- readRDS("DESeq/DESeq_res.RDS")

DEG = ExtractDEGs(DEseq_res = DEseq, FC_thresh = 1, padj_thresh = 5)
SLIDE_LFs = ReadAllLFs("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/ActiveFibro/121223/results/SLIDE_Results")
overlap_gene = FindOverlapDEGs(SLIDE_LFs, DEG)
