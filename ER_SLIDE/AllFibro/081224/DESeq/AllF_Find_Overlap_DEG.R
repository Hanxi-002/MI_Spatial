source("/ix/djishnu/Hanxi/MI_Spatial/DEG_LF_overlap.R")
setwd("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllFibro/081224")

DEseq <- readRDS("DESeq/DESeq_res.RDS")

DEG = ExtractDEGs(DEseq_res = DEseq, FC_thresh = 1, padj_thresh = 5)
SLIDE_LFs = ReadAllLFs("results/SLIDE_results/0.01_0.5_out")
overlap_gene = FindOverlapDEGs(SLIDE_LFs, DEG)

write.table(overlap_gene, "DESeq/overlap_gene.txt")
