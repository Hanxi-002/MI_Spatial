setwd("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423")
library(EnhancedVolcano)
library(ggplot2)
library(DESeq2)
library(airway)
library(magrittr)

count <- t(read.csv("x.csv", row.names = 1))
meta <- read.csv("y_index.csv", row.names = 1)
meta$Status <- NULL


dds <- DESeqDataSetFromMatrix(countData = round(count), colData = meta, design = ~y, tidy = FALSE)
dds <- DESeq(dds, betaPrior = FALSE)
res <- results(dds)
res <- lfcShrink(dds, coef=2, res=res)
saveRDS(res, "DESeq_res.RDS")

ggsave(plot = EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                xlim = c(-5, 5),
                pCutoff = 10e-6),
       file = "toy.pdf"
)

setwd('/ocean/projects/cis240075p/hxiao2/MI_Spatial/ER_SLIDE/AllCell/022423')
source('/ocean/projects/cis240075p/hxiao2/MI_Spatial/ER_SLIDE/Find_DEGs_Helper.R')

count <- t(read.csv("x.csv", row.names = 1))
meta <- read.csv("y_index.csv", row.names = 1)
meta$Status <- NULL
DEGs = calculate_DEGs(count, meta, file_name = 'DESeq2_res_20250417.RDS')


custom_color = c("grey30", "grey30", "grey30", "#8FBC8F")
file_name = "volcano_20250417.pdf"
plot_volcano(DEGs = DEGs, custom_color = custom_color, file_name = file_name)
