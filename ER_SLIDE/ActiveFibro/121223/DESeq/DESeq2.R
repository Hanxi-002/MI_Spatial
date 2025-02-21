setwd("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/ActiveFibro/121223/")
library(EnhancedVolcano)
library(ggplot2)
library(DESeq2)
library(airway)
library(magrittr)

count <- t(read.csv("Data/x.csv", row.names = 1))
meta <- read.csv("Data/y.csv", row.names = 1)
meta$V1 <- factor(meta$V1)

dds <- DESeqDataSetFromMatrix(countData = round(count), colData = meta, design = ~V1, tidy = FALSE)
dds <- DESeq(dds, betaPrior = FALSE)
res <- results(dds)
res <- lfcShrink(dds, coef=2, res=res)
saveRDS(res, "DESeq/DESeq_res.RDS")

ggsave(plot = EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1,
                pCutoff = 10e-6,
                xlim = c(-6, 6)),
       file = "DESeq/Volcano_symm.pdf"
)