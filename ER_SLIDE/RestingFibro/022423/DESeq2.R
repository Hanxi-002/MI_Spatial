library(EnhancedVolcano)
library(ggplot2)
library(DESeq2)
library(airway)
library(magrittr)

getwd()
count <- t(read.csv("x.csv", row.names = 1))
meta <- read.csv("y.csv", row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = round(count), colData = meta, design = ~x, tidy = FALSE)
dds <- DESeq(dds, betaPrior = FALSE)
res <- results(dds)
res <- lfcShrink(dds, coef=2, res=res)

ggsave(plot = EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                pCutoff = 10e-6),
       file = "Volcano.png"
)

