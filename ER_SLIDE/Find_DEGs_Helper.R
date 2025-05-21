library(EnhancedVolcano)
library(ggplot2)
library(DESeq2)
library(airway)
library(magrittr)

calculate_DEGs <- function(count, labels, file_name){
  
  cat('Current working directory is:', getwd())
  if (ncol(labels) != 1){stop('labels contain more than 1 column...\n')}
  column_name <- colnames(labels)[1]
  labels <- data.frame(condition = labels[[column_name]], row.names = rownames(labels))
  
  labels$condition <- factor(labels$condition)
  dds <- DESeqDataSetFromMatrix(countData = round(count), colData = labels, design = ~condition, tidy = FALSE)
  dds <- DESeq(dds, betaPrior = FALSE)
  res <- results(dds)
  res <- lfcShrink(dds, coef=2, res=res)
  saveRDS(res, file_name)
  return (res)
}


plot_volcano <- function(DEGs, custom_color, file_name){
  # Create the volcano plot with custom colors
  ggsave(plot = EnhancedVolcano(DEGs,
                                lab = rownames(DEGs),
                                x = 'log2FoldChange',
                                y = 'padj',
                                FCcutoff = 1,
                                xlim = c(-6, 6),
                                pCutoff = 10e-6,
                                colAlpha = 0.7,
                                # Four colors for the four possible regions of the plot
                                col = custom_color,
                                # Names for the legend
                                legendLabels = c("Not sig", "Log2FC", "p-value", "p-value & Log2FC"),
                                # Additional plot customization
                                title = "Volcano Plot",
                                subtitle = "Differential Expression Results",
                                caption = "",
                                labSize = 5,
                                pointSize = 2.0
  ),
  file = file_name
  )
}

# example usage
# setwd("/CD68/121423")
# count <- t(read.csv("Data/x.csv", row.names = 1))
# meta <- read.csv("Data/y.csv", row.names = 1)
# rds_name <- 'DESeq/DESeq_res.RDS'
# DEGs <- calculate_DEGs(count, meta, rds_name)
# 
# custom_color = c("grey30", "grey30", "grey30", "#8FBC8F")
# file_name = "DESeq/toytoy.pdf"
# plot_volcano(DEGs, custom_color, file_name)