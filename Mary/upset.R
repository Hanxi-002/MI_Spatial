{\rtf1\ansi\ansicpg1252\cocoartf2821
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 upset.R\
\
# MDC 2025 #\
\
# =====================\
## LOAD LIBRARIES ##\
# =====================\
library(ggplot2)\
library(ComplexUpset)\
library(patchwork)\
\
\
# ------- FIBROBLASTS -------#\
# ---------------------------#\
\
## PLOT 1 ##\
# Define data frames with lists of DEGs\
IZ_vs_Control_df <- read.csv("./SM_DEGs_fibroblasts_IZ_CTRL_AF.csv")\
BZ_vs_Control_df <- read.csv("./SM_DEGs_fibroblasts_BZ_CTRL_AF.csv")\
AF_LFs_df <- read.csv("./AF_121223_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(IZ_vs_Control_df))\
print(colnames(BZ_vs_Control_df))\
print(colnames(AF_LFs_df))\
\
# Define the features as 'Genes'\
IZ_vs_Control <- as.list(IZ_vs_Control_df[[1]])  # Assuming the genes are in the first column\
BZ_vs_Control <- as.list(BZ_vs_Control_df[[1]])\
AF_LFs <- as.list((AF_LFs_df[[2]]))\
\
# Combine gene lists and create a binary matrix for upset plot\
combine_gene_lists <- function(gene_lists, column_names) \{\
  # Ensure the number of column names matches the number of gene lists\
  if (length(gene_lists) != length(column_names)) \{\
    stop("The number of column names must match the number of gene lists.")\
  \}\
  # Get unique genes from all lists\
  unique_genes <- unique(unlist(gene_lists))\
  # Initialize a data frame with the unique genes\
  df <- data.frame(Genes = unique_genes)\
  # Add columns for each gene list indicating presence (1) or absence (0)\
  for (i in seq_along(gene_lists)) \{\
    df[[column_names[i]]] <- as.integer(df$Genes %in% gene_lists[[i]])\
  \}\
  return(df)\
\}\
# Prepare data for upset plot\
gene_lists <- list(IZ_vs_Control, BZ_vs_Control, AF_LFs)\
colnames <- c("IZ vs. Control (size-matched)","BZ vs. Control (size-matched)","HF vs. Control Latent Factors")\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_1 <- upset(\
  df, c("IZ vs. Control (size-matched)","BZ vs. Control (size-matched)","HF vs. Control Latent Factors"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('IZ vs. Control (size-matched)', 'BZ vs. Control (size-matched)'),\
    c('IZ vs. Control (size-matched)', 'HF vs. Control Latent Factors'),\
    c('BZ vs. Control (size-matched)', 'HF vs. Control Latent Factors')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'IZ vs. Control (size-matched)'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='HF vs. Control Latent Factors',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'BZ vs. Control (size-matched)'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    )\
  )\
) + \
  ggtitle("Size Matched to Latent Factors") + \
  theme(\
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), # Customize title appearance\
    axis.text.y = element_text(size = 20), # Larger y-axis labels\
    axis.title = element_text(size = 18), # Larger axis titles\
    axis.text.y.top = element_text(size = 40),\
    axis.text.y.top = element_text(size = 40)\
  )\
\
\
## PLOT TWO ##\
# Define data frames with lists of DEGs\
IZ_vs_Control_df <- read.csv("./filtered_DEGs_fibroblasts_IZ_CTRL.csv")\
BZ_vs_Control_df <- read.csv("./filtered_DEGs_fibroblasts_BZ_CTRL.csv")\
AF_LFs_df <- read.csv("./AF_121223_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(IZ_vs_Control_df))\
print(colnames(BZ_vs_Control_df))\
print(colnames(AF_LFs_df))\
\
# Define the features as 'Genes'\
IZ_vs_Control <- as.list(IZ_vs_Control_df[[1]])  # Assuming the genes are in the first column\
BZ_vs_Control <- as.list(BZ_vs_Control_df[[1]])\
AF_LFs <- as.list((AF_LFs_df[[2]]))\
\
# Combine gene lists and create a binary matrix for upset plot\
combine_gene_lists <- function(gene_lists, column_names) \{\
  # Ensure the number of column names matches the number of gene lists\
  if (length(gene_lists) != length(column_names)) \{\
    stop("The number of column names must match the number of gene lists.")\
  \}\
  # Get unique genes from all lists\
  unique_genes <- unique(unlist(gene_lists))\
  # Initialize a data frame with the unique genes\
  df <- data.frame(Genes = unique_genes)\
  # Add columns for each gene list indicating presence (1) or absence (0)\
  for (i in seq_along(gene_lists)) \{\
    df[[column_names[i]]] <- as.integer(df$Genes %in% gene_lists[[i]])\
  \}\
  return(df)\
\}\
# Prepare data for upset plot\
# Prepare data for upset plot\
gene_lists <- list(IZ_vs_Control, BZ_vs_Control, AF_LFs)\
colnames <- c("IZ vs. Control","BZ vs. Control","HF vs. Control Latent Factors")\
\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_2 <- upset(\
  df, c("IZ vs. Control","BZ vs. Control","HF vs. Control Latent Factors"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('IZ vs. Control', 'BZ vs. Control'),\
    c('IZ vs. Control', 'HF vs. Control Latent Factors'),\
    c('BZ vs. Control', 'HF vs. Control Latent Factors')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'IZ vs. Control'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='HF vs. Control Latent Factors',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'BZ vs. Control'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    )\
  )\
) + \
  ggtitle("Full Sets") + \
  theme(\
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), # Customize title appearance\
    axis.text.y = element_text(size = 20), # Larger y-axis labels\
    axis.title = element_text(size = 18), # Larger axis titles\
    axis.text.y.top = element_text(size = 40),\
    axis.text.y.top = element_text(size = 40)\
  )\
\
## PLOT THREE ##\
# Define data frames with lists of DEGs\
IZ_vs_Control_df <- read.csv("./SM_DEGs_fibroblasts_IZ_CTRL_AFCCR2.csv")\
BZ_vs_Control_df <- read.csv("./SM_DEGs_fibroblasts_BZ_CTRL_AFCCR2.csv")\
AFCCR2_LFs_df <- read.csv("./AFCCR2_091124_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(IZ_vs_Control_df))\
print(colnames(BZ_vs_Control_df))\
print(colnames(AFCCR2_LFs_df))\
\
# Define the features as 'Genes'\
IZ_vs_Control <- as.list(IZ_vs_Control_df[[1]])  # Assuming the genes are in the first column\
BZ_vs_Control <- as.list(BZ_vs_Control_df[[1]])\
AFCCR2_LFs <- as.list((AFCCR2_LFs_df[[2]]))\
\
# Combine gene lists and create a binary matrix for upset plot\
combine_gene_lists <- function(gene_lists, column_names) \{\
  # Ensure the number of column names matches the number of gene lists\
  if (length(gene_lists) != length(column_names)) \{\
    stop("The number of column names must match the number of gene lists.")\
  \}\
  # Get unique genes from all lists\
  unique_genes <- unique(unlist(gene_lists))\
  # Initialize a data frame with the unique genes\
  df <- data.frame(Genes = unique_genes)\
  # Add columns for each gene list indicating presence (1) or absence (0)\
  for (i in seq_along(gene_lists)) \{\
    df[[column_names[i]]] <- as.integer(df$Genes %in% gene_lists[[i]])\
  \}\
  return(df)\
\}\
# Prepare data for upset plot\
gene_lists <- list(IZ_vs_Control, BZ_vs_Control, AFCCR2_LFs)\
colnames <- c("IZ vs. Control (size-matched)","BZ vs. Control (size-matched)","HF vs. Control Latent Factors (CCR2)")\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_3 <- upset(\
  df, c("IZ vs. Control (size-matched)","BZ vs. Control (size-matched)","HF vs. Control Latent Factors (CCR2)"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('IZ vs. Control (size-matched)', 'BZ vs. Control (size-matched)'),\
    c('IZ vs. Control (size-matched)', 'HF vs. Control Latent Factors (CCR2)'),\
    c('BZ vs. Control (size-matched)', 'HF vs. Control Latent Factors (CCR2)')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors (CCR2)', 'IZ vs. Control (size-matched)'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='HF vs. Control Latent Factors (CCR2)',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors (CCR2)', 'BZ vs. Control (size-matched)'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    )\
  )\
) + \
  ggtitle("Size Matched to Latent Factors") + \
  theme(\
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), # Customize title appearance\
    axis.text.y = element_text(size = 20), # Larger y-axis labels\
    axis.title = element_text(size = 18), # Larger axis titles\
    axis.text.y.top = element_text(size = 40),\
    axis.text.y.top = element_text(size = 40)\
  )\
\
## PLOT FOUR ##\
# Define data frames with lists of DEGs\
IZ_vs_Control_df <- read.csv("./filtered_DEGs_fibroblasts_IZ_CTRL.csv")\
BZ_vs_Control_df <- read.csv("./filtered_DEGs_fibroblasts_BZ_CTRL.csv")\
AFCCR2_LFs_df <- read.csv("./AFCCR2_091124_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(IZ_vs_Control_df))\
print(colnames(BZ_vs_Control_df))\
print(colnames(AFCCR2_LFs_df))\
\
# Define the features as 'Genes'\
IZ_vs_Control <- as.list(IZ_vs_Control_df[[1]])  # Assuming the genes are in the first column\
BZ_vs_Control <- as.list(BZ_vs_Control_df[[1]])\
AFCCR2_LFs <- as.list((AFCCR2_LFs_df[[2]]))\
\
# Combine gene lists and create a binary matrix for upset plot\
combine_gene_lists <- function(gene_lists, column_names) \{\
  # Ensure the number of column names matches the number of gene lists\
  if (length(gene_lists) != length(column_names)) \{\
    stop("The number of column names must match the number of gene lists.")\
  \}\
  # Get unique genes from all lists\
  unique_genes <- unique(unlist(gene_lists))\
  # Initialize a data frame with the unique genes\
  df <- data.frame(Genes = unique_genes)\
  # Add columns for each gene list indicating presence (1) or absence (0)\
  for (i in seq_along(gene_lists)) \{\
    df[[column_names[i]]] <- as.integer(df$Genes %in% gene_lists[[i]])\
  \}\
  return(df)\
\}\
# Prepare data for upset plot\
gene_lists <- list(IZ_vs_Control, BZ_vs_Control, AFCCR2_LFs)\
colnames <- c("IZ vs. Control","BZ vs. Control","HF vs. Control Latent Factors (CCR2)")\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_4 <- upset(\
  df, c("IZ vs. Control","BZ vs. Control","HF vs. Control Latent Factors (CCR2)"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('IZ vs. Control', 'BZ vs. Control'),\
    c('IZ vs. Control', 'HF vs. Control Latent Factors (CCR2)'),\
    c('BZ vs. Control', 'HF vs. Control Latent Factors (CCR2)')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors (CCR2)', 'IZ vs. Control'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='HF vs. Control Latent Factors (CCR2)',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors (CCR2)', 'BZ vs. Control'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    )\
  )\
) + \
  ggtitle("Full Sets") + \
  theme(\
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), # Customize title appearance\
    axis.text.y = element_text(size = 20), # Larger y-axis labels\
    axis.title = element_text(size = 18), # Larger axis titles\
    axis.text.y.top = element_text(size = 40),\
    axis.text.y.top = element_text(size = 40)\
  )\
\
## ----- COMBINE PLOTS ----- ##\
# Combine the plots and save\
combined_plot <- upset_plot_1 / upset_plot_2 | upset_plot_3 / upset_plot_4\
## ------------------------- ##}