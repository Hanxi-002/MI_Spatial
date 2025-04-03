{\rtf1\ansi\ansicpg1252\cocoartf2821
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 upset_Amrute.R\
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
########### ICM/AMI-CTRL vs AF/CCR2 ############\
# ------- FIBROBLASTS -------#\
# ---------------------------#\
\
## PLOT 1 ##\
# Define data frames with lists of DEGs\
AF_AMI_vs_Control_df <- read.csv("./AF_AMI_vs_Control_DEGs_SizeMatched.csv")\
AF_ICM_vs_Control_df <- read.csv("./AF_ICM_vs_Control_DEGs_SizeMatched.csv")\
AF_LFs_df <- read.csv("./AF_121223_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(AF_AMI_vs_Control_df))\
print(colnames(AF_ICM_vs_Control_df))\
print(colnames(AF_LFs_df))\
\
# Define the features as 'Genes'\
AF_AMI_vs_Control <- as.list(AF_AMI_vs_Control_df[[1]])  # Assuming the genes are in the first column\
AF_ICM_vs_Control <- as.list(AF_ICM_vs_Control_df[[1]])\
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
gene_lists <- list(AF_AMI_vs_Control, AF_ICM_vs_Control, AF_LFs)\
colnames <- c("AMI vs. Control (size-matched)","ICM vs. Control (size-matched)","HF vs. Control Latent Factors")\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_1 <- upset(\
  df, c("AMI vs. Control (size-matched)","ICM vs. Control (size-matched)","HF vs. Control Latent Factors"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('AMI vs. Control (size-matched)', 'ICM vs. Control (size-matched)'),\
    c('AMI vs. Control (size-matched)', 'HF vs. Control Latent Factors'),\
    c('ICM vs. Control (size-matched)', 'HF vs. Control Latent Factors')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'AMI vs. Control (size-matched)'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='HF vs. Control Latent Factors',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'ICM vs. Control (size-matched)'),\
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
\
## PLOT TWO ##\
# Define data frames with lists of DEGs\
setwd("/ix/djishnu/Mary/MI_data/Lavine_2024/plots+figures/AF_upset_plot_lists/not_size_matched")\
AF_AMI_vs_Control_df <- read.csv("./AF_AMI_vs_Control_DEGs.csv")\
AF_ICM_vs_Control_df <- read.csv("./AF_ICM_vs_Control_DEGs.csv")\
AF_LFs_df <- read.csv("./AF_121223_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(AF_AMI_vs_Control_df))\
print(colnames(AF_ICM_vs_Control_df))\
print(colnames(AF_LFs_df))\
\
# Define the features as 'Genes'\
AF_AMI_vs_Control <- as.list(AF_AMI_vs_Control_df[[1]])  # Assuming the genes are in the first column\
AF_ICM_vs_Control <- as.list(AF_ICM_vs_Control_df[[1]])\
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
gene_lists <- list(AF_AMI_vs_Control, AF_ICM_vs_Control, AF_LFs)\
colnames <- c("AMI vs. Control","ICM vs. Control","HF vs. Control Latent Factors")\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_2 <- upset(\
  df, c("AMI vs. Control","ICM vs. Control","HF vs. Control Latent Factors"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('AMI vs. Control', 'ICM vs. Control'),\
    c('AMI vs. Control', 'HF vs. Control Latent Factors'),\
    c('ICM vs. Control', 'HF vs. Control Latent Factors')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'AMI vs. Control'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='HF vs. Control Latent Factors',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'ICM vs. Control'),\
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
setwd("/ix/djishnu/Mary/MI_data/Lavine_2024/plots+figures/AF_CCR2_upset_plot_lists/size_matched")\
AFCCR2_AMI_vs_Control_df <- read.csv("./AFCCR2_AMI_vs_Control_DEGs_SizeMatched.csv")\
AFCCR2_ICM_vs_Control_df <- read.csv("./AFCCR2_ICM_vs_Control_DEGs_SizeMatched.csv")\
AFCCR2_LFs_df <- read.csv("./AFCCR2_091124_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(AFCCR2_AMI_vs_Control_df))\
print(colnames(AFCCR2_ICM_vs_Control_df))\
print(colnames(AFCCR2_LFs_df))\
\
# Define the features as 'Genes'\
AFCCR2_AMI_vs_Control <- as.list(AFCCR2_AMI_vs_Control_df[[1]])  # Assuming the genes are in the first column\
AFCCR2_ICM_vs_Control <- as.list(AFCCR2_ICM_vs_Control_df[[1]])\
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
gene_lists <- list(AFCCR2_AMI_vs_Control, AFCCR2_ICM_vs_Control, AFCCR2_LFs)\
colnames <- c("CCR2 AMI vs. Control (size-matched)","CCR2 ICM vs. Control (size-matched)","HF vs. Control Latent Factors (CCR2)")\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_3 <- upset(\
  df, c("CCR2 AMI vs. Control (size-matched)","CCR2 ICM vs. Control (size-matched)","HF vs. Control Latent Factors (CCR2)"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('CCR2 AMI vs. Control (size-matched)', 'CCR2 ICM vs. Control (size-matched)'),\
    c('CCR2 AMI vs. Control (size-matched)', 'HF vs. Control Latent Factors (CCR2)'),\
    c('CCR2 ICM vs. Control (size-matched)', 'HF vs. Control Latent Factors (CCR2)')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors (CCR2)', 'CCR2 AMI vs. Control (size-matched)'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='HF vs. Control Latent Factors (CCR2)',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors (CCR2)', 'CCR2 ICM vs. Control (size-matched)'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    )\
  )\
) + \
  ggtitle("Size-Matched to Latent Factors") + \
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
setwd("/ix/djishnu/Mary/MI_data/Lavine_2024/plots+figures/AF_CCR2_upset_plot_lists/not_size_matched")\
AFCCR2_AMI_vs_Control_df <- read.csv("./AFCCR2_AMI_vs_Control_DEGs.csv")\
AFCCR2_ICM_vs_Control_df <- read.csv("./AFCCR2_ICM_vs_Control_DEGs.csv")\
AFCCR2_LFs_df <- read.csv("./AFCCR2_091124_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(AFCCR2_AMI_vs_Control_df))\
print(colnames(AFCCR2_ICM_vs_Control_df))\
print(colnames(AFCCR2_LFs_df))\
\
# Define the features as 'Genes'\
AFCCR2_AMI_vs_Control <- as.list(AFCCR2_AMI_vs_Control_df[[1]])  # Assuming the genes are in the first column\
AFCCR2_ICM_vs_Control <- as.list(AFCCR2_ICM_vs_Control_df[[1]])\
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
gene_lists <- list(AFCCR2_AMI_vs_Control, AFCCR2_ICM_vs_Control, AFCCR2_LFs)\
colnames <- c("CCR2 AMI vs. Control","CCR2 ICM vs. Control","HF vs. Control Latent Factors (CCR2)")\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_4 <- upset(\
  df, c("CCR2 AMI vs. Control","CCR2 ICM vs. Control","HF vs. Control Latent Factors (CCR2)"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('CCR2 AMI vs. Control', 'CCR2 ICM vs. Control'),\
    c('CCR2 AMI vs. Control', 'HF vs. Control Latent Factors (CCR2)'),\
    c('CCR2 ICM vs. Control', 'HF vs. Control Latent Factors (CCR2)')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors (CCR2)', 'CCR2 AMI vs. Control'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='HF vs. Control Latent Factors (CCR2)',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors (CCR2)', 'CCR2 ICM vs. Control'),\
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
\
## ----- COMBINE PLOTS ----- ##\
setwd("/ix/djishnu/Mary/MI_data/Lavine_2024/plots+figures")\
# Combine the plots and save\
combined_plot <- upset_plot_1 / upset_plot_2 | upset_plot_3 / upset_plot_4\
ggsave("AF_combined_upset_inclusive_color_all.pdf", plot = combined_plot, width = 30, height = 14, dpi = 300)\
## ------------------------- ##\
\
\
\
# ------- MACROPHAGES -------#\
# ---------------------------#\
\
## PLOT 1 ##\
# Define data frames with lists of DEGs\
setwd("/ix/djishnu/Mary/MI_data/Lavine_2024/plots+figures/Mac_upset_plot_lists/size_matched")\
Mac_AMI_vs_Control_df <- read.csv("./Mac_AMI_vs_Control_DEGs_SizeMatched.csv")\
Mac_ICM_vs_Control_df <- read.csv("./Mac_ICM_vs_Control_DEGs_SizeMatched.csv")\
Mac_LFs_df <- read.csv("./Mac_121423_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(Mac_AMI_vs_Control_df))\
print(colnames(Mac_ICM_vs_Control_df))\
print(colnames(Mac_LFs_df))\
\
# Define the features as 'Genes'\
Mac_AMI_vs_Control <- as.list(Mac_AMI_vs_Control_df[[1]])  # Assuming the genes are in the first column\
Mac_ICM_vs_Control <- as.list(Mac_ICM_vs_Control_df[[1]])\
Mac_LFs <- as.list((Mac_LFs_df[[2]]))\
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
gene_lists <- list(Mac_AMI_vs_Control, Mac_ICM_vs_Control, Mac_LFs)\
colnames <- c("AMI vs. Control (size-matched)","ICM vs. Control (size-matched)","HF vs. Control Latent Factors")\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_1 <- upset(\
  df, c("AMI vs. Control (size-matched)","ICM vs. Control (size-matched)","HF vs. Control Latent Factors"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18)) +\
      coord_cartesian(ylim = c(0, 10))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('AMI vs. Control (size-matched)', 'ICM vs. Control (size-matched)'),\
    c('AMI vs. Control (size-matched)', 'HF vs. Control Latent Factors'),\
    c('ICM vs. Control (size-matched)', 'HF vs. Control Latent Factors')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'AMI vs. Control (size-matched)'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='HF vs. Control Latent Factors',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'ICM vs. Control (size-matched)'),\
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
setwd("/ix/djishnu/Mary/MI_data/Lavine_2024/plots+figures/Mac_upset_plot_lists/not_size_matched")\
Mac_AMI_vs_Control_df <- read.csv("./Mac_AMI_vs_Control_DEGs.csv")\
Mac_ICM_vs_Control_df <- read.csv("./Mac_ICM_vs_Control_DEGs.csv")\
Mac_LFs_df <- read.csv("./Mac_121423_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(Mac_AMI_vs_Control_df))\
print(colnames(Mac_ICM_vs_Control_df))\
print(colnames(Mac_LFs_df))\
\
# Define the features as 'Genes'\
Mac_AMI_vs_Control <- as.list(Mac_AMI_vs_Control_df[[1]])  # Assuming the genes are in the first column\
Mac_ICM_vs_Control <- as.list(Mac_ICM_vs_Control_df[[1]])\
Mac_LFs <- as.list((Mac_LFs_df[[2]]))\
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
gene_lists <- list(Mac_AMI_vs_Control, Mac_ICM_vs_Control, Mac_LFs)\
colnames <- c("AMI vs. Control","ICM vs. Control","HF vs. Control Latent Factors")\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_2 <- upset(\
  df, c("AMI vs. Control","ICM vs. Control","HF vs. Control Latent Factors"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18)) +\
      coord_cartesian(ylim = c(0, 1000))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('AMI vs. Control', 'ICM vs. Control'),\
    c('AMI vs. Control', 'HF vs. Control Latent Factors'),\
    c('ICM vs. Control', 'HF vs. Control Latent Factors')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'AMI vs. Control'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='HF vs. Control Latent Factors',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'ICM vs. Control'),\
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
setwd("/ix/djishnu/Mary/MI_data/Lavine_2024/plots+figures/Mac_distance_upset_plot_lists/size_matched")\
Mac_AMI_vs_Control_df <- read.csv("./Mac_AMI_vs_Control_DEGs_SizeMatched.csv")\
Mac_ICM_vs_Control_df <- read.csv("./Mac_ICM_vs_Control_DEGs_SizeMatched.csv")\
Mac_LFs_df <- read.csv("./Mac_dist_121223_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(Mac_AMI_vs_Control_df))\
print(colnames(Mac_ICM_vs_Control_df))\
print(colnames(Mac_LFs_df))\
\
# Define the features as 'Genes'\
Mac_AMI_vs_Control <- as.list(Mac_AMI_vs_Control_df[[1]])  # Assuming the genes are in the first column\
Mac_ICM_vs_Control <- as.list(Mac_ICM_vs_Control_df[[1]])\
Mac_LFs <- as.list((Mac_LFs_df[[2]]))\
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
gene_lists <- list(Mac_AMI_vs_Control, Mac_ICM_vs_Control, Mac_LFs)\
colnames <- c("AMI vs. Control (size-matched)","ICM vs. Control (size-matched)","dist HF vs. Control Latent Factors")\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_3 <- upset(\
  df, c("AMI vs. Control (size-matched)","ICM vs. Control (size-matched)","dist HF vs. Control Latent Factors"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18)) +\
      coord_cartesian(ylim = c(0, 10))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('AMI vs. Control (size-matched)', 'ICM vs. Control (size-matched)'),\
    c('AMI vs. Control (size-matched)', 'dist HF vs. Control Latent Factors'),\
    c('ICM vs. Control (size-matched)', 'dist HF vs. Control Latent Factors')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('dist HF vs. Control Latent Factors', 'AMI vs. Control (size-matched)'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='dist HF vs. Control Latent Factors',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('dist HF vs. Control Latent Factors', 'ICM vs. Control (size-matched)'),\
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
setwd("/ix/djishnu/Mary/MI_data/Lavine_2024/plots+figures/Mac_distance_upset_plot_lists/not_size_matched")\
Mac_AMI_vs_Control_df <- read.csv("./Mac_AMI_vs_Control_DEGs.csv")\
Mac_ICM_vs_Control_df <- read.csv("./Mac_ICM_vs_Control_DEGs.csv")\
Mac_LFs_df <- read.csv("./Mac_dist_121223_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(Mac_AMI_vs_Control_df))\
print(colnames(Mac_ICM_vs_Control_df))\
print(colnames(Mac_LFs_df))\
\
# Define the features as 'Genes'\
Mac_AMI_vs_Control <- as.list(Mac_AMI_vs_Control_df[[1]])  # Assuming the genes are in the first column\
Mac_ICM_vs_Control <- as.list(Mac_ICM_vs_Control_df[[1]])\
Mac_LFs <- as.list((Mac_LFs_df[[2]]))\
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
gene_lists <- list(Mac_AMI_vs_Control, Mac_ICM_vs_Control, Mac_LFs)\
colnames <- c("AMI vs. Control","ICM vs. Control","dist HF vs. Control Latent Factors")\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_4 <- upset(\
  df, c("AMI vs. Control","ICM vs. Control","dist HF vs. Control Latent Factors"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18)) +\
      coord_cartesian(ylim = c(0, 1000))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('AMI vs. Control', 'ICM vs. Control'),\
    c('AMI vs. Control', 'dist HF vs. Control Latent Factors'),\
    c('ICM vs. Control', 'dist HF vs. Control Latent Factors')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('dist HF vs. Control Latent Factors', 'AMI vs. Control'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='dist HF vs. Control Latent Factors',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('dist HF vs. Control Latent Factors', 'ICM vs. Control'),\
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
\
\
## ----- COMBINE PLOTS ----- ##\
setwd("/ix/djishnu/Mary/MI_data/Lavine_2024/plots+figures")\
# Combine the plots and save\
combined_plot <- upset_plot_1 / upset_plot_2 | upset_plot_3 / upset_plot_4\
ggsave("Mac_combined_upset_inclusive_color_all.pdf", plot = combined_plot, width = 30, height = 14, dpi = 300)\
# ------------------------- ##\
\
\
\
## ------- SUPPLEMENTARY PLOTS ------- ##\
\
## PLOT 1 ##\
# Define data frames with lists of DEGs\
setwd("/ix/djishnu/Mary/MI_data/Lavine_2024/plots+figures/Mac_upset_plot_lists/size_matched")\
MacCCR2_AMI_vs_Control_df <- read.csv("./MacCCR2+_AMI_vs_Control_DEGs_SizeMatched.csv")\
MacCCR2_ICM_vs_Control_df <- read.csv("./MacCCR2+_ICM_vs_Control_DEGs_SizeMatched.csv")\
Mac_LFs_df <- read.csv("./Mac_121423_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(MacCCR2_AMI_vs_Control_df))\
print(colnames(MacCCR2_ICM_vs_Control_df))\
print(colnames(Mac_LFs_df))\
\
# Define the features as 'Genes'\
MacCCR2_AMI_vs_Control <- as.list(MacCCR2_AMI_vs_Control_df[[1]])\
MacCCR2_ICM_vs_Control <- as.list(MacCCR2_ICM_vs_Control_df[[1]])\
Mac_LFs <- as.list((Mac_LFs_df[[2]]))\
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
gene_lists <- list(MacCCR2_AMI_vs_Control, \
                   MacCCR2_ICM_vs_Control, Mac_LFs)\
colnames <- c("CCR2 AMI vs. Control (size-matched)",\
              "CCR2 ICM vs. Control (size-matched)","HF vs. Control Latent Factors")\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_1 <- upset(\
  df, c("CCR2 AMI vs. Control (size-matched)",\
        "CCR2 ICM vs. Control (size-matched)","HF vs. Control Latent Factors"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18)) +\
      coord_cartesian(ylim = c(0, 10))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('CCR2 AMI vs. Control (size-matched)', 'CCR2 ICM vs. Control (size-matched)'),\
    c('CCR2 AMI vs. Control (size-matched)', 'HF vs. Control Latent Factors'),\
    c('CCR2 ICM vs. Control (size-matched)', 'HF vs. Control Latent Factors')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'CCR2 AMI vs. Control (size-matched)'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='HF vs. Control Latent Factors',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'CCR2 ICM vs. Control (size-matched)'),\
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
setwd("/ix/djishnu/Mary/MI_data/Lavine_2024/plots+figures/Mac_upset_plot_lists/not_size_matched")\
MacCCR2_ICM_vs_Control_df <- read.csv("./MacCCR2+_ICM_vs_Control_DEGs.csv")\
MacCCR2_ICM_vs_Control_df <- read.csv("./MacCCR2+_ICM_vs_Control_DEGs.csv")\
Mac_LFs_df <- read.csv("./Mac_121423_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(MacCCR2_ICM_vs_Control_df))\
print(colnames(MacCCR2_ICM_vs_Control_df))\
print(colnames(Mac_LFs_df))\
\
# Define the features as 'Genes'\
MacCCR2_AMI_vs_Control <- as.list(MacCCR2_AMI_vs_Control_df[[1]])\
MacCCR2_ICM_vs_Control <- as.list(MacCCR2_ICM_vs_Control_df[[1]])\
Mac_LFs <- as.list((Mac_LFs_df[[2]]))\
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
gene_lists <- list(MacCCR2_AMI_vs_Control, \
                   MacCCR2_ICM_vs_Control, Mac_LFs)\
colnames <- c("CCR2 AMI vs. Control",\
              "CCR2 ICM vs. Control","HF vs. Control Latent Factors")\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_2 <- upset(\
  df, c("CCR2 AMI vs. Control",\
        "CCR2 ICM vs. Control","HF vs. Control Latent Factors"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18)) +\
      coord_cartesian(ylim = c(0, 10))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('CCR2 AMI vs. Control', 'CCR2 ICM vs. Control'),\
    c('CCR2 AMI vs. Control', 'HF vs. Control Latent Factors'),\
    c('CCR2 ICM vs. Control', 'HF vs. Control Latent Factors')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'CCR2 AMI vs. Control'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='HF vs. Control Latent Factors',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('HF vs. Control Latent Factors', 'CCR2 ICM vs. Control'),\
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
setwd("/ix/djishnu/Mary/MI_data/Lavine_2024/plots+figures/Mac_distance_upset_plot_lists/size_matched")\
MacCCR2_AMI_vs_Control_df <- read.csv("./MacCCR2+_AMI_vs_Control_DEGs_SizeMatched.csv")\
MacCCR2_ICM_vs_Control_df <- read.csv("./MacCCR2+_ICM_vs_Control_DEGs_SizeMatched.csv")\
Mac_LFs_df <- read.csv("./Mac_dist_121223_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(MacCCR2_AMI_vs_Control_df))\
print(colnames(MacCCR2_ICM_vs_Control_df))\
print(colnames(Mac_LFs_df))\
\
# Define the features as 'Genes'\
MacCCR2_AMI_vs_Control <- as.list(MacCCR2_AMI_vs_Control_df[[1]])\
MacCCR2_ICM_vs_Control <- as.list(MacCCR2_ICM_vs_Control_df[[1]])\
Mac_LFs <- as.list((Mac_LFs_df[[2]]))\
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
gene_lists <- list(MacCCR2_AMI_vs_Control,\
                   MacCCR2_ICM_vs_Control, Mac_LFs)\
colnames <- c("CCR2 AMI vs. Control (size-matched)",\
              "CCR2 ICM vs. Control (size-matched)","dist HF vs. Control Latent Factors")\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_3 <- upset(\
  df, c("CCR2 AMI vs. Control (size-matched)",\
        "CCR2 ICM vs. Control (size-matched)","dist HF vs. Control Latent Factors"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18)) +\
      coord_cartesian(ylim = c(0, 10))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('CCR2 AMI vs. Control (size-matched)', 'CCR2 ICM vs. Control (size-matched)'),\
    c('CCR2 AMI vs. Control (size-matched)', 'dist HF vs. Control Latent Factors'),\
    c('CCR2 ICM vs. Control (size-matched)', 'dist HF vs. Control Latent Factors')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('dist HF vs. Control Latent Factors', 'CCR2 AMI vs. Control (size-matched)'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='dist HF vs. Control Latent Factors',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('dist HF vs. Control Latent Factors', 'CCR2 ICM vs. Control (size-matched)'),\
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
setwd("/ix/djishnu/Mary/MI_data/Lavine_2024/plots+figures/Mac_distance_upset_plot_lists/not_size_matched")\
MacCCR2_AMI_vs_Control_df <- read.csv("./MacCCR2+_AMI_vs_Control_DEGs.csv")\
MacCCR2_ICM_vs_Control_df <- read.csv("./MacCCR2+_ICM_vs_Control_DEGs.csv")\
Mac_LFs_df <- read.csv("./Mac_dist_121223_LFs.csv")\
\
# Inspect column names to ensure correct referencing\
print(colnames(MacCCR2_AMI_vs_Control_df))\
print(colnames(MacCCR2_ICM_vs_Control_df))\
print(colnames(Mac_LFs_df))\
\
# Define the features as 'Genes'\
MacCCR2_AMI_vs_Control <- as.list(MacCCR2_AMI_vs_Control_df[[1]])\
MacCCR2_ICM_vs_Control <- as.list(MacCCR2_ICM_vs_Control_df[[1]])\
Mac_LFs <- as.list((Mac_LFs_df[[2]]))\
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
gene_lists <- list(MacCCR2_AMI_vs_Control, \
                   MacCCR2_ICM_vs_Control, Mac_LFs)\
colnames <- c("CCR2 AMI vs. Control",\
              "CCR2 ICM vs. Control","dist HF vs. Control Latent Factors")\
result <- combine_gene_lists(gene_lists, colnames)\
head(result)\
# Ensure column names and remove Genes column\
print(colnames(result))\
# Remove the 'Genes' column using base R\
df <- result[, -which(colnames(result) == "Genes")]\
\
upset_plot_4 <- upset(\
  df, c("CCR2 AMI vs. Control",\
        "CCR2 ICM vs. Control","dist HF vs. Control Latent Factors"),\
  set_sizes=(\
    upset_set_size()\
    + ylab('Set Size') + theme(axis.title.x = element_text(size = 18)) +\
      coord_cartesian(ylim = c(0, 10))\
  ),\
  width_ratio = 0.2,\
  themes = upset_default_themes(text = element_text(size = 20)),\
  mode = 'inclusive_intersection',  # Set the mode to inclusive intersection\
  intersections = list(\
    c('CCR2 AMI vs. Control', 'CCR2 ICM vs. Control'),\
    c('CCR2 AMI vs. Control', 'dist HF vs. Control Latent Factors'),\
    c('CCR2 ICM vs. Control', 'dist HF vs. Control Latent Factors')\
  ),\
  sort_intersections= 'descending', #FALSE\
  queries=list(\
    upset_query(\
      intersect = c('dist HF vs. Control Latent Factors', 'CCR2 AMI vs. Control'),\
      color='blue',\
      fill='blue',\
      only_components=c('intersections_matrix', 'Intersection size')\
    ),\
    upset_query(\
      set='dist HF vs. Control Latent Factors',\
      fill='blue'\
    ),\
    upset_query(\
      intersect = c('dist HF vs. Control Latent Factors', 'CCR2 ICM vs. Control'),\
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
\
\
## ----- COMBINE PLOTS ----- ##\
# Combine the plots and save\
combined_plot <- upset_plot_1 / upset_plot_2 | upset_plot_3 / upset_plot_4\
# ------------------------- ##}