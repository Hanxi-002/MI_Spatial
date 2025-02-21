library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

######################## load cell count information for CellChat ########################
Geomx_V3 <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/Geomx_V3.RDS")
norm <- as.matrix(Geomx_V3@assayData$q_norm)
count_HF <- norm[, colnames(norm) %in% row.names(Geomx_V3@protocolData@data[Geomx_V3@protocolData@data$Status == 'HF', ])]
dim(count_HF)
count_Control <- norm[, colnames(norm) %in% row.names(Geomx_V3@protocolData@data[Geomx_V3@protocolData@data$Status == 'Control', ])]
dim(count_Control)

######################## load cell label information for CellChat ########################
meta <- Geomx_V3@phenoData@data
meta$new_column <- apply(meta, 1, function(row) {
  if (row["CD 68"]) {
    "CD68"
  } else if (row["resting fibroblast"]) {
    "RF"
  } else if (row["active fibroblast"]) {
    "AF"
  } else {
    NA
  }
})

# read in the CD68 distance labels
proximity <- read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/Data/y.csv", row.names = 1)
proximity_far <- proximity[proximity$y.row.names.y...in..row.names.CD68_HF.... == 0, , drop = FALSE]
proximity_close <- proximity[proximity$y.row.names.y...in..row.names.CD68_HF.... == 1, , drop = FALSE]

meta$new_column[row.names(meta) %in% row.names(proximity_far)] <- 'CD68_far'
meta$new_column[row.names(meta) %in% row.names(proximity_close)] <- 'CD68_close'

# now split meta into HF and Control, this is because we will do CCC on these 2 separately
meta_HF <- meta[row.names(meta) %in% row.names(Geomx_V3@protocolData@data[Geomx_V3@protocolData@data$Status == 'HF', ]), ]
meta_Control <- meta[row.names(meta) %in% row.names(Geomx_V3@protocolData@data[Geomx_V3@protocolData@data$Status == 'Control', ]), ]

# these are regions that are missing a AOI, either RF or AF. They were excluded from the SLIDE analysis. 
meta_HF[row.names(meta_HF) == 'DSP-1001660011816-D-A11.dcc', ]$new_column <-'CD68_far'
meta_HF[row.names(meta_HF) == 'DSP-1001660012268-C-G05.dcc', ]$new_column <-'CD68_close'
meta_HF[row.names(meta_HF) == 'DSP-1001660012268-C-D03.dcc', ]$new_column <-'CD68_far'
meta_HF[row.names(meta_HF) == 'DSP-1001660012268-C-D05.dcc', ]$new_column <-'CD68_far'
meta_HF[row.names(meta_HF) == 'DSP-1001660011816-D-E04.dcc', ]$new_column <-'CD68_far'
meta_HF[row.names(meta_HF) == 'DSP-1001660011816-D-E06.dcc', ]$new_column <-'CD68_far'
meta_HF[row.names(meta_HF) == 'DSP-1001660011816-D-E08.dcc', ]$new_column <-'CD68_far'

# we borrow the above code to output the newly formated meta data for multinichenet
# meta[row.names(meta) == 'DSP-1001660011816-D-A11.dcc', ]$new_column <-'CD68_far'
# meta[row.names(meta) == 'DSP-1001660012268-C-G05.dcc', ]$new_column <-'CD68_close'
# meta[row.names(meta) == 'DSP-1001660012268-C-D03.dcc', ]$new_column <-'CD68_far'
# meta[row.names(meta) == 'DSP-1001660012268-C-D05.dcc', ]$new_column <-'CD68_far'
# meta[row.names(meta) == 'DSP-1001660011816-D-E04.dcc', ]$new_column <-'CD68_far'
# meta[row.names(meta) == 'DSP-1001660011816-D-E06.dcc', ]$new_column <-'CD68_far'
# meta[row.names(meta) == 'DSP-1001660011816-D-E08.dcc', ]$new_column <-'CD68_far'
# 
# meta['Status'] <- Geomx_V3@protocolData@data$Status
# meta <- as.matrix(meta)
# write.csv(meta, "/ix/djishnu/Hanxi/MI_Spatial/CCC/MultiNicheNet/CCC_meta.csv", row.names = T)
########################### run cell chat for HF######################################
# this line does FDR control for the count matrix
#library(SLIDE)
# Function to analyze covariance structure and identify low-variance features
analyze_covariance <- function(data_matrix, 
                               variance_threshold = 0.01,
                               correlation_threshold = 0.9,
                               plot = TRUE) {
  # Convert input to matrix if it's not already
  data_matrix <- as.matrix(data_matrix)
  
  # Calculate variance for each feature
  feature_variance <- apply(data_matrix, 2, var)
  
  # Calculate covariance and correlation matrices
  cov_matrix <- cov(data_matrix)
  cor_matrix <- cor(data_matrix)
  
  # Identify low variance features
  low_var_features <- names(feature_variance[feature_variance < variance_threshold])
  
  # Find highly correlated feature pairs
  high_cor_pairs <- which(abs(cor_matrix) > correlation_threshold & 
                            abs(cor_matrix) < 1, arr.ind = TRUE)
  high_cor_pairs <- high_cor_pairs[high_cor_pairs[,1] < high_cor_pairs[,2], ]
  
  # Create summary statistics
  summary_stats <- list(
    feature_variance = sort(feature_variance),
    mean_variance = mean(feature_variance),
    median_variance = median(feature_variance),
    low_variance_features = low_var_features,
    high_correlation_pairs = high_cor_pairs,
    correlation_matrix = cor_matrix,
    covariance_matrix = cov_matrix
  )
  
  # Generate visualizations if requested
  if (plot) {
    # Variance distribution plot
    par(mfrow = c(2, 2))
    
    # Histogram of variances
    hist(feature_variance, 
         main = "Distribution of Feature Variances",
         xlab = "Variance",
         breaks = 30)
    abline(v = variance_threshold, col = "red", lty = 2)
    
    # Boxplot of variances
    boxplot(feature_variance,
            main = "Feature Variance Boxplot",
            ylab = "Variance")
    
    # Heatmap of correlation matrix
    image(cor_matrix,
          main = "Correlation Matrix Heatmap",
          xaxt = "n", yaxt = "n")
    
    # Reset plotting parameters
    par(mfrow = c(1, 1))
  }
  
  return(summary_stats)
}

# Function to remove low variance features
remove_low_variance <- function(data_matrix, 
                                variance_threshold = 0.01,
                                return_features = FALSE) {
  # Calculate variances
  feature_variance <- apply(data_matrix, 2, var)
  
  # Identify features to keep
  keep_features <- feature_variance >= variance_threshold
  
  # Create filtered dataset
  filtered_data <- data_matrix[, keep_features, drop = FALSE]
  
  if (return_features) {
    return(list(
      filtered_data = filtered_data,
      removed_features = names(feature_variance[!keep_features]),
      kept_features = names(feature_variance[keep_features]),
      feature_variances = feature_variance
    ))
  } else {
    return(filtered_data)
  }
}

# Function to identify and remove highly correlated features
remove_correlated_features <- function(data_matrix, 
                                       correlation_threshold = 0.9,
                                       method = "variance") {
  # Calculate correlation matrix
  cor_matrix <- cor(data_matrix)
  
  # Find highly correlated pairs
  high_cor <- which(abs(cor_matrix) > correlation_threshold & 
                      abs(cor_matrix) < 1, arr.ind = TRUE)
  
  # If no highly correlated features found, return original data
  if (nrow(high_cor) == 0) {
    return(data_matrix)
  }
  
  # Calculate feature importance (variance by default)
  if (method == "variance") {
    importance <- apply(data_matrix, 2, var)
  }
  
  # Initialize features to remove
  features_to_remove <- c()
  
  # For each correlated pair, remove the one with lower importance
  for (i in 1:nrow(high_cor)) {
    pair <- high_cor[i, ]
    if (importance[pair[1]] < importance[pair[2]]) {
      features_to_remove <- c(features_to_remove, colnames(data_matrix)[pair[1]])
    } else {
      features_to_remove <- c(features_to_remove, colnames(data_matrix)[pair[2]])
    }
  }
  
  # Remove duplicates
  features_to_remove <- unique(features_to_remove)
  
  # Return filtered dataset
  return(data_matrix[, !colnames(data_matrix) %in% features_to_remove])
}

analysis_results <- analyze_covariance(t(count_HF), 
                                       variance_threshold = 0.02,
                                       correlation_threshold = 0.9)

# Remove low variance features
filtered_data <- remove_low_variance(t(count_HF), 
                                     variance_threshold = 0.02,
                                     return_features = TRUE)

# Further remove highly correlated features
final_data <- remove_correlated_features(filtered_data$filtered_data,
                                         correlation_threshold = 0.8)
dim(final_data)
dim(count_HF)
count_HF <- t(final_data)

cellchat_HF <- createCellChat(object = count_HF, meta = meta_HF, group.by = 'new_column')
cellchat_HF <- addMeta(cellchat_HF, meta = meta_HF)
cellchat_HF <- setIdent(cellchat_HF, ident.use = "new_column") # set "labels" as default cell identity
levels(cellchat_HF@idents)
groupSize <- as.numeric(table(cellchat_HF@idents))

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat_HF@DB <- CellChatDB.use

cellchat_HF <- subsetData(cellchat_HF)
cellchat_HF <- identifyOverExpressedGenes(cellchat_HF)
cellchat_HF <- identifyOverExpressedInteractions(cellchat_HF)


cellchat_HF <- computeCommunProb(cellchat_HF)
cellchat_HF <- computeCommunProbPathway(cellchat_HF)
cellchat_HF <- aggregateNet(cellchat_HF)

groupSize <- as.numeric(table(cellchat_HF@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_HF@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F)
title(main = "Number of interactions", line = 0.1)
netVisual_circle(cellchat_HF@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F)
title(main = "Interaction weights/strength", line = 0.1)


mat <- cellchat_HF@net$weight
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat))
  title(main = rownames(mat)[i], line = 2)
}


#df.net <- subsetCommunication(cellchat_HF)
#write.csv(df.net, '/ix/djishnu/Hanxi/MI_Spatial/CCC/CellChat/HF_LR_Pairs.csv')
########################### run cell chat for Control######################################
cellchat <- createCellChat(object = count_Control, meta = meta_Control, group.by = 'new_column')

cellchat <- addMeta(cellchat, meta = meta_Control)
cellchat <- setIdent(cellchat, ident.use = "new_column") # set "labels" as default cell identity
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F)
title(main = "Number of interactions", line = 0.1)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F)
title(main = "Interaction weights/strength", line = 0.1)


mat <- cellchat@net$weight
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat))
  title(main = rownames(mat)[i], line = 2)
}

df.net <- subsetCommunication(cellchat)
write.csv(df.net, '/ix/djishnu/Hanxi/MI_Spatial/CCC/CellChat/Controls_LR_Pairs.csv')