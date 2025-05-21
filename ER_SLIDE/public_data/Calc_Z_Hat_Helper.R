library(Seurat)
library(ggplot2)
library(SLIDE)

####################################################################################################################################
# Observe the sign order of the latent factor in the original SLIDE run
find_LF_sign <- function(y_path, z_path, z_score_column, custom_labels = NULL){
  # plot the LF of interest's z score with the original y. 
  y <- read.csv(y_path, row.names = 1)
  SLIDE_z <- read.csv(z_path, row.names = 1)
  z_column_df <- as.data.frame(SLIDE_z[, z_score_column])
  names(z_column_df)[1] = 'score'
  row.names(z_column_df) <- row.names(SLIDE_z)
  if (sum(row.names(y) != row.names(z_column_df)) != 0 ) {warning("The row order of y and z is different. ")}
  z_column_df$condition = as.factor(y[ , 1])
  
  # Convert condition to factor first
  p <- ggplot(z_column_df, aes(x = condition, y = score)) +
    geom_boxplot() +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          legend.key.size = unit(1, "cm"))+
    labs(title = paste0("Z", z_score_column, " in GeoMx"),
         x = "Condition",
         y = "Score")
  
  if(!is.null(custom_labels))  {
    p <- p + scale_x_discrete(labels = custom_labels)
  }
  p
}


####################################################################################################################################
# find new Z

# Function to preprocess scRNA-seq data and ER results
preprocess_data <- function(x, er_results_path) {
  
  er_res <- readRDS(er_results_path)
  
  intersecting_cols <- intersect(colnames(x), row.names(er_res$A))
  x_intersected <- x[, intersecting_cols, drop = FALSE]
  
  missing_elements <- row.names(er_res$A)[!(row.names(er_res$A) %in% colnames(x))]
  missing_zeros <- matrix(0, 
                          nrow = nrow(x), 
                          ncol = length(missing_elements), 
                          dimnames = list(NULL, missing_elements))
  
  x_complete <- cbind(x_intersected, missing_zeros)
  
  return(list(
    x_complete = x_complete,
    er_res = er_res
  ))
}


extract_z_scores <- function(preprocessed_data, z_score_column) {
  # construct the new z for the new dataset
  z_hat <- predZ(preprocessed_data$x_complete, preprocessed_data$er_res)
  
  # only look at the LF we are interested in
  z_hat_99 <- as.data.frame(z_hat[, z_score_column])
  colnames(z_hat_99) <- "score"
  
  # subset the seurat object
  z_hat_99_seurat <- subset(x = preprocessed_data$seurat_obj, 
                            cells = rownames(z_hat_99))
  
  
  z_hat_99$condition <- z_hat_99_seurat@meta.data$major_labl[match(rownames(z_hat_99), 
                                                                   colnames(z_hat_99_seurat))]
  
  return(z_hat_99)
}



extract_z_scores_Lavine <- function(preprocessed_data, z_score_column) {
  # construct the new z for the new dataset
  z_hat <- predZ(preprocessed_data$x_complete, preprocessed_data$er_res)
  
  # only look at the LF we are interested in
  z_hat_99 <- as.data.frame(z_hat[, z_score_column])
  colnames(z_hat_99) <- "score"
  
  # subset the seurat object
  z_hat_99_seurat <- subset(x = preprocessed_data$seurat_obj, 
                            cells = rownames(z_hat_99))
  
  
  z_hat_99$condition <- z_hat_99_seurat$HF.etiology[match(rownames(z_hat_99), 
                                                          colnames(z_hat_99_seurat))]
  
  return(z_hat_99)
}


# Function to create boxplot of z-scores
plot_z_score_boxplot <- function(z_hat_99, title, custom_order = NULL) {
  
  if (!is.null(custom_order)) {
    z_hat_99$condition <- factor(z_hat_99$condition, levels = custom_order)
  }
  
  ggplot(z_hat_99, aes(x = condition, y = score)) +
    geom_boxplot() + 
    labs(
      title = title,
      x = "Condition",
      y = "Z_Score"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          legend.key.size = unit(1, "cm"))
}

# Main script execution
main <- function(seurat_path, er_results_path, z_score_column, plot_title, custom_order) {
  
  seurat_obj <- readRDS(seurat_path)
  x <- as.matrix(seurat_obj@assays$SCT@data)
  x <- t(x)
  
  preprocessed_data <- preprocess_data(x, er_results_path)
  
  preprocessed_data$seurat_obj <- seurat_obj
  z_hat_99 <- extract_z_scores(preprocessed_data, z_score_column)
  
  if (length(custom_order) != length(unique(z_hat_99$condition))){
    cat('subsetting z scores to only keep conditions specified in custom_order...')
    z_hat_99 <- z_hat_99[z_hat_99$condition %in% custom_order, ]
  }
  
  plot <- plot_z_score_boxplot(z_hat_99, title = plot_title, custom_order = custom_order)
  print(plot)
  
  return(list(
    preprocessed_data = preprocessed_data,
    z_hat_99 = z_hat_99,
    plot = plot
  ))
}


main_Lavine <- function(seurat_path, er_results_path, z_score_column, plot_title, custom_order) {
  
  seurat_obj <- readRDS(seurat_path)
  x <- as.matrix(seurat_obj@assays$RNA$data)
  x <- t(x)
  
  preprocessed_data <- preprocess_data(x, er_results_path)
  
  preprocessed_data$seurat_obj <- seurat_obj
  z_hat_99 <- extract_z_scores_Lavine(preprocessed_data, z_score_column)
  
  plot <- plot_z_score_boxplot(z_hat_99, title = plot_title, custom_order = custom_order)
  print(plot)
  
  return(list(
    preprocessed_data = preprocessed_data,
    z_hat_99 = z_hat_99,
    plot = plot
  ))
}


####################################################################################################################################
# Mann Whitney
perform_mw_tests <- function(df, condition_col = "condition", score_col = "score") {
  
  required_cols <- c(condition_col, score_col)
  if (!all(required_cols %in% colnames(df))) {
    stop("Missing required columns in dataframe")
  }
  
  # Get unique conditions
  conditions <- unique(df[[condition_col]])
  n_conditions <- length(conditions)
  
  # Create empty vectors to store results
  comparisons <- c()
  p_values <- c()
  
  if(n_conditions == 2) {
    formula_str <- paste(score_col, "~", condition_col)
    result <- wilcox.test(as.formula(formula_str), data = df)
    
    comparisons <- paste(conditions[1], "vs", conditions[2])
    p_values <- result$p.value
    
  } else if(n_conditions > 2) {
    # Perform specific pairwise comparisons
    for(i in 1:(n_conditions-1)) {
      for(j in (i+1):n_conditions) {
        # Subset data for the two conditions being compared
        test_data <- df[df[[condition_col]] %in% c(conditions[i], conditions[j]), ]
        
        # Mann-Whitney U test
        formula_str <- paste(score_col, "~", condition_col)
        result <- wilcox.test(as.formula(formula_str), data = test_data)
        
        # Store results
        comparisons <- c(comparisons, paste(conditions[i], "vs", conditions[j]))
        p_values <- c(p_values, result$p.value)
      }
    }
  } else {
    stop("Error: Need at least 2 conditions to perform the test")
  }
  
  # Create results dataframe
  results_df <- data.frame(
    p_val = p_values,
    row.names = comparisons
  )
  
  return(results_df)
}


####################################################################################################################################
# Cliff's Delta

# Function to perform pairwise Cliff's Delta calculations and save results
perform_cliffs_delta <- function(df) {
  # Load required package
  library(effsize)
  
  # Get unique conditions
  conditions <- unique(df$condition)
  n_conditions <- length(conditions)
  
  # Create empty vectors to store results
  comparisons <- c()
  delta_values <- c()
  
  if(n_conditions == 2) {
    # For exactly 2 conditions
    group1 <- df$score[df$condition == conditions[1]]
    group2 <- df$score[df$condition == conditions[2]]
    
    # Calculate Cliff's Delta
    result <- cliff.delta(group1, group2)
    
    # Store results
    comparisons <- paste(conditions[1], "vs", conditions[2])
    delta_values <- result$estimate
    
  } else if(n_conditions > 2) {
    # Perform specific pairwise comparisons
    for(i in 1:(n_conditions-1)) {
      for(j in (i+1):n_conditions) {
        # Get scores for each condition
        group1 <- df$score[df$condition == conditions[i]]
        group2 <- df$score[df$condition == conditions[j]]
        
        # Calculate Cliff's Delta
        result <- cliff.delta(group1, group2)
        
        # Store results
        comparisons <- c(comparisons, paste(conditions[i], "vs", conditions[j]))
        delta_values <- c(delta_values, result$estimate)
      }
    }
  } else {
    stop("Error: Need at least 2 conditions to calculate Cliff's Delta")
  }
  
  # Create results dataframe
  results_df <- data.frame(
    cliffs_delta = delta_values,
    row.names = comparisons
  )
  
  return(results_df)
}


####################################################################################################################################
# Q3 Box Plot

# baseline , string, is a condition in the condition column that used as the baseline for the ratio. 
counts_above_q3 <- function(df, baseline = 'CTRL') {
  
  # the input df has to have a score column and a condition column. 
  # can just use the output of the main function above. 
  q3_value <- quantile(df[df$condition == baseline, ]$score, 0.75)
  
  # Count total points for each condition
  total_counts <- aggregate(score ~ condition, data = df, FUN = length)
  names(total_counts)[2] <- "total_count"
  
  # Count points above Q3 for each condition
  counts_above <- aggregate(score ~ condition, 
                            data = df[df$score > q3_value,], 
                            FUN = length)
  names(counts_above)[2] <- "count_above_q3"
  
  # Manual merge using base R
  result <- total_counts
  result$count_above_q3 <- 0  # Initialize with zeros
  
  # Fill in the counts for conditions that have points above q3
  for(i in 1:nrow(counts_above)) {
    cond <- counts_above$condition[i]
    idx <- which(result$condition == cond)
    if(length(idx) > 0) {
      result$count_above_q3[idx] <- counts_above$count_above_q3[i]
    }
  }
  
  # Calculate the ratio
  result$ratio_above_q3 <- result$count_above_q3 / result$total_count
  
  return(result)
}


plot_ratios <- function(counts_df, custom_order = NULL) {
  # If custom order is provided, convert condition to factor with that order
  if (!is.null(custom_order)) {
    # Verify all conditions are in the custom order
    if (!all(counts_df$condition %in% custom_order)) {
      stop("All conditions must be present in custom_order")
    }
    if (!all(custom_order %in% counts_df$condition)) {
      warning("Some values in custom_order are not present in the data")
    }
    
    counts_df$condition <- factor(counts_df$condition, 
                                  levels = custom_order)
  }
  
  # Check if CI columns exist
  has_ci <- all(c("ci_lower", "ci_upper") %in% colnames(counts_df))
  
  # Create the base plot
  p <- ggplot(counts_df, aes(x = condition, y = ratio_above_q3)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "navy", size = 0.7, alpha = 0.7) +
    theme_minimal() +
    labs(title = "Ratio of Points Above Q3 by Condition",
         x = "Condition",
         y = "Ratio of Points > Q3") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          axis.line = element_line(color = "black"))
  
  # Add confidence intervals if they exist
  if (has_ci) {
    p <- p + 
      geom_errorbar(
        aes(ymin = ci_lower, ymax = ci_upper),
        width = 0.2,
        color = "black",
        size = 1
      ) +
      labs(title = "Ratio of Points Above Q3 by Condition",
           subtitle = "With 95% Clopper-Pearson confidence intervals",
           x = "Condition",
           y = "Ratio of Points > Q3")
  }
  
  return(p)
}

####################################################################################################################################
# Q3 Box Plot, calculate exact binomial test and CI for each group

calc_Clopper_Pearson_CIs <- function(df, count_col = "count_above_q3", total_col = "total_count", conf_level = 0.95) {
  # df is the output of counts_above_q3 ratio
  
  required_cols <- c(count_col, total_col)
  if (!all(required_cols %in% colnames(df))) {
    stop("Missing required columns in dataframe")
  }
  
  
  df$ci_lower <- NA
  df$ci_upper <- NA
  
  # Calculate CI for each condition
  for (i in 1:nrow(df)) {
    successes <- df[[count_col]][i]
    trials <- df[[total_col]][i]
    
    ci <- binom.test(
      x = successes,
      n = trials,
      conf.level = conf_level
    )$conf.int
    
    df$ci_lower[i] <- ci[1]
    df$ci_upper[i] <- ci[2]
  }
  
  return(df)
}


calc_SEM_CIs <- function(df, ratio_col = "ratio_above_q3", count_col = "count_above_q3", total_col = "total_count", condition_col = "condition", conf_level = 0.95) {
  # Check if required columns exist
  required_cols <- c(count_col, total_col)
  if (!all(required_cols %in% colnames(df))) {
    stop("Missing required columns in dataframe")
  }
  
  
  #df$ci_lower <- NA
  #df$ci_upper <- NA
  
  # Get z-score for the specified confidence level
  #z_score <- qnorm(1 - (1 - conf_level) / 2)
  
  
  df$sem <- sqrt(df[[ratio_col]] * (1 - df[[ratio_col]]) / df[[total_col]])
  
  df$ci_lower <- df[[ratio_col]]  - df$sem
  df$ci_upper <- df[[ratio_col]]  + df$sem
  
  # Ensure confidence intervals don't go below 0 or above 1
  df$ci_lower <- pmax(0, df$ci_lower)
  df$ci_upper <- pmin(1, df$ci_upper)
  
  # Return both the plot and the data with CIs
  return(df)
}


exact_binomial_test <- function(df, count_col = "count_above_q3", total_col = "total_count", condition_col = "condition", baseline = "IZ") {
  # df is the output of counts_above_q3
  # counts: number of counts of sucess
  # ns: number of trails
  # contidiont: different conditions
  # baseline: the baseline to comapre to
  # find the index of the baseline condition
  required_cols <- c(count_col, total_col, condition_col)
  if (!all(required_cols %in% colnames(df))) {
    stop("Missing required columns in dataframe")
  }
  
  conditions = df[, condition_col]
  counts = df[, count_col]
  ns = df[, total_col]
  
  baseline_idx <- which(conditions == baseline)
  if(length(baseline_idx) == 0) {
    stop(paste("Baseline condition '", baseline, "' not found in condition names."))
  }
  
  # get all other condition indices
  test_indices <- setdiff(1:length(counts), baseline_idx)
  
  results <- data.frame(
    condition = character(),
    condition_prop = numeric(),
    baseline_prop = numeric(),
    diff = numeric(),
    p_value = numeric(),
    test_type = character(),
    stringsAsFactors = FALSE
  )
  
  # get baseline values
  p_baseline <- counts[baseline_idx] / ns[baseline_idx]
  
  for (idx in test_indices) {
    successes <- counts[idx]
    trials <- ns[idx]
    p_test <- successes / trials
    
    # 2 tailed p value
    p_val <- binom.test(
      x = successes,
      n = trials,
      p = p_baseline,
      alternative = "two.sided"
    )$p.value
    test_type <- "two.sided"
    
    results <- rbind(results, data.frame(
      condition = conditions[idx],
      condition_prop = p_test,
      baseline_prop = p_baseline,
      diff = p_test - p_baseline,
      p_value = p_val,
      test_type = test_type
    ))
  }
  
  # Add significance stars based on p-values
  results$significance <- ""
  results$significance[results$p_value < 0.05] <- "*"
  results$significance[results$p_value < 0.01] <- "**"
  results$significance[results$p_value < 0.001] <- "***"
  
  return(results)
}



####################################################################################################################################
# Q3 Box Plot, calculate paired proportional z test

control_prop_test <- function(df, props_col = 'ratio_above_q3', total_col = 'total_count', condition_col = 'condition', baseline = "Control") {
  
  required_cols <- c(props_col, total_col, condition_col)
  if (!all(required_cols %in% colnames(df))) {
    stop("Missing required columns in dataframe")
  }
  
  props = df[, props_col]
  conditions = df[, condition_col]
  ns = df[, total_col]
  
  baseline_idx <- which(conditions == baseline)
  if(length(baseline_idx) == 0) {
    stop(paste("Baseline condition '", baseline, "' not found in condition names."))
  }
  
  # get all other condition indices
  test_indices <- setdiff(1:length(props), baseline_idx)
  
  # initiate an empty df to hold results
  results <- data.frame(
    condition = character(),
    condition_prop = numeric(),
    baseline_prop = numeric(),
    diff = numeric(),
    z_statistic = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # get baseline values
  p_baseline <- props[baseline_idx]
  n_baseline <- ns[baseline_idx]
  
  # calcualte the test statiscitcs
  for (idx in test_indices) {
    p_test <- props[idx]
    n_test <- ns[idx]
    
    
    p_pooled <- (p_baseline * n_baseline + p_test * n_test) / (n_baseline + n_test)
    
    
    se <- sqrt(p_pooled * (1 - p_pooled) * (1/n_baseline + 1/n_test))
    
    
    z <- (p_test - p_baseline) / se
    
    
    p_val <- 2 * pnorm(-abs(z))
    
    results <- rbind(results, data.frame(
      condition = conditions[idx],
      condition_prop = p_test,
      baseline_prop = p_baseline,
      diff = p_test - p_baseline,
      z_statistic = z,
      p_value = p_val
    ))
  }
  
  # add significance stars based on raw p-values
  results$significance <- ""
  results$significance[results$p_value < 0.05] <- "*"
  results$significance[results$p_value < 0.01] <- "**"
  results$significance[results$p_value < 0.001] <- "***"
  
  return(results)
}

