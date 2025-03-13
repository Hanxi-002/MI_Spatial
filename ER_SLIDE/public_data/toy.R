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