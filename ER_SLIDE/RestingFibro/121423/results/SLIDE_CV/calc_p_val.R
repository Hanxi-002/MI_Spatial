library(ggplot2)

parent_folder_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/RestingFibro/121423/results/SLIDE_CV"

# List all RDS files in the parent folder and its subdirectories
rds_files <- list.files(path = parent_folder_path, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)

# Read each RDS file and store the dataframes in a list
dataframes <- lapply(rds_files, function(file) {
  data <- readRDS(file)
  if ("final_corr" %in% names(data)) {
    return(data$final_corr)
  } else {
    warning(paste("final_corr not found in file:", file))
    return(NULL)
  }
})

big_dataframe <- do.call(rbind, dataframes)

lambda_boxplot = ggpubr::ggboxplot(data = big_dataframe, x = "method", y = "auc", palette = "aaas",
                                   fill = "method" ) + ggpubr::stat_compare_means(label = "p.signif")

lambda_boxplot

p_val <- ggpubr::compare_means(auc ~ method, data = big_dataframe)
