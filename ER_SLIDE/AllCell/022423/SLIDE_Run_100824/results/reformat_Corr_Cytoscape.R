library(dplyr)


convert_to_long_format <- function(cor_matrix) {
  features <- colnames(cor_matrix)

  cor_long <- as.data.frame(as.table(as.matrix(cor_matrix)))
  colnames(cor_long) <- c("target", "source", "correlation")
  

  cor_long <- cor_long %>% filter(target != source)
  cor_long <- cor_long %>%
    filter(as.numeric(factor(target)) < as.numeric(factor(source)))
  
  return(cor_long)
}


read_convert_and_save_txt <- function(folder_path, output_folder) {
  
  csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  
  for (file in csv_files) {
    cor_matrix <- read.csv(file, row.names = 1)
    cor_long <- convert_to_long_format(cor_matrix)
    file_name <- tools::file_path_sans_ext(basename(file))
    output_file <- file.path(output_folder, paste0(file_name, "_transformed.txt"))
    write.table(cor_long, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}


folder_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423/SLIDE_Run_100824/results/0.01_0.5_out/corr_networks/threh_0_5"        # Specify your input folder path here
output_folder <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423/SLIDE_Run_100824/results/0.01_0.5_out/corr_networks/cytoscape"    # Specify the output folder path here


if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}


read_convert_and_save_txt(folder_path, output_folder)


