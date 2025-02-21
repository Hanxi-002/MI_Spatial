library(dplyr)

# Function to convert a correlation matrix into long format
convert_to_long_format <- function(cor_matrix) {
  # Get the row and column names (assume square matrix)
  features <- colnames(cor_matrix)
  
  # Create a long-format data frame with target, source, and correlation
  cor_long <- as.data.frame(as.table(as.matrix(cor_matrix)))
  colnames(cor_long) <- c("target", "source", "correlation")
  
  # Remove self-correlations (diagonal elements)
  cor_long <- cor_long %>% filter(target != source)
  
  # Ensure no repeated pairs (keep only the lower triangle of the matrix)
  cor_long <- cor_long %>%
    filter(as.numeric(factor(target)) < as.numeric(factor(source)))
  
  return(cor_long)
}

# Function to read, convert, and save each network as a txt file
read_convert_and_save_txt <- function(folder_path, output_folder) {
  # Get a list of all CSV files in the folder
  csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Loop through each CSV file, convert, and save as a new txt file
  for (file in csv_files) {
    # Read the correlation matrix from the CSV
    cor_matrix <- read.csv(file, row.names = 1)
    
    # Convert the correlation matrix to long format
    cor_long <- convert_to_long_format(cor_matrix)
    
    # Get the file name without extension for naming the new txt file
    file_name <- tools::file_path_sans_ext(basename(file))
    
    # Define the output file path for the txt file
    output_file <- file.path(output_folder, paste0(file_name, "_transformed.txt"))
    
    # Save the long-format data frame as a tab-delimited txt file
    write.table(cor_long, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}


folder_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllFibro/081224/results/SLIDE_results/0.1_0.5_out_final/corr_network"        # Specify your input folder path here
output_folder <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllFibro/081224/results/SLIDE_results/0.1_0.5_out_final/corr_network"    # Specify the output folder path here


if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

read_convert_and_save_txt(folder_path, output_folder)

