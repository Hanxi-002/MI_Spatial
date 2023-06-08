x_path <- "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/CD68/050123/Data/x.csv"
x <- as.matrix(utils::read.csv(x_path, row.names = 1))
x <- scale(x, T, T)
er_res <- readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/CD68/050123/ER_Results/final_delta_0.01_lambda_0.5.rds")

z_matrix <- predZ(x, er_res)
colnames(z_matrix) <- paste0("Z", c(1:ncol(z_matrix)))
write.csv(z_matrix, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/CD68/050123/ER_Results/z_matrix.csv", row.names = TRUE, col.names = TRUE)
