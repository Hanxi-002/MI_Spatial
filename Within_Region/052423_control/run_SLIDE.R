
#devtools::install_github("Hanxi-002/SLIDEHelper", ref = "main")
library(doParallel)
library(SLIDE)
library(SLIDEHelper)

x_path <- "x.csv"
er_path <- "results/final_delta_0.01_lambda_0.5.rds"
out_path <- "SLIDE_Result/"

z_matrix <- CalcZMatrix(x_path, er_path, out_path)

y_path <- "y.csv"
SLIDE_res <- runSLIDE(y_path, z_path = NULL, z_matrix, er_path, do_interacts = TRUE, niter = 500, spec = 0.1)

num_top_feats <- 10
condition <- "auc"
Final_res <- SLIDEHelper::GetTopFeatures(x_path, y_path, er_path, out_path, SLIDE_res, num_top_feats = 10, condition)

saveRDS(Final_res, "SLIDE_Result/SLIDE_res.rds")

