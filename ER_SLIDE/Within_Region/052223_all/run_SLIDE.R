
#devtools::install_github("Hanxi-002/SLIDEHelper", ref = "main")
library(doParallel)
library(SLIDE)
library(SLIDEHelper)

x_path <- "/ix/djishnu/Hanxi/MI_Spatial/Within_Region/052223_all/x.csv"
er_path <- "/ix/djishnu/Hanxi/MI_Spatial/Within_Region/052223_all/results/final_delta_0.01_lambda_0.5.rds"
out_path <- "/ix/djishnu/Hanxi/MI_Spatial/Within_Region/052223_all/SLIDE_results/"

z_matrix <- CalcZMatrix(x_path, er_path, out_path)

y_path <- "/ix/djishnu/Hanxi/MI_Spatial/Within_Region/052223_all/y.csv"
SLIDE_res <- runSLIDE(y_path, z_path = NULL, z_matrix, er_path, do_interacts = TRUE, niter = 500, spec = 0.1)

num_top_feats <- 10
condition <- "corr"
Final_res <- SLIDEHelper::GetTopFeatures(x_path, y_path, er_path, out_path, SLIDE_res, num_top_feats = 10, condition = "corr")

saveRDS(Final_res, "/ix/djishnu/Hanxi/MI_Spatial/Within_Region/052223_all/SLIDE_results/SLIDE_res.rds")
