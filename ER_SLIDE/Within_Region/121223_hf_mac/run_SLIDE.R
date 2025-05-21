#devtools::install_github("Hanxi-002/SLIDEHelper", ref = "main")
setwd("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac")
library(doParallel)
library(ggplot2)
library(SLIDE)
library(SLIDEHelper)
library(ggraph)

x_path <- "Data/x.csv"
er_path <- "results/ER_Results/final_delta_0.001_lambda_0.5.rds"
out_path <- "results/SLIDE_Results/"

z_matrix <- CalcZMatrix(x_path, er_path, out_path)

y_path <- "Data/y.csv"
SLIDE_res <- runSLIDE(y_path, z_path = NULL, z_matrix, er_path, do_interacts = TRUE, niter = 500, spec = 0.1)

num_top_feats <- 10
condition <- "auc"
Final_res <- SLIDEHelper::GetTopFeatures(x_path, y_path, er_path, out_path, SLIDE_res, num_top_feats = 10, condition)

saveRDS(Final_res, "results/SLIDE_Results/SLIDE_res.rds")
CalcControlPerformance(z_matrix = z_matrix, y_path, Final_res, niter = 1000, condition, out_path)
plotSigGenes(Final_res, out_path, plot_interaction = TRUE)
