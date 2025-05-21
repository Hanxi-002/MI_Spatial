#devtools::install_github("Hanxi-002/SLIDEHelper", ref = "main")
setwd("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423")
library(doParallel)
library(ggplot2)
library(SLIDE)
library(SLIDEHelper)
library(ggraph)

x_path <- "x.csv"
er_path <- "results/final_delta_0.01_lambda_0.5.rds"
out_path <- "SLIDE_Results/"

z_matrix <- CalcZMatrix(x_path, er_path, out_path)

y_path <- "y.csv"
SLIDE_res <- runSLIDE(y_path, z_path = NULL, z_matrix, er_path, do_interacts = TRUE, niter = 300, spec = 0.2)

num_top_feats <- 10
condition <- "auc"
Final_res <- SLIDEHelper::GetTopFeatures(x_path, y_path, er_path, out_path, SLIDE_res, num_top_feats = 10, condition)

saveRDS(Final_res, "SLIDE_Results/SLIDE_res.rds")

plotSigGenes(Final_res, out_path, plot_interaction = TRUE)
