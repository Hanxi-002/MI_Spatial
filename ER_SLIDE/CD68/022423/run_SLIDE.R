#devtools::install_github("Hanxi-002/SLIDEHelper", ref = "main")
setwd("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/CD68/022423")
library(doParallel)
library(ggplot2)
library(SLIDE)
library(SLIDEHelper)
library(ggraph)

getwd()
x_path <- "x.csv"
er_path <- "Results/ER_Results/final_er_output.rds"
out_path <- "Results/SLIDE_Results/"

z_matrix <- CalcZMatrix(x_path, er_path, out_path)

y_path <- "y.csv"
SLIDE_res <- runSLIDE(y_path, z_path = NULL, z_matrix, er_path, do_interacts = TRUE, niter = 500, spec = 0.1)

num_top_feats <- 10
condition <- "auc"
Final_res <- SLIDEHelper::GetTopFeatures(x_path, y_path, er_path, out_path, SLIDE_res, num_top_feats = 10, condition)

saveRDS(Final_res, "Results/SLIDE_Res/SLIDE_res.rds")

plotSigGenes(Final_res, out_path, plot_interaction = TRUE)
CalcControlPerformance(z_matrix = z_matrix, y_path, Final_res, niter = 1000, condition, out_path)
