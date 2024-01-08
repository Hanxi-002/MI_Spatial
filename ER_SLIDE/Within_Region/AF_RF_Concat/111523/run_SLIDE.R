
#devtools::install_github("Hanxi-002/SLIDEHelper", ref = "main")
library(doParallel)
library(SLIDE)
library(SLIDEHelper)

x_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/AF_RF_Concat/111523/Data/x.csv"
er_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/AF_RF_Concat/111523/results/final_er_output.rds"
out_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/AF_RF_Concat/111523/SLIDE_results/"

z_matrix <- CalcZMatrix(x_path, er_path, out_path)

y_path <-  "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/AF_RF_Concat/111523/Data/y.csv"
SLIDE_res <- runSLIDE(y_path, z_path = NULL, z_matrix, er_path, do_interacts = TRUE, niter = 1000, spec = 0.1)

num_top_feats <- 10
condition <- "corr"
Final_res <- SLIDEHelper::GetTopFeatures(x_path, y_path, er_path, out_path, SLIDE_res, num_top_feats = 10, condition = "corr")
plotSigGenes(Final_res, out_path, plot_interaction = TRUE)

saveRDS(Final_res, paste0(out_path, "SLIDE_res.rds"))
