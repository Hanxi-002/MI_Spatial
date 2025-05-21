source("/ix/djishnu/Hanxi/Common_R/PlotCorrNetwork_Helper.R")

x <- as.matrix(read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/Data/ccr2_af_x.csv", row.names = 1))
LF <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/plotSigGenes_data.RDS")
LF[94, ]$names <- "HLA.B"
threshold = 0.5

#setwd("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/corr_networks_v2")
PlotCorNetwork(x, LF, threshold, output_path = "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/corr_networks_v2")
