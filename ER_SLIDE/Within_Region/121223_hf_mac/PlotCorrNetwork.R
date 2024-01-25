source("/ix/djishnu/Hanxi/Common_R/PlotCorrNetwork.R")

x <- as.matrix(read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/Data/x.csv", row.names = 1))
LF <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/results/SLIDE_Results/plotSigGenes_data.RDS")
threshold = 0.9
setwd("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/results/corr_network")

PlotCorNetwork(x, LF, threshold)
