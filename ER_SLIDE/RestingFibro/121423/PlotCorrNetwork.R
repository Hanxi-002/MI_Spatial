source("/ix/djishnu/Hanxi/Common_R/PlotCorrNetwork.R")

x <- as.matrix(read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/RestingFibro/121423/Data/x.csv", row.names = 1))
LF <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/RestingFibro/121423/results/SLIDE_Results/plotSigGenes_data.RDS")
threshold = 0.3
setwd("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/RestingFibro/121423/results/corr_network")

PlotCorNetwork(x, LF, threshold)