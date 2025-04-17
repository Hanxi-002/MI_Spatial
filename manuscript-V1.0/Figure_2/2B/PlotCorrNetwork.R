source("/ix/djishnu/Hanxi/Common_R/PlotCorrNetwork_Helper.R")

x <- as.matrix(read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/CD68/121423/Data/x.csv", row.names = 1))
LF <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/CD68/121423/results/SLIDE_Results/plotSigGenes_data.RDS")
threshold = 0.5
#setwd("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/CD68/121423/results/corr_network")

output_path = "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/CD68/121423/results/corr_network/thresh_0_5/circular"
layout = 'circular'
PlotCorNetwork(x, LF, threshold, layout = layout, output_path = output_path)


output_path = "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/CD68/121423/results/corr_network/thresh_0_5/spring"
layout = 'spring'
PlotCorNetwork(x, LF, threshold, layout = layout, output_path = output_path)