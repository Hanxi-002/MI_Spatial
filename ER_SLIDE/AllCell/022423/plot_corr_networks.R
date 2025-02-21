source("/ix/djishnu/Hanxi/Common_R/PlotCorrNetwork_Helper.R")

x <- as.matrix(read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423/x.csv", row.names = 1))
LF <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423/results/SLIDE_Results_071223/plotSigGenes_data.RDS")
threshold = 0.3 # set to whatever you want
output_path = "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423/results/SLIDE_Results_071223/corr_networks"
PlotCorNetwork(x, LF, threshold, output_path = output_path)
