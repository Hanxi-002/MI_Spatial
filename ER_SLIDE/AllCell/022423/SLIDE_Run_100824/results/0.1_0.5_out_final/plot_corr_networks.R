source('/ix/djishnu/Hanxi/Common_R/PlotCorrNetwork_Helper.R')

x <- as.matrix(read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423/x.csv", row.names = 1))
LF <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423/SLIDE_Run_100824/results/0.1_0.5_out_final/plotSigGenes_data.RDS")
threshold = 0.5 # set to whatever you want
#setwd("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423/SLIDE_Run_100824/results/0.01_0.5_out/corr_networks")
output_path = "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllCell/022423/SLIDE_Run_100824/results/0.1_0.5_out_final/corr_networks"
PlotCorNetwork(x, LF, threshold, layout = "spring", output_path)