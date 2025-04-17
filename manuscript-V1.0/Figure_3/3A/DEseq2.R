setwd("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllFibro/081224/DESeq")
source('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Find_DEGs_Helper.R')


count <- t(read.csv("../Data/x.csv", row.names = 1))
meta <- read.csv("../Data/y.csv", row.names = 1)
DEGs = calculate_DEGs(count, meta, file_name = 'DESeq2_res_20250319.RDS')


custom_color = c("grey30", "grey30", "grey30", "#8FBC8F")
file_name = "volcano_20250319.pdf"
plot_volcano(DEGs = DEGs, custom_color = custom_color, file_name = file_name)

