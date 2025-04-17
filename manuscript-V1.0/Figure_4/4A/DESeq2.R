setwd("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/DESeq")
source('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Find_DEGs_Helper.R')

count <- t(read.csv("../Data/x.csv", row.names = 1))
meta <- read.csv("../Data/y.csv", row.names = 1)
DEGs = calculate_DEGs(count, meta, file_name = 'DESeq2_res_20250319.RDS')


custom_color = c("grey30", "grey30", "grey30", "#FC4E2A")
file_name = "volcano_20250324.pdf"
plot_volcano(DEGs = DEGs, custom_color = custom_color, file_name = file_name)


