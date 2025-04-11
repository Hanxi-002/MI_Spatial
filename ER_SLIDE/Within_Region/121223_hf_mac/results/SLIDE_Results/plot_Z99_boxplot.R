library(ggplot2)

z_mat <- read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/results/SLIDE_Results/z_matrix.csv", row.names = 1)
y <- read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/Data/y.csv", row.names = 1)

y$y.row.names.y...in..row.names.CD68_HF.... <- ifelse(y$y.row.names.y...in..row.names.CD68_HF....==0, "Control", "HF")

plot_data <- data.frame(group = as.factor(y$y.row.names.y...in..row.names.CD68_HF....),
                        value = z_mat$Z99)

ggpubr::ggboxplot(data = plot_data, x = "group", y = "value", palette = "aaas", fill = "group" , main = "Z99 Stratifying HF and Control", ylab = "Z99 Value") 
                  + ggpubr::stat_compare_means(label = "p.signif")
