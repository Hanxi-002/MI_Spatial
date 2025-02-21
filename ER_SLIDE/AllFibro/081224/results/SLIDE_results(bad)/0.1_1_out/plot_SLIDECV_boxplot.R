perRes <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllFibro/081224/results/SLIDE_results/0.01_1_out/SLIDECV_boxplot_data.rds")

lambda_boxplot = ggpubr::ggboxplot(data = perRes, x = "method", y = "auc", palette = "aaas",
                                   fill = "method" ) +
  ggpubr::stat_compare_means(label = "p.signif")
ggplot2::ggsave(plot = lambda_boxplot, filename = paste0("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllFibro/081224/results/SLIDE_results/0.1_1_out", "/SLIDECV_boxplot.pdf"), height = 6, width = 6)
