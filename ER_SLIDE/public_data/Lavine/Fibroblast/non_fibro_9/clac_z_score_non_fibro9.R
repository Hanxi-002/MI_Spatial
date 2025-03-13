library(Seurat)
library(ggplot2)
source('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Calc_Z_Hat_Helper.R')

# ----------------------------------------------------------
# Calculate Z values for the LF of interest in the new data
# ----------------------------------------------------------
seurat_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/non_fibro_9/non_fibro_9.RDS"
er_results_path <- '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/AllLatentFactors.rds'
z_score_column = 37
plot_title = "Z37 Score by Condition (Lavine Fibroblast_9)"

custom_order = c("Donor", "AMI", "ICM")
results <- main_Lavine(seurat_path, er_results_path, z_score_column, plot_title, custom_order = custom_order)
write.csv(results$plot$data,'/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/non_fibro_9/box_plot/Z37_Score_by_Condition_(Lavine_NonFibro_9).csv')
ggsave('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/non_fibro_9/box_plot/Z37_Score_by_Conition_(Lavine_NonFibro_9).pdf', results$plot, width = 4, height = 6)

# --------------------------------------------------------------
# Mann-Whitney U test on the groups of the new z scores (z_hat)
# --------------------------------------------------------------
score <-read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/non_fibro_9/box_plot/Z37_Score_by_Condition_(Lavine_NonFibro_9).csv', row.names = 1)
df = perform_mw_tests(score)

# --------------------------------------------------------------
# Cliff's Delta Calculation
# --------------------------------------------------------------

results <- perform_cliffs_delta(score)

sig <- cbind(df, results)
write.csv(sig, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/non_fibro_9/box_plot/score_box_significance.csv')

####################################################################################################################################
# ------------------------------------------------------------------
# Replotting box plot to bar plot which shows # of samples above Q3
# ------------------------------------------------------------------
score <-read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/non_fibro_9/box_plot/Z37_Score_by_Condition_(Lavine_NonFibro_9).csv', row.names = 1)
q3_ratio <- counts_above_q3(score, baseline = 'Donor')

# ------------------------------------------------------------------
# Clopper Pearson Confidence Interval
# ------------------------------------------------------------------

q3_ratio <- calc_SEM_CIs(q3_ratio)
#q3_ratio <- calc_Clopper_Pearson_CIs(q3_ratio,count_col = 'count_above_q3', total_col = 'total_count')
write.csv(q3_ratio, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/non_fibro_9/bar_plot/q3_ratio.csv')

# ------------------------------------------------------------------
# Proportional Test
# ------------------------------------------------------------------

#res = exact_binomial_test(df = q3_ratio, count_col = 'count_above_q3', condition_col = 'condition', baseline = 'IZ', total_col = 'total_count')
res = control_prop_test(df = q3_ratio, props_col = 'ratio_above_q3', total_col = 'total_count', condition_col = 'condition', baseline = 'Donor')
write.csv(res, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/non_fibro_9/bar_plot/q3_significance.csv')

# ------------------------------------------------------------------
# Plot Q3 Bar Plot
# ------------------------------------------------------------------
q3_bar <- plot_ratios(q3_ratio, custom_order = c('Donor', 'AMI', 'ICM'))
ggsave('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/non_fibro_9/bar_plot/q3_bar_plot.pdf', q3_bar, width = 4, height = 6)
