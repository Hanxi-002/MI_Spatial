library(Seurat)
library(ggplot2)
source('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Calc_Z_Hat_Helper.R')

####################################################### Z 37 #############################################################################
####################################################################################################################################
# --------------------------------------------
# Check Z37 sign with GeoMx Data
# --------------------------------------------

y_path = '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/Data/ccr2_af_y.csv'
z_path = '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/z_matrix.csv'
z_score_column = 37
custom_labels = c('close to CCR2- macs', 'close to CCR2+ macs')
find_LF_sign(y_path, z_path, z_score_column, custom_labels)

####################################################################################################################################
# ----------------------------------------------------------
# Calculate Z values for the LF of interest in the new data
# ----------------------------------------------------------
seurat_path <- '/ix/djishnu/Mary/MI_data/Kramann_2023_Visium/data/Kramann_fibroblasts.rds'
er_results_path <- '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/AllLatentFactors.rds'
z_score_column = 37
plot_title = "Z37 Score by Condition (Kramann Fibroblast)"
# custom_order = c("CTRL", "BZ", "FZ", "IZ", "RZ")
custom_order = c("FZ", "IZ")

results <- main(seurat_path, er_results_path, z_score_column, plot_title, custom_order = custom_order)
write.csv(results$plot$data, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Fibro/Z37/Z37_Score_FZ_IZ_(Kramann_Fibroblast).csv')
ggsave('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Fibro/Z37/Z37_Score_FZ_IZ_(Kramann_Fibroblast).pdf', results$plot, width = 4, height = 6)

# --------------------------------------------------------------
# Mann-Whitney U test on the groups of the new z scores (z_hat)
# --------------------------------------------------------------
score <-read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Fibro/Z37/Z37_Score_FZ_IZ_(Kramann_Fibroblast).csv', row.names = 1)
df = perform_mw_tests(score)
#write.csv(df, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Fibro/Z37/FZ_IZ_Mann_Whitney_P_val.csv')

# --------------------------------------------------------------
# Cliff's Delta Calculation
# --------------------------------------------------------------
#score <-read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Fibro/Z37/Z37_Score_by_Condition_(Kramann_Fibroblast).csv', row.names = 1)
score <-read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Fibro/Z37/Z37_Score_FZ_IZ_(Kramann_Fibroblast).csv', row.names = 1)
results <- perform_cliffs_delta(score)

sig <- cbind(df, results)
write.csv(sig, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Fibro/Z37/score_box_significance.csv')


####################################################################################################################################
# ------------------------------------------------------------------
# Replotting box plot to bar plot which shows # of samples above Q3
# ------------------------------------------------------------------
score <-read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Fibro/Z37/Z37_Score_FZ_IZ_(Kramann_Fibroblast).csv', row.names = 1)
dim(score[score$condition == 'FZ', ] ) # 6096
dim(score[score$condition == 'IZ', ] ) # 5868

q3_ratio <- counts_above_q3(score, baseline = 'IZ')

# ------------------------------------------------------------------
# Clopper Pearson Confidence Interval
# ------------------------------------------------------------------

q3_ratio <- Calc_Clopper_Pearson_CIs(q3_ratio,count_col = 'count_above_q3', total_col = 'total_count')
write.csv(q3_ratio, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Fibro/Z37/q3_ratio.csv')

# ------------------------------------------------------------------
# Exact Binomial Test
# ------------------------------------------------------------------
res = exact_binomial_test(df = q3_ratio, count_col = 'count_above_q3', condition_col = 'condition', baseline = 'IZ', total_col = 'total_count')
write.csv(res, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Fibro/Z37/q3_significance.csv')

# ------------------------------------------------------------------
# Plot Q3 Bar Plot
# ------------------------------------------------------------------
q3_bar <- plot_ratios(q3_ratio)
ggsave('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Fibro/Z37/q3_bar_plot.pdf', q3_bar, width = 4, height = 6)


############################################################# Z 17 #######################################################################
####################################################################################################################################
seurat_path <- '/ix/djishnu/Mary/MI_data/Kramann_2023_Visium/data/Kramann_fibroblasts.rds'
er_results_path <- '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/AllLatentFactors.rds'
z_score_column = 17
plot_title = "Z17 Score by Condition (Kramann Fibroblast)"
results <- main(seurat_path, er_results_path, z_score_column, plot_title)
write.csv(results$plot$data, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Fibro/Z17_Score_by_Condition_(Kramann_Fibroblast).csv')




