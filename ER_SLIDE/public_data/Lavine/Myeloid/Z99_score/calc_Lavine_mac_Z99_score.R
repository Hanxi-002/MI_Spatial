library(Seurat)
library(ggplot2)
source('/ocean/projects/cis240075p/hxiao2/MI_Spatial/ER_SLIDE/public_data/Calc_Z_Hat_Helper.R')

####################################################################################################################################
# --------------------------------------------
# Check Z99 sign with GeoMx Data
# -------------------------------------------
y_path = '/ocean/projects/cis240075p/hxiao2/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/Data/y.csv'
z_path = '/ocean/projects/cis240075p/hxiao2/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/results/SLIDE_Results/z_matrix.csv'
z_score_column = 99
custom_labels = c('close to RF', 'close to AF')
find_LF_sign(y_path, z_path, z_score_column, custom_labels)

####################################################################################################################################
# seurat_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Macrophage_Seurat_obj.RDS"
# data <- readRDS(seurat_path)
# data <- subset(x = data, subset = HF.etiology %in% c('Donor', 'ICM', 'AMI'))
# saveRDS(data,  "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Mac_Seurat_obj_no_NICM.RDS")

# ----------------------------------------------------------
# Calculate Z values for the LF of interest in the new data
# ----------------------------------------------------------
seurat_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Z99_score/Mac_Seurat_obj_no_NICM.RDS"
er_results_path <- '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/results/ER_Results/final_delta_0.001_lambda_0.5.rds'
z_score_column = 99
plot_title = "Z99 Score by Condition (Lavine Macrophage)"


custom_order = c("Donor", "AMI", "ICM")
results <- main_Lavine(seurat_path, er_results_path, z_score_column, plot_title, custom_order = custom_order)
write.csv(results$plot$data, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Z99_Score_by_Condition_(Lavine_Macrophage).csv')
ggsave('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Z99_score/Z99_Score_by_Condition_(Lavine_Mac).pdf', results$plot, width = 4, height = 6)
# --------------------------------------------------------------
# Mann-Whitney U test on the groups of the new z scores (z_hat)
# --------------------------------------------------------------
score <-read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Z99_score/Z99_Score_by_Condition_(Lavine_Macrophage).csv', row.names = 1)
df = perform_mw_tests(score)

# --------------------------------------------------------------
# Cliff's Delta Calculation
# --------------------------------------------------------------
score <-read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Z99_score/Z99_Score_by_Condition_(Lavine_Macrophage).csv', row.names = 1)
results <- perform_cliffs_delta(score)

sig <- cbind(df, results)
write.csv(sig, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Z99_score/score_box_significance.csv')

####################################################################################################################################
# ------------------------------------------------------------------
# Replotting box plot to bar plot which shows # of samples above Q3
# ------------------------------------------------------------------
score <-read.csv('//ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Z99_score/box_plot/Z99_Score_by_Condition_(Lavine_Macrophage).csv', row.names = 1)
q3_ratio <- counts_above_q3(score, baseline = 'Donor')

# ------------------------------------------------------------------
# Clopper Pearson Confidence Interval
# ------------------------------------------------------------------

q3_ratio <- calc_SEM_CIs(q3_ratio)
#q3_ratio <- calc_Clopper_Pearson_CIs(q3_ratio,count_col = 'count_above_q3', total_col = 'total_count')
write.csv(q3_ratio, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Z99_score/bar_plot/q3_ratio.csv')

# ------------------------------------------------------------------
# Proportional Test
# ------------------------------------------------------------------

#res = exact_binomial_test(df = q3_ratio, count_col = 'count_above_q3', condition_col = 'condition', baseline = 'IZ', total_col = 'total_count')
res = control_prop_test(df = q3_ratio, props_col = 'ratio_above_q3', total_col = 'total_count', condition_col = 'condition', baseline = 'Donor')
write.csv(res, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Z99_score/q3_significance.csv')

# ------------------------------------------------------------------
# Plot Q3 Bar Plot
# ------------------------------------------------------------------
q3_bar <- plot_ratios(q3_ratio, custom_order = c('Donor', 'AMI', 'ICM'))
ggsave('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Z99_score/bar_plot/q3_bar_plot.pdf', q3_bar, width = 4, height = 6)






