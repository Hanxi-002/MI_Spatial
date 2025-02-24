library(Seurat)
library(ggplot2)
source('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Calc_Z_Hat_Helper.R')


####################################################################################################################################
# --------------------------------------------
# Check Z99 sign with GeoMx Data
# --------------------------------------------
y_path = '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/Data/y.csv'
z_path = '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/results/SLIDE_Results/z_matrix.csv'
z_score_column = 99
custom_labels = c('close to RF', 'close to AF')
find_LF_sign(y_path, z_path, z_score_column, custom_labels)


####################################################################################################################################
# ----------------------------------------------------------
# Calculate Z values for the LF of interest in the new data
# ----------------------------------------------------------
seurat_path <- "/ix/djishnu/Mary/MI_data/Kramann_2023_Visium/data/Kramann_macrophages.rds"
er_results_path <- '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/results/ER_Results/final_delta_0.001_lambda_0.5.rds'
z_score_column = 99
plot_title = "Z99 Score by Condition (Kramann Macrophages)"


custom_order = c("IZ", "FZ")
results <- main(seurat_path, er_results_path, z_score_column, plot_title, custom_order = custom_order)
write.csv(results$plot$data, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Macrophages/Z99_Score_by_Condition_(Kramann_Mac).csv')


# --------------------------------------------------------------
# Mann-Whitney U test on the groups of the new z scores (z_hat)
# --------------------------------------------------------------
score <-read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Macrophages/Z99_Score_by_Condition_(Kramann_Mac).csv', row.names = 1)
df = perform_mw_tests(score)
write.csv(df, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Macrophages/Mann_Whitney_P_val.csv')

# --------------------------------------------------------------
# Cliff's Delta Calculation
# --------------------------------------------------------------
score <-read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Macrophages/Z99_Score_by_Condition_(Kramann_Mac).csv', row.names = 1)
results <- perform_cliffs_delta(score)
write.csv(results, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Macrophages/Cliffs_Delta.csv')


####################################################################################################################################
# ------------------------------------------------------------------
# Replotting box plot to bar plot which shows # of samples above Q3
# ------------------------------------------------------------------
# Create the bar plot

score <-read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Macrophages/Z99_Score_by_Condition_(Kramann_Mac).csv', row.names = 1)
# nrow(score[score$condition == 'IZ', ]) # FZ v IZ =  29 v 4346

set.seed(123) 
FZ_score <- score[score$condition == 'FZ', ]
IZ_score <- score[score$condition == 'IZ', ][sample(sum(score$condition == 'IZ'), nrow(FZ_score)), ]
score_new <- rbind(FZ_score, IZ_score)

q3_ratio <- counts_above_q3(score_new, baseline = 'IZ') # which one should I use for baseline.
plot_ratios(q3_ratio, custom_order = c('FZ', 'IZ'))
write.csv(q3_ratio, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Macrophages/q3_ratio.csv')


# ------------------------------------------------------------------
# Calculate exact binomial test
# ------------------------------------------------------------------
# Pairwised comparison, but only with each condition vs control. 
# p values not adjusted

counts = q3_ratio$count_above_q3
ns = q3_ratio$total_count
conditions = q3_ratio$condition
res = exact_binomial_test(counts, ns, conditions, baseline = 'IZ')

write.csv(res, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Macrophages/q3_ratio_significance.csv')

