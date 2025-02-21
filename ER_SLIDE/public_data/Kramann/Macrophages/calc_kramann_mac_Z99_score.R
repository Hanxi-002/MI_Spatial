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
write.csv(df, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Macrophages/Cliffs_Delta.csv')


####################################################################################################################################
# ------------------------------------------------------------------
# Replotting box plot to bar plot which shows # of samples above Q3
# ------------------------------------------------------------------
# Create the bar plot

score <-read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Kramann/Macrophages/Z99_Score_by_Condition_(Kramann_Mac).csv', row.names = 1)
results <- counts_above_q3(score)
plot_counts(results, custom_order = NULL)
