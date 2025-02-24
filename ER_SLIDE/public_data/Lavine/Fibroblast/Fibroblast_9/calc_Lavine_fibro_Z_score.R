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

# seurat_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/fibro_9.RDS"
# data <- readRDS(seurat_path)
# data <- subset(x = data, subset = HF.etiology %in% c('Donor', 'ICM', 'AMI'))
# saveRDS(data,  "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/fibro_9_no_NICM.RDS")

# ----------------------------------------------------------
# Calculate Z values for the LF of interest in the new data
# ----------------------------------------------------------
seurat_path <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/fibro_9_no_NICM.RDS"
er_results_path <- '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/AllLatentFactors.rds'
z_score_column = 37
plot_title = "Z37 Score by Condition (Lavine Fibroblast_9)"

custom_order = c("Donor", "AMI", "ICM")
results <- main_Lavine(seurat_path, er_results_path, z_score_column, plot_title, custom_order = custom_order)
write.csv(results$plot$data,'/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/Z37/Z37_Score_by_Condition_(Lavine_Fibroblast_9).csv')

# --------------------------------------------------------------
# Mann-Whitney U test on the groups of the new z scores (z_hat)
# --------------------------------------------------------------
score <-read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/Z37/Z37_Score_by_Condition_(Lavine_Fibroblast_9).csv', row.names = 1)
df = perform_mw_tests(score)
write.csv(df, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/Z37/Mann_Whitney_P_val.csv')

# --------------------------------------------------------------
# Cliff's Delta Calculation
# --------------------------------------------------------------
score <-read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/Z37/Z37_Score_by_Condition_(Lavine_Fibroblast_9).csv', row.names = 1)
results <- perform_cliffs_delta(score)
write.csv(results, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/Z37/Cliffs_Delta.csv')


####################################################################################################################################
# ------------------------------------------------------------------
# Replotting box plot to bar plot which shows # of samples above Q3
# ------------------------------------------------------------------
score <-read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/Z37/Z37_Score_by_Condition_(Lavine_Fibroblast_9).csv', row.names = 1)
q3_ratio <- counts_above_q3(score, baseline = 'Donor')
plot_ratios(q3_ratio, custom_order = c('Donor', 'AMI', 'ICM'))
write.csv(q3_ratio, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/Z37/q3_ratio.csv')

# ------------------------------------------------------------------
# Calculate exact binomal test
# ------------------------------------------------------------------
# Pairwised comparison, but only with each condition vs control. 
# p values not adjusted
q3_ratio = read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/Z37/q3_ratio.csv', row.names = 1)

counts = q3_ratio$count_above_q3
ns = q3_ratio$total_count
conditions = q3_ratio$condition
res = exact_binomial_test(counts, ns, conditions, baseline = 'Donor')
write.csv(res, '/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Fibroblast/Fibroblast_9/Z37/q3_ratio_significance.csv')


