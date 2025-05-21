source("/ix/djishnu/Hanxi/Common_R/CalcCliffDelta_Helper.R")

z_matrix = read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllFibro/081224/results/SLIDE_results/0.1_0.5_out_final/z_matrix.csv", row.names = 1)
interacting_z = z_matrix[ , colnames(z_matrix)%in% c("Z1", "Z4", "Z3")]

y = read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/AllFibro/081224/Data/y.csv", row.names = 1)

comb = GetPairwiseComb(y[, 1])

all_cds = CalcCliffDelta(interacting_z, y, comb, sig_idx = NULL)

# from all_cds, we see that LF1 has the highest Cliff's Delta
# plot Z1 with Z2

z_mat = z_matrix

z_mat = z_mat[ ,colnames(z_mat) %in% c('Z1', 'Z2')]
z_mat['color'] = y$V1
z_mat$color <- ifelse(z_mat$color == 0, 'Control', "HF")

#ggplot(z_mat, aes(x = Z39, y = Z44, color = color)) + geom_point() + labs(title = 'Activated Fibroblasts', y = 'Z44', x = 'Z39') + theme_minimal()

ggplot(z_mat, aes(x = Z1, y = Z2, color = color)) + 
  geom_point() + 
  labs(title = 'All Fibroblasts', y = 'Z2', x = 'Z1') + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),         # Remove grid lines
    panel.background = element_blank(),   # Remove gray background
    axis.line = element_line(color = "black") # Add axis lines
  )
