source("/ix/djishnu/Hanxi/Common_R/CalcCliffDelta_Helper.R")

z_matrix = read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/results/SLIDE_Results/z_matrix.csv", row.names = 1)
interacting_z = z_matrix[ , colnames(z_matrix)%in% c("Z45", "Z131")]

y = read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/Data/y.csv", row.names = 1)

comb = GetPairwiseComb(y[, 1])

all_cds = CalcCliffDelta(interacting_z, y, comb, sig_idx = c('Z45', "Z131"))

# from all_cds, we see that LF1 has the highest Cliff's Delta
# plot Z1 with Z2

z_mat = z_matrix
z_mat = z_mat[ ,colnames(z_mat) %in% c('Z45', 'Z99')]
z_mat['color'] = y[, 1]
z_mat$color <- ifelse(z_mat$color == 0, 'Control', "HF")

#ggplot(z_mat, aes(x = Z39, y = Z44, color = color)) + geom_point() + labs(title = 'Activated Fibroblasts', y = 'Z44', x = 'Z39') + theme_minimal()

ggplot(z_mat, aes(x = Z99, y = Z45, color = color)) + 
  geom_point() + 
  labs(title = 'Macs distal v proximal to AFs', y = 'Z45', x = 'Z99') + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),         # Remove grid lines
    panel.background = element_blank(),   # Remove gray background
    axis.line = element_line(color = "black") # Add axis lines
  )
