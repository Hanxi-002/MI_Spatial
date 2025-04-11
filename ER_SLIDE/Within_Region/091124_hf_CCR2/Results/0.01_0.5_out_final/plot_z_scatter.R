library(ggplot2)

z_mat = read.csv("/ocean/projects/cis240075p/hxiao2/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/z_matrix.csv", row.names = 1)
y = read.csv('/ocean/projects/cis240075p/hxiao2/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/Data/ccr2_af_y.csv', row.names = 1)


z_mat = z_mat[ ,colnames(z_mat) %in% c('Z37', 'Z17')]
z_mat['color'] = y$ccr2
z_mat$color <- ifelse(z_mat$color == 0, 'Control', "HF")

#ggplot(z_mat, aes(x = Z39, y = Z44, color = color)) + geom_point() + labs(title = 'Activated Fibroblasts', y = 'Z44', x = 'Z39') + theme_minimal()

p <- ggplot(z_mat, aes(x = Z37, y = Z17, color = color)) + 
  geom_point() + 
  labs(title = 'AF proximal v distal to CCR2+ Macs', y = 'Z17', x = 'Z37') + 
  scale_color_manual(values = c("Control" = "#0000FF", "HF" = "#FF0000")) +
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),         # Remove grid lines
    panel.background = element_blank(),   # Remove gray background
    axis.line = element_line(color = "black") # Add axis lines
  )

p <- ggsave("/ocean/projects/cis240075p/hxiao2/MI_Spatial/ER_SLIDE/Within_Region/091124_hf_CCR2/Results/0.01_0.5_out_final/Z37_Z17_scatter_20250411.pdf", plot = p, width = 4, height = 5, units = "in")
