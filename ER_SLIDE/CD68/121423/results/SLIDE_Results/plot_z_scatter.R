library(ggplot2)

z_mat = read.csv("/ocean/projects/cis240075p/hxiao2/MI_Spatial/ER_SLIDE/CD68/121423/results/SLIDE_Results/z_matrix.csv", row.names = 1)
y = read.csv('/ocean/projects/cis240075p/hxiao2/MI_Spatial/ER_SLIDE/CD68/121423/Data/y.csv', row.names = 1)


z_mat = z_mat[ ,colnames(z_mat) %in% c('Z12', 'Z17')]
z_mat['color'] = y$V1
z_mat$color <- ifelse(z_mat$color == 0, 'Control', "HF")

p <- ggplot(z_mat, aes(x = Z17, y = Z12, color = color)) + 
  geom_point() + 
  labs(title = 'Macrophages', y = 'Z17', x = 'Z12') + 
  scale_color_manual(values = c("Control" = "#0000FF", "HF" = "#FF0000")) +
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),         # Remove grid lines
    panel.background = element_blank(),   # Remove gray background
    axis.line = element_line(color = "black") # Add axis lines
  )

ggsave("/ocean/projects/cis240075p/hxiao2/MI_Spatial/ER_SLIDE/CD68/121423/results/SLIDE_Results/ Z12_Z17_scatter_V2.pdf", plot = p, width = 4, height = 5, units = "in")