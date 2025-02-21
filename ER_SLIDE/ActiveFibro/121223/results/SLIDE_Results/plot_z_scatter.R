library(ggplot2)

z_mat = read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/ActiveFibro/121223/results/SLIDE_Results/z_matrix.csv", row.names = 1)
y = read.csv('//ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/ActiveFibro/121223/Data/y.csv', row.names = 1)


z_mat = z_mat[ ,colnames(z_mat) %in% c('Z44', 'Z39')]
z_mat['color'] = y$V1
z_mat$color <- ifelse(z_mat$color == 0, 'Control', "HF")

#ggplot(z_mat, aes(x = Z39, y = Z44, color = color)) + geom_point() + labs(title = 'Activated Fibroblasts', y = 'Z44', x = 'Z39') + theme_minimal()

ggplot(z_mat, aes(x = Z39, y = Z44, color = color)) + 
  geom_point() + 
  labs(title = 'Activated Fibroblasts', y = 'Z44', x = 'Z39') + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),         # Remove grid lines
    panel.background = element_blank(),   # Remove gray background
    axis.line = element_line(color = "black") # Add axis lines
  )
