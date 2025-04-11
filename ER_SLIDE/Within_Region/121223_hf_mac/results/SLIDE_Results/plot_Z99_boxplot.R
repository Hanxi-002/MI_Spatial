library(ggplot2)

z_mat <- read.csv("/ocean/projects/cis240075p/hxiao2/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/results/SLIDE_Results/z_matrix.csv", row.names = 1)
y <- read.csv("/ocean/projects/cis240075p/hxiao2/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/Data/y.csv", row.names = 1)

y$y.row.names.y...in..row.names.CD68_HF.... <- ifelse(y$y.row.names.y...in..row.names.CD68_HF....==0, "Control", "HF")

plot_data <- data.frame(group = as.factor(y$y.row.names.y...in..row.names.CD68_HF....),
                        value = z_mat$Z99)


test_result <- t.test(value ~ group, data = plot_data)
p_value <- test_result$p.value

# Create significance label
signif_label <- if(p_value < 0.001) {
  "***"
} else if(p_value < 0.01) {
  "**"
} else if(p_value < 0.05) {
  "*"
} else {
  "ns"
}

# Get the y-position for the significance bar
y_max <- max(plot_data$value, na.rm = TRUE)
y_pos <- y_max + (y_max * 0.1)  # Place it 10% above the max value

# Create the plot
p <- ggplot(plot_data, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#0000FF", "#FF0000")) +  # AAAS palette colors
  labs(title = "Z99 Stratifying HF and Control", y = "Z99 Value") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),         
    panel.background = element_blank(),   
    axis.line = element_line(color = "black")
    ) +
  # Add significance annotation
  annotate("text", x = 1.5, y = y_pos, label = signif_label) +
  annotate("segment", x = 1, xend = 2, y = y_pos * 0.95, yend = y_pos * 0.95) +
  theme(legend.position = "right")

print(p)
ggsave('/ocean/projects/cis240075p/hxiao2/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/results/SLIDE_Results/Z99_box_plot_20250411.pdf', p, width = 4, height = 5, units = 'in')

