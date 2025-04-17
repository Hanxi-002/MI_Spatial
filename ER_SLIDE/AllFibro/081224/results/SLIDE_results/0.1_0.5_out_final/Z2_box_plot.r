library(ggplot2)

z_mat <- read.csv("ER_SLIDE/AllFibro/081224/results/SLIDE_results/0.1_0.5_out_final/z_matrix.csv", row.names = 1)
y <- read.csv("ER_SLIDE/AllFibro/081224/Data/y.csv", row.names = 1)

y$V1 <- ifelse(y$V1==0, "Control", "HF")

plot_data <- data.frame(group = as.factor(y$V1),
                        value = z_mat$Z2)


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
  labs(title = "Z2 Stratifying HF and Control", y = "Z2 Value") +
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
ggsave('ER_SLIDE/AllFibro/081224/results/SLIDE_results/0.1_0.5_out_final/Z2_box_plot_20250415.pdf', p, width = 4, height = 5, units = 'in')
