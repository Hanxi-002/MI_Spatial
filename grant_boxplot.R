SLIDE_res <- readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/AllCell/022423/Results/SLIDE_res.rds")
x <- as.matrix(utils::read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/AllCell/022423/Data/x.csv", row.names = 1))
y <- as.matrix(utils::read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/AllCell/022423/Data/y.csv", row.names = 1))
ER_res <- readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/AllCell/022423/Results/final_delta_0.01_lambda_0.5.rds")
#idx <- which(A[, "Z12"] != 0)
#A_loading <- abs(A[, "Z12"][-which(A[, "Z12"] == 0)])
#gene_names <- colnames(x)
#names <- gene_names[idx]
#x <- x[ , idx]

gene_names <- c("PABPC1", "MARCKS", "SAFB", "GRN", "RNASE1", "AP5MG", "COL1A1", "PLXND1", "MCL1", "CTSC", "SPARC", "VIM", "FOSL2", "CD74", "RBM5", 
"LRP1", "PCBP2")

x <- x[ , colnames(x) %in% gene_names]
pos_idx <- which(y == 1)
neg_idx <- which(y == 0)

pos_x <- x[pos_idx, ]
neg_x <- x[neg_idx, ]

pos_x_exp <- apply(pos_x,2,function(x) median(x[x!=0]))
neg_x_exp <- apply(neg_x,2,function(x) median(x[x!=0]))

final_df <- cbind(neg_x_exp, pos_x_exp)
toy <- melt(final_df)

ggp <- ggplot(toy, aes(X2, X1)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low="white", high="blue") +
  theme_bw()
ggp                                                               # Print heatmap


