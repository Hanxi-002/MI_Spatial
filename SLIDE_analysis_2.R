library(doParallel)
library(dplyr)
library(pROC)
library(SLIDE)

final_res <- NULL
x <- as.matrix(utils::read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/HF/RestVActive/022823/Data/x.csv", row.names = 1))
z <- as.matrix(utils::read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/HF/CdVRest/022823/Results/z_matrix.csv", row.names = 1)) ## standardized
y <- as.matrix(utils::read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/HF/RestVActive/022823/Data/y.csv", row.names = 1))
er_res <- readRDS('/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/HF/CdVRest/022823/Results/final_delta_0.01_lambda_0.5.rds')


SLIDE_res <- SLIDE(z, y, method = 4, do_interacts = TRUE, betas = NULL, top_prop = NULL, marginals = NULL,
                   spec = 0.1, fdr = 0.1, niter = 500, elbow = FALSE, f_size = 15, parallel = TRUE, ncore = 10)

print(SLIDE_res$marginal_vars)
print(SLIDE_res$interaction_vars)
#feature_res <- findFeats(x, y, z, er_res, thresh = 0.01, p_thresh = 0.1)
SLIDE_param <- c(4, 0.1, 0.1, 500, 15)
names(SLIDE_param) <- c("method", "spec", "fdr", "niter", "f_size")

ks <- c(2, 8, 10, 12, 4, 7, 9, 12)

A <- er_res$A[, ks]
gene_names <- colnames(x)
temp <- NULL
for (i in 1:ncol(A)){
  AUCs <- c()
  signs <- c()
  corrs <- c()
  idx <- which(A[, i] != 0)
  A_loading <- abs(A[, i][-which(A[, i] == 0)])
  names <- gene_names[idx]
  for (j in 1:length(idx)) { ## loop through variables with big enough loadings
    corr <- cor(x[,idx[j]],y,method = "spearman")
    sign <- sign(cor(x[,idx[j]],y,method = "spearman"))
    AUC <- auc(y, x[, idx[j]])
    AUCs <- c(AUCs, AUC)
    corrs <- c(corrs, corr)
    signs <- c(signs, sign)
  }
  color <- recode(signs, "-1" = "Blue", "1"= "Red")
  df <- data.frame(names, A_loading, AUCs, corrs, color)
  df <- df[order(-df$A_loading), ]
  top <- df[1:10, ]
  df <- df[order(-df$AUCs), ]
  bot <- df[1:10, ]
  final <- unique(rbind(top, bot))
  temp[[i]] <- final
  write.table(final, file = paste0("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/HF/CdVRest/022823/Results/_gene_list_",
                                   colnames(A)[i], ".txt"), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
}

names(temp) <- colnames(A)
final_res$SLIDE_res <- SLIDE_res
final_res$feature_res <- temp
final_res$SLIDE_param <- SLIDE_param
#saveRDS(final_res, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/HF/CdVRest/022823/Results/SLIDE_res.rds")


toy_1 <- c()
toy_2 <- c()
for (i in 1:ncol(x)){
  toy_1 <- append(toy_1, sign(cor(x[ ,i], y, method = 'spearman')))
  toy_2 <- append(toy_2, sign(cor.test(x[,i], y)$estimate))
}
sum(toy_2)
length(which(toy_1 == 1))
length(which(toy_2 == 1))
which(toy_1 != toy_2)


