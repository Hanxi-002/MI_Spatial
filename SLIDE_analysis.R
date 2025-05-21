library(doParallel)
library(dplyr)
library(pROC)

files <- list.files("/Users/xiaoh/Desktop/Research/Hierarchical_ER/R")
source_files <- paste0("/Users/xiaoh/Desktop/Research/Hierarchical_ER/R/", files)
sapply(source_files, source)


final_res <- NULL
x <- as.matrix(utils::read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/CD68/022423/Data/x.csv", row.names = 1))
z <- as.matrix(utils::read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/CD68/022423/Results/z_matrix.csv", row.names = 1)) ## standardized
y <- as.matrix(utils::read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/CD68/022423/Data/y.csv", row.names = 1))
er_res <- readRDS('/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/CD68/022423/Results/final_er_output.rds')


SLIDE_res <- SLIDE(z, y, method = 4, do_interacts = TRUE, betas = NULL, top_prop = NULL, marginals = NULL,
                   spec = 0.1, fdr = 0.1, niter = 500, elbow = FALSE, f_size = 22, parallel = TRUE, ncore = 10)

print(SLIDE_res$marginal_vars)
print(SLIDE_res$interaction_vars)
#feature_res <- findFeats(x, y, z, er_res, thresh = 0.01, p_thresh = 0.1)
SLIDE_param <- c(4, 0.1, 0.1, 500, 22)
names(SLIDE_param) <- c("method", "spec", "fdr", "niter", "f_size")


ks <- c(12, 17, 4, 6, 22)
# temp <- NULL
# #names <-
# for (i in 1:length(ks)){
#   feature_res <-SigGenes_A(er_res, ks[[i]], as.matrix(x), matrix(y,ncol = 1), thre=0.05, negcol="blue", posCol="red",orderbycol=T)
#   feature_res <- feature_res[order(feature_res[, 1], decreasing = TRUE), ]
#   temp[[i]] <- feature_res
# }
# names(temp) <- ks


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
  # write.table(final, file = paste0("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/CD68/011123/Results2/_gene_list_",
  #                                  colnames(A)[i], ".txt"), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
}

names(temp) <- colnames(A)
final_res$SLIDE_res <- SLIDE_res
final_res$feature_res <- temp
final_res$SLIDE_param <- SLIDE_param
#saveRDS(final_res, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/HF/RestVActive/011322/Restuls2/SLIDE_res.rds")

