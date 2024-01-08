setwd("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Cross_Pred_HF_Control/RF_AF/CrossPred_103023")
#library(glmnet)
library(reshape)
source("/ix/djishnu/Hanxi/Common_R/LASSO.R")

#load the LFs that are used for cross prediction
SLIDE_res <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/RestingFibro/022423/Results/SLIDE_Results_071223/plotSigGenes_data.RDS")
#load the data that are being predicted
x <- read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/ActiveFibro/022423/x.csv", row.names = 1)
y <- read.csv("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/ActiveFibro/022423/y.csv", row.names = 1)

#subset data
LF = SLIDE_res[SLIDE_res$lf_num %in% c(22, 18, 2), ]
LF_x = x[, colnames(x) %in% LF$names]
dim(x)
dim(LF_x)
#colnames(x)

#LASSO
data <- as.data.frame(scale(LF_x, T, T))
data["Y"] <- scale(y, T, T)
lasso_res <- VanillaLasso(30, 5, data)
outpath <- "/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Cross_Pred_HF_Control/RF_AF/CrossPred_103023/"
write.csv(lasso_res[[1]], paste0(outpath, "auc.csv"))
write.csv(lasso_res[[2]], paste0(outpath, "pred_y.csv"))
write.csv(lasso_res[[3]], paste0(outpath, "features.csv"))

#plotting
auc <- lasso_res[[1]]
meltData <- melt(auc)
p <- ggplot(meltData, aes(x = variable, y = value)) + geom_boxplot()
p
ggsave(paste0(outpath, "performance.pdf"))
