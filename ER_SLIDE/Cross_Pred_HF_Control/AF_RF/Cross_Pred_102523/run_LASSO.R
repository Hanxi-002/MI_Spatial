setwd("/ix/djishnu/Hanxi/CK_COPD_6MonB6/CrossPred")
#library(glmnet)
library(reshape)
source("LASSO.R")

#load the LFs that are used for cross prediction
SLIDE_res <- readRDS("/ix/djishnu/Hanxi/CK_COPD_6MonB6/AT2/JS_072023/SLIDE_Results/spec_0.2_/plotSigGenes_data.RDS")
#load the data that are being predicted
x <- read.csv("/ix/djishnu/Hanxi/CK_COPD_6MonB6/CC/JS_072023/Data/CC_JS_x.csv", row.names = 1)
y <- read.csv("/ix/djishnu/Hanxi/CK_COPD_6MonB6/CC/JS_072023/Data/CC_JS_y.csv", row.names = 1)

#subset data
LF = SLIDE_res[SLIDE_res$lf_num == 2, ]
LF_x = x[, colnames(x) %in% LF$names]
dim(x)
dim(LF_x)
#colnames(x)

#LASSO
data <- as.data.frame(scale(LF_x, T, T))
data["Y"] <- scale(y, T, T)
lasso_res <- VanillaLasso(30, 5, data)
write.csv(lasso_res[[1]], "/ix/djishnu/Hanxi/CK_COPD_6MonB6/CrossPred/AT2_CC/Justin/LASSO_101823/auc.csv")
write.csv(lasso_res[[2]], "/ix/djishnu/Hanxi/CK_COPD_6MonB6/CrossPred/AT2_CC/Justin/LASSO_101823/pred_y.csv")
write.csv(lasso_res[[3]], "/ix/djishnu/Hanxi/CK_COPD_6MonB6/CrossPred/AT2_CC/Justin/LASSO_101823/features.csv")
#plotting
auc <- lasso_res[[1]]
meltData <- melt(auc)
p <- ggplot(meltData, aes(x = variable, y = value)) + geom_boxplot()
p
ggsave("/ix/djishnu/Hanxi/CK_COPD_6MonB6/CrossPred/AT2_CC/Justin/LASSO_101823/performance.pdf")
