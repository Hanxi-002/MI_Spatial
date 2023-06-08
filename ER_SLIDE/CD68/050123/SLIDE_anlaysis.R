library(SLIDEHelper)
y_path <- "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/CD68/050123/Data/y.csv"
z_path <- "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/CD68/050123/ER_Results/z_matrix.csv"
er_path <- "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/CD68/050123/ER_Results/final_delta_0.01_lambda_0.5.rds"
method <- 4
do_interacts <- TRUE
spec <- 0.1
fdr <- 0.1
niter <- 500



Z_matrix <- CalcZMatrix(x_path, er_path, out_path)

SLIDE_res <- runSLIDE(y_path = y_path,
                      z_path = z_path,
                      er_path = er_path, 
                      do_interacts = do_interacts,
                      method = method,
                      spec = spec,
                      fdr = fdr,
                      niter = niter)


x_path <- "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/CD68/050123/Data/x.csv"
out_path <- "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/CD68/050123/SLIDE_Results"
SLIDE_res <- SLIDE_res
num_top_feats <- 10
condition <- "corr"

SLIDE_res <- GetTopFeatures(x_path, y_path, er_path, out_path, SLIDE_res, num_top_feats = 10, condition = "corr")


