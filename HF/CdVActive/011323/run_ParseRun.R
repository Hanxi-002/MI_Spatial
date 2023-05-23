library(EssReg)
library(doParallel)


cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
if(is.na(cores)) cores <- detectCores()
# cores <- 6
registerDoParallel(cores)
cat('number of cores using', cores, '. . .\n')

yaml_path = '/ix/djishnu/Hanxi/Dutta/HF/CdVActive/011323/plainER.yaml'
parseRun(yaml_path)
