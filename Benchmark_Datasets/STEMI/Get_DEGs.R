library(Seurat)

STEMI <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/Benchmark_Datasets/STEMI/Hanxi-Stemi.rds")

# only extract the first 3 elements of each type (CAD, CON, STE)
conditions <- as.character(sapply(STEMI@meta.data$orig.ident, function(x) substr(x, 1, 3)))
STEMI[['conditions']] <- conditions
Idents(STEMI) <- STEMI@meta.data$conditions


conditions <- unique(STEMI@meta.data$conditions)
deg_results <- list()
# Perform pairwise differential expression analysis
for (i in 1:(length(conditions)-1)) {
  for (j in (i+1):length(conditions)) {
    # Define current pair
    cond1 <- conditions[i]
    cond2 <- conditions[j]
    
    # Find DEGs between current pair
    deg_name <- paste(cond1, "_vs_", cond2, sep="")
    deg_results[[deg_name]] <- FindMarkers(STEMI, ident.1 = cond1, ident.2 = cond2, 
                                           min.pct = 0.1, logfc.threshold = 0.25)
  }
}

################################## compare the SLIDE results with the DEG results ###############################
tested_genes = list()
tested_genes[["murine_BMDM"]] <- c("AASDHPPT", "TVP23B", "SERBP1", "IGGAP1", "CIITA", "GPR107", "DDX5", "MBD2", "AEBP1", "TDRD6", "B3GALNT1")
tested_genes[['LF_genes']] <- c("AASDHPPT", "TVP23B", "SERBP1", "IGGAP1", "GPR107", "DDX5", "TDRD6", "B3GALNT1", "NPIPB6", "HNRNPH3") 
tested_genes[['TF_genes']] <- c("ATF3", "EGR1", "EP300", "KLF3", "MAFK", "MAZ")
tested_genes[['in_vivo']] <- c("B3GALNT1", "TDRD6", "TVP23B", "MAFK", "MAZ")

for (l in tested_genes) {
  cat(l[!(l %in% row.names(STEMI))])
}
# these are the genes not in the dataset
#IGGAP1 AEBP1 TDRD6 B3GALNT1 IGGAP1 TDRD6 B3GALNT1 NPIPB6 ATF3 B3GALNT1 TDRD6

cnt = 0
for (l in tested_genes){
  for (j in deg_results) {
    cat(l[l %in% row.names(j)])
    cnt = cnt+1
    cat(cnt)
  }
}

#EP300 KLF3



