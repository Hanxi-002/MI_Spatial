## Lavine 2024 DEG with SLIDE LF comparison

### Background
  * Comparing SLIDE LF genes with Lavine 2024 DEGs to see what SLIDE captures

### Input
  * In Mary's directory: Lavine2024 processed scRNA-seq seurat object RDS. 
  * SLIDE LFs: different SLIDE LFs from different SLIDE analyses

### Scripts
  * Lavine_Get_DEGs.R : read in the seurat object and save different DEG analyses to folder as csvs.
  * get_overlap_DEGs_helper.R : helper functions to ouput number of DEGs, number SLIDE LF genes

### Output
  * celltype_DEGs.csv: the DEGs for different cell types in Lavine 2024.
  * 1_vs_all_DEGs.csv: one for each cell type. DEG calculation of 1 vs all for conditions (3 different HF and 1 control)
  * AMI_vs_Control_DEGs.csv: one for each cell type. DEG calculations of Lavine 2024 of AMI vs Donor. 
  
### Notes


