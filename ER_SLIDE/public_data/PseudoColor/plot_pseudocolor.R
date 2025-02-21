library(Seurat)
library(ggplot2)

genes <- c('TVP23B', 'B3GALNT1', 'TDRD6', 'MAFK', 'MAZ')

################################################################################################################
#---------------------------------------------------------
#pseudo color plot for macrophages in GeoMX
#---------------------------------------------------------
umap <- read.csv('/ix/djishnu/Hanxi/MI_Spatial/Cell_Oracle/CD68/020124/CD68_umap_coordinates.csv', row.names = 1)
umap <- as.data.frame(umap)

geomx <- readRDS('/ix/djishnu/Hanxi/MI_Spatial/Geomx_V3.RDS')
geomx <- subset(geomx, select = geomx@phenoData@data$`CD 68` == TRUE)
dim(geomx)
dim(umap)

for (gene in genes){
  ggplot(data=umap, aes(x=UMAP1, y=UMAP2, 
                        color=geomx@assayData$q_norm[gene,])) + 
    geom_point(aes(shape=condition, size=condition), stroke=1.2) +  # Add size aesthetic
    scale_color_gradient(name=gene, low="lightgrey", high="purple") +
    scale_shape_manual(name="Condition", 
                       values=c("HF"=18, "Control"=19)) +
    scale_size_manual(name="Condition",
                      values=c("HF"=2, "Control"=1)) +  # Separate scale for size
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          legend.key.size = unit(1, "cm"))  # Larger legend keys
  ggsave(paste0("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/PseudoColor/Geomx_", gene, '.pdf'), width=8, height=6)
}
################################################################################################################
#---------------------------------------------------------
#Pseudo plot for Lavine Macrophages. 
#---------------------------------------------------------
Lavine <- readRDS('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Macrophage_Seurat_obj.RDS')
Lavine <- subset(x = Lavine, subset =HF.etiology %in% c('Donor', 'AMI'))
umap <- as.data.frame(Lavine@reductions$umap@cell.embeddings)
colnames(umap) <- c('UMAP1', 'UMAP2')
# check the order for merging
#sum(names(Lavine$HF.etiology) != row.names(umap))

umap['condition'] = Lavine$HF.etiology
for (gene in genes){
  ggplot(data=umap, aes(x=UMAP1, y=UMAP2, 
                        color=Lavine@assays$RNA$data[gene, ])) + 
    geom_point(aes(shape=condition, size=condition), stroke=1.2) +  # Add size aesthetic
    scale_color_gradient(name=gene, low="lightgrey", high="purple") +
    scale_shape_manual(name="Condition", 
                       values=c("AMI"=18, "Donor"=19)) +
    scale_size_manual(name="Condition",
                      values=c("AMI"=2, "Donor"=1)) +  # Separate scale for size
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          legend.key.size = unit(1, "cm"))  # Larger legend keys
  
  ggsave(paste0("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/PseudoColor/Lavine_", gene, '.pdf'), width=8, height=6)
}


################################################################################################################
#---------------------------------------------------------
#re run Umap in Macrophages after subsetting Donor and AMI 
#---------------------------------------------------------
Lavine <- readRDS('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Lavine/Myeloid/Macrophage_Seurat_obj.RDS')
Donor_AMI_Mac <- subset(x = Lavine, subset =HF.etiology %in% c('Donor', 'AMI'))

# Calculate UMAP from normalized Seurat object
#Donor_AMI_Mac <- FindVariableFeatures(Donor_AMI_Mac, selection.method = "vst", nfeatures = 2000)
#Donor_AMI_Mac <- ScaleData(Donor_AMI_Mac)
Donor_AMI_Mac <- RunPCA(Donor_AMI_Mac, features = VariableFeatures(object = Donor_AMI_Mac))
Donor_AMI_Mac <- RunUMAP(Donor_AMI_Mac, dims = 1:20)

# Basic UMAP plot
p <- DimPlot(Donor_AMI_Mac, reduction = "umap", group.by = "orig.ident") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"))
p
# If you want different shapes by condition (assuming you have a "condition" column)
p <- DimPlot(Donor_AMI_Mac, 
             reduction = "umap", 
             group.by = "HF.etiology",
             pt.size = 1,) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"))

print(p)


umap <- as.data.frame(Donor_AMI_Mac@reductions$umap@cell.embeddings)
colnames(umap) <- c('UMAP1', 'UMAP2')
umap['condition'] = Donor_AMI_Mac$HF.etiology

for (gene in genes){
  ggplot(data=umap, aes(x=UMAP1, y=UMAP2, 
                        color=Donor_AMI_Mac@assays$RNA$data[gene, ])) + 
    geom_point(aes(shape=condition, size=condition), stroke=1.2) +  # Add size aesthetic
    scale_color_gradient(name=gene, low="lightgrey", high="purple") +
    scale_shape_manual(name="Condition", 
                       values=c("AMI"=18, "Donor"=19)) +
    scale_size_manual(name="Condition",
                      values=c("AMI"=2, "Donor"=1)) +  # Separate scale for size
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          legend.key.size = unit(1, "cm"))  # Larger legend keys
  ggsave(paste0("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/PseudoColor/NewUmap_Lavine_", gene, '.pdf'), width=8, height=6)
}
################################################################################################################
# --------------------------------
# pseudo color plot for Kramnn 10X Visium, Macrophages
# --------------------------------
Kramann <- readRDS('/ix/djishnu/Mary/MI_data/Kramann_2023_Visium/data/Kramann_macrophages_PP.rds')
umap <- as.data.frame(Kramann@reductions$umap@cell.embeddings) 
colnames(umap) <- c('UMAP1', 'UMAP2')

umap['condition'] = Kramann$major_labl
for (gene in genes){
  ggplot(data=umap, aes(x=UMAP1, y=UMAP2, 
                        color=Kramann@assays$SCT$data[gene, ])) + 
    geom_point(aes(shape=condition, size=condition), stroke=1.2) +  # Add size aesthetic
    scale_color_gradient(name=gene, low="lightgrey", high="purple") +
    scale_shape_manual(name="Condition", 
                       values=c("FZ"=18, "IZ"=19)) +
    scale_size_manual(name="Condition",
                      values=c("FZ"=2, "IZ"=1)) +  # Separate scale for size
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          legend.key.size = unit(1, "cm"))  # Larger legend keys
  ggsave(paste0("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/PseudoColor/Kramann_", gene, '.pdf'), width=8, height=6)
}

################################################################################################################
# --------------------------------
# re-run umap for Kramann Mcrophages.
# --------------------------------
Kramann <- readRDS('/ix/djishnu/Mary/MI_data/Kramann_2023_Visium/data/Kramann_macrophages_PP.rds')
Kramann <- ScaleData(Kramann)
Kramann <- RunUMAP(Kramann, dims = 1:20)

# Basic UMAP plot
p <- DimPlot(Kramann, reduction = "umap", group.by = "major_labl") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"))
p

p <- DimPlot(Kramann, 
             reduction = "umap", 
             group.by = "major_labl",
             pt.size = 1,) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"))

print(p)

for (gene in genes){
  ggplot(data=umap, aes(x=UMAP1, y=UMAP2, 
                        color=Kramann@assays$SCT$data[gene, ])) + 
    geom_point(aes(shape=condition, size=condition), stroke=1.2) +  # Add size aesthetic
    scale_color_gradient(name=gene, low="lightgrey", high="purple") +
    scale_shape_manual(name="Condition", 
                       values=c("FZ"=18, "IZ"=19)) +
    scale_size_manual(name="Condition",
                      values=c("FZ"=2, "IZ"=1)) +  # Separate scale for size
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          legend.key.size = unit(1, "cm"))  # Larger legend keys
  ggsave(paste0("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/PseudoColor/NewUmap_Kramann_", gene, '.pdf'), width=8, height=6)
}
