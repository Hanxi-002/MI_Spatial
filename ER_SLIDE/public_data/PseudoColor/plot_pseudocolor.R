library(Seurat)
library(ggplot2)
library(dplyr)

genes <- c('TVP23B', 'B3GALNT1', 'TDRD6', 'MAFK', 'MAZ')


######################################### plot in the umap space #########################################
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
################################################ Plot as ################################################################
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




############################################### Plot as Bar Plot ###############################################
################################################################################################################
#---------------------------------------------------------
#pseudo color plot for macrophages in GeoMX
#---------------------------------------------------------
source('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/Calc_Z_Hat_Helper.R')
genes <- c('TVP23B', 'B3GALNT1', 'TDRD6', 'MAFK', 'MAZ')


geomx <- readRDS('/ix/djishnu/Hanxi/MI_Spatial/Geomx_V3.RDS')
geomx <- subset(geomx, select = geomx@phenoData@data$`CD 68` == TRUE)
# the label for distance
dist_label <- read.csv('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/Within_Region/121223_hf_mac/Data/y.csv', row.names = 1)

# create data for plotting
#data = as.data.frame(geomx@assayData$q_norm[gene,])
data = as.data.frame(geomx@protocolData@data$Status)
colnames(data) = 'condition'
row.names(data) = row.names(geomx@protocolData@data)

# not all HF sampels got a distance label
plot_df <- data %>% filter( condition == "Control" | (condition == "HF" & row.names(data) %in% row.names(dist_label)))

for (i in 1:nrow(plot_df)){
  if (plot_df$condition[i] == 'HF'){
    label_idx <- which(row.names(dist_label) == row.names(plot_df)[i])
    dist <- dist_label[label_idx, 1]
    if (dist == 0){
      dist = 'close to RF'
    } else if (dist == 1){
      dist = 'close to AF'
    }
    plot_df$condition[i] <- dist
  }
}

# subset geomx to keep all the samples in the plot_df
geomx <- geomx[, colnames(geomx) %in% rownames(plot_df)]



for (gene in genes) {
  plot_df['expr'] = geomx@assayData$q_norm[gene,]
  corr_df <- plot_df
  corr_df$condition <- ifelse(corr_df$condition == "Control", 0, 1)
  cat('Correlation for gene', gene, ' is ', cor(corr_df$condition, corr_df$expr))
  
  
  # Mann-Whitney U Test
  non_controls = plot_df[plot_df$condition %in% c("close to RF", "close to AF"), ]
  p_vals = perform_mw_tests(non_controls, condition_col = "condition", score_col ="expr")
  
  custom_order <- c("close to RF", "close to AF")
  # Convert the condition column to a factor with the custom order
  non_controls$condition <- factor(non_controls$condition, levels = custom_order)
  # Then use your original plot code
  ggplot(non_controls, aes(x = condition, y = expr)) +
    geom_boxplot(width = 0.7, outlier.shape = 16, outlier.size = 2) +
    # Add individual data points with jitter for better visualization
    # geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    theme_classic() +
    labs(
      title = paste0(gene, " Expression by Condition"),
      x = "Condition",
      y = "Expression"
    ) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
      axis.title = element_text(size = 14)
    )
  ggsave(paste0("/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/PseudoColor/GeoMX/", gene, '_box_plot.pdf'), width=8, height=6)
  write.csv(non_controls, paste0('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/PseudoColor/GeoMX/', gene, '_box_plot.csv'))
  write.csv(p_vals, paste0('/ix/djishnu/Hanxi/MI_Spatial/ER_SLIDE/public_data/PseudoColor/GeoMX/', gene, '_p_vals.csv'))
}

##########################################
# mac_geomx <- readRDS('/ix/djishnu/Hanxi/MI_Spatial/Geomx_V3.RDS')
# mac_geomx <- subset(mac_geomx, select = mac_geomx@phenoData@data$`CD 68` == TRUE)
# mac_y = as.data.frame(mac_geomx@protocolData@data$Status)
# row.names(mac_y) <- row.names(mac_geomx@protocolData@data)
# 
# expr <- mac_geomx@assayData$q_norm['TVP23B', ]
# mac_y[, 1] <- ifelse(mac_y[, 1] == "HF", 1, 0)
# 
# cor(expr, mac_y[, 1]) # -0.3500334

