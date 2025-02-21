library(nichenetr)
library(multinichenetr)
library(SingleCellExperiment)
library(dplyr)


######################## create a single cell experiment object ########################
Geomx_V3 <- readRDS("/ix/djishnu/Hanxi/MI_Spatial/Geomx_V3.RDS")
counts <- as.matrix(Geomx_V3@assayData$exprs)
sce <- SingleCellExperiment(assays = list(counts = counts))
# logcounts <- as.matrix(Geomx_V3@assayData$q_norm)
# sce <- SingleCellExperiment(assays = list(logcounts = logcounts))
new_meta <- read.csv('/ix/djishnu/Hanxi/MI_Spatial/CCC/MultiNicheNet/CCC_meta.csv', row.names = 1)
sce@metadata <- new_meta
metadata_as_list <- as.list(new_meta)
sce@metadata <- metadata_as_list

######################## load multinichenet LR database ########################
organism = "human"
options(timeout = 300)

if(organism == "human"){
  lr_network_all = 
    readRDS('/ix/djishnu/Hanxi/MI_Spatial/CCC/MultiNicheNet/lr_network_human_allInfo_30112033.rds') %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS('/ix/djishnu/Hanxi/MI_Spatial/CCC/MultiNicheNet/ligand_target_matrix_nsga2r_final.rds')
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
} else if(organism == "mouse"){
  
  lr_network_all = readRDS(url(
    "https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds"
  )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
  ))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
}


######################## run MultiNicheNet ########################
sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

sce@metadata$sample_id <- row.names(new_meta)
sce@colData <- as(sce@metadata, 'DataFrame')
sample_id = "sample_id"
group_id = "Status"
celltype_id = "new_column"

covariates = NA
batches = NA

# define contrasting groups
contrasts_oi = c("'HF-C','C-HF'")
contrast_tbl = tibble(contrast = c("HF-C","C-HF"), group = c("H","C"))

# define cell types
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi)]

conditions_keep = c("HF", "Control")
sce = sce[, SummarizedExperiment::colData(sce)[,group_id] %in% conditions_keep]


SummarizedExperiment::colData(sce)[,sample_id] = make.names(SummarizedExperiment::colData(sce)$sample_id)
################################ skipping cell filtering ################################

################################ gene filtering ################################
min_cells = 1
min_sample_prop = 0.1
fraction_cutoff = 0.01

frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

abundance_expression_info = process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)





