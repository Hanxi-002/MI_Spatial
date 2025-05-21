library(umap)
library(Rtsne)

target_demoData <- Geomx_V3

# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42
# run UMAP
umap_out <-
  umap(t(log2(assayDataElement(target_demoData , elt = "q_norm"))),  
       config = custom_umap)
pData(target_demoData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]

graph_df <- pData(target_demoData)
row.names(graph_df) == row.names(target_demoData@protocolData@data)
graph_df['Status'] <- target_demoData@protocolData@data$Status0p
graph_df['SegmentLabels'] <- target_demoData@protocolData@data$SegmentLabel

ggplot(graph_df,
       aes(x = UMAP1, y = UMAP2, color = Status, shape = SegmentLabels)) +
  geom_point(size = 3) +
  theme_bw()

set.seed(42) # set the seed for tSNE as well
tsne_out <-
  Rtsne(t(log2(assayDataElement(target_demoData , elt = "q_norm"))),
        perplexity = ncol(target_demoData)*.15)
pData(target_demoData)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
ggplot(pData(target_demoData),
       aes(x = tSNE1, y = tSNE2, color = region, shape = class)) +
  geom_point(size = 3) +
  theme_bw()