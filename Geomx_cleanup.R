library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(knitr)
library(dplyr)
library(ggforce)
library(ggplot2)

#https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html

#http://www.bioconductor.org/packages/release/bioc/vignettes/GeomxTools/inst/doc/Developer_Introduction_to_the_NanoStringGeoMxSet.html

#https://www.bioconductor.org/packages/release/bioc/manuals/GeomxTools/man/GeomxTools.pdf
##################################################################
##                          Load Files                          ##
##################################################################

datadir <- file.path("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Dutta_NanoString/")

#DCCs files - expression count data and sequencing quality metadata
DCCFiles <- dir(file.path(datadir, "DCC-20230127"), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)

#PKCs file(s) - probe assay metadata describing the gene targets present in the data
PKCFiles <- dir(file.path(datadir), pattern = ".pkc$", full.names = TRUE, recursive = TRUE)

#SampleAnnotationFile <- dir(file.path(datadir, "annotation"), pattern = ".xlsx$", full.names = TRUE, recursive = TRUE)
SampleAnnotationFile <- "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/raw_data/Dcc_Initial_Dataset_Dutta02.xlsx"

demoData <- readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                   pkcFiles = PKCFiles,
                                   phenoDataFile = SampleAnnotationFile,
                                   phenoDataSheet = "Sheet1",
                                   phenoDataDccColName = "DccNames",
                                   protocolDataColNames = c("ScanLabel", "SegmentLabel", "ROILabel", 'Status'),
                                   experimentDataColNames = c("Status"))
pkcs <- annotation(demoData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))
##################################################################
##                 assess quality for every ROI                 ##
##################################################################
demoData <- shiftCountsOne(demoData, useDALogic=TRUE)
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 5,   # Minimum negative control counts (10)
       maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100)
       minArea = 1000)
demoData <- setSegmentQCFlags(demoData, qcCutoffs=QC_params)
QCResults <- protocolData(demoData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

kable(QC_Summary, caption = "QC Summary Table for each Segment")

warning_segments <- demoData@phenoData@data[which(QCResults$QCStatus == "WARNING"), ]
status <- demoData@protocolData@data$Status[which(QCResults$QCStatus == "WARNING")]
#warning_segments <- cbind(warning_segments, tmp)
##################################################################
##                         set probe QC                         ##
##################################################################

demoData <- setBioProbeQCFlags(demoData,
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20),
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(demoData)[["QCFlags"]]


# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))

#################################################################
##                    Region to Gene Matrix                    ##
#################################################################
length(unique(featureData(demoData)[["TargetName"]]))
target_demoData <- aggregateCounts(demoData)
dim(target_demoData)
#exprs(target_demoData)[1:5, 1:2]

##################################################################
##                   Limit of Quantifications                   ##
##################################################################
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_demoData))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_demoData)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_demoData)[, vars[1]] *
             pData(target_demoData)[, vars[2]] ^ cutoff)
  }
}

pData(target_demoData)$LOQ <- LOQ

LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_demoData)$Module == module
  Mat_i <- t(esApply(target_demoData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_demoData)$TargetName, ]

pData(target_demoData)$GenesDetected <-
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_demoData)$GeneDetectionRate <-
  pData(target_demoData)$GenesDetected / nrow(target_demoData)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_demoData)$DetectionThreshold <-
  cut(pData(target_demoData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

##################################################################
##                        Gene Filtering                        ##
##################################################################
library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_demoData)]
fData(target_demoData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_demoData)$DetectionRate <-
  fData(target_demoData)$DetectedSegments / nrow(pData(target_demoData))

# Gene of interest detection table
goi <- c("PDCD1", "CD274", "IFNG", "CD8A", "CD68", "EPCAM",
         "KRT18", "NPHS1", "NPHS2", "CALB1", "CLDN8")
goi_df <- data.frame(
  Gene = goi,
  Number = fData(target_demoData)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(target_demoData)[goi, "DetectionRate"]))


# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 3, 5, 10, 20, 30))

# plot_detect$Number <-
#   unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
#                 function(x) {sum(fData(target_demoData)$DetectionRate >= x)}))

plot_detect$Number <-
  unlist(lapply(c(0.01, 0.03, 0.05, 0.1, 0.2, 0.3),
                function(x) {sum(fData(target_demoData)$DetectionRate >= x)}))


plot_detect$Rate <- plot_detect$Number / nrow(fData(target_demoData))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

negativeProbefData <- subset(fData(target_demoData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)

#5% detection rate
target_demoData <-
  target_demoData[fData(target_demoData)$DetectionRate >= 0.03 |
                    fData(target_demoData)$TargetName %in% neg_probes, ]
dim(target_demoData)

##################################################################
##                    Quantile Normalization                    ##
##################################################################
target_demoData <- normalize(target_demoData , norm_method="quant",
                             desiredQuantile = .75, toElt = "q_norm")

#saveRDS(target_demoData, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Geomx_v3.RDS")

##################################################################
##                        Zero Filtering                        ##
##################################################################
data_mat <- as.matrix(exprs(target_demoData))
dim(data_mat)
y <- recode(target_demoData@protocolData@data[["Status"]], 'Control' = 0, 'HF' = 1)
write.csv(t(data_mat), '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/AllCell/022423/Data/x.csv')
write.csv(y, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/AllCell/022423/Data/y.csv')


subset_demoData <- subset(target_demoData, select = phenoData(target_demoData)[["Segment (Name/ Label)"]] == 'CD 68')
data_mat <- as.matrix(exprs(subset_demoData))
y <- recode(subset_demoData@protocolData@data[["Status"]], 'Control' = 0, 'HF' = 1)
if(dim(data_mat)[2] != length(y)) {stop("length doesn't match between x and y...")}
write.csv(t(data_mat), '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/CD68/022423/Data/x.csv')
write.csv(y, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/CD68/022423/Data/y.csv')

subset_demoData <- subset(target_demoData, select = phenoData(target_demoData)[["Segment (Name/ Label)"]] == 'active fibroblast')
data_mat <- as.matrix(exprs(subset_demoData))
y <- recode(subset_demoData@protocolData@data[["Status"]], 'Control' = 0, 'HF' = 1)
if(dim(data_mat)[2] != length(y)) {stop("length doesn't match between x and y...")}
write.csv(t(data_mat), '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/ActiveFibro/022423/Data/x.csv')
write.csv(y, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/ActiveFibro/022423/Data/y.csv')

subset_demoData <- subset(target_demoData, select = phenoData(target_demoData)[["Segment (Name/ Label)"]] == 'resting fibroblast')
data_mat <- as.matrix(exprs(subset_demoData))
y <- recode(subset_demoData@protocolData@data[["Status"]], 'Control' = 0, 'HF' = 1)
if(dim(data_mat)[2] != length(y)) {stop("length doesn't match between x and y...")}
write.csv(t(data_mat), '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/RestingFibro/022423/Data/x.csv')
write.csv(y, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/RestingFibro/022423/Data/y.csv')


##################################################################
##                     cell type comparison                     ##
##################################################################
#data_mat <- as.matrix(exprs(target_demoData))
#dim(data_mat)
subset_demoData <- subset(target_demoData, select = phenoData(target_demoData)[["Segment (Name/ Label)"]] != 'CD 68')
subset_demoData <- subset_demoData[ , subset_demoData@protocolData@data[['Status']] == "Control"]
dim(subset_demoData)
data_mat <- as.matrix(exprs(subset_demoData))
dim(data_mat)
y <- recode(subset_demoData@phenoData@data[["Segment (Name/ Label)"]], 'resting fibroblast' = 0, 'active fibroblast' = 1)

write.csv(t(data_mat), '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/HF/RestVActive/022823/Data/x.csv')
write.csv(y, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/HF/RestVActive/022823/Data/y.csv')

subset_demoData <- subset(target_demoData, select = phenoData(target_demoData)[["Segment (Name/ Label)"]] != 'resting fibroblast')
subset_demoData <- subset_demoData[ , subset_demoData@protocolData@data[['Status']] == "Control"]
dim(subset_demoData)
data_mat <- as.matrix(exprs(subset_demoData))
dim(data_mat)
y <- recode(subset_demoData@phenoData@data[["Segment (Name/ Label)"]], 'CD 68' = 0, 'active fibroblast' = 1)

write.csv(t(data_mat), '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/HF/CdVActive/022823/Data/x.csv')
write.csv(y, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/HF/CdVActive/022823/Data/y.csv')


subset_demoData <- subset(target_demoData, select = phenoData(target_demoData)[["Segment (Name/ Label)"]] != 'active fibroblast')
subset_demoData <- subset_demoData[ , subset_demoData@protocolData@data[['Status']] == "Control"]
dim(subset_demoData)
data_mat <- as.matrix(exprs(subset_demoData))
dim(data_mat)
y <- recode(subset_demoData@phenoData@data[["Segment (Name/ Label)"]], 'CD 68' = 1, 'resting fibroblast' = 0)

write.csv(t(data_mat), '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/HF/CdVRest/022823/Data/x.csv')
write.csv(y, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MultiOmic/Dutta_Spatial/ER_SLIDE/HF/CdVRest/022823/Data/y.csv')

##################################################################
##                   Within Region Comparison                   ##
##################################################################
y <- as.matrix(read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/052223/Data/y.csv", row.names = 1))

GetMatrix <- function(status, data_object){
  dim(data_object)
  subset_demoData <- data_object[ , data_object@protocolData@data[['Status']] == status]
  data_mat <- as.matrix(exprs(subset_demoData))
  names <- paste0(row.names(y), ".dcc")
  final_x <- t(data_mat[ ,which(colnames(data_mat) %in% names)])
  cat(dim(final_x))
  return(final_x)
}

CD68_control <- GetMatrix("Control", target_demoData)
control_y <- y[which(row.names(y) %in% colnames(CD68_control))  , ]
CD68_HF <- GetMatrix("HF", target_demoData)

write.csv(t(CD68), "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/MI_Spatial/ER_SLIDE/Within_Region/052223_all/Data/x.csv")



