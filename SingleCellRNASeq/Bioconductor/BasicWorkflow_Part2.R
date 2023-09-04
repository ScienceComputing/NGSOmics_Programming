library(scRNAseq)
sce <- SegerstolpePancreasData()

#####Perform quality control#####
library(scuttle)
qcstats <- perCellQCMetrics(sce)
filtered <- quickPerCellQC(x=qcstats, percent_subsets="altexps_ERCC_percent")
sce <- sce[, !filtered$discard]

#####Normalize the counts#####
sce <- logNormCounts(sce)

#####Select important features#####
library(scran)
dec <- modelGeneVar(sce, block=sce$individual) # A factor specifying the blocking levels for each cell in x. If specified, variance modelling is performed separately within each block and statistics are combined across blocks.
hvg <- getTopHVGs(stats=dec, prop=0.1)

#####Correct batch effects#####
library(batchelor)
set.seed(1234)
sce <- correctExperiments(sce, batch=sce$individual, 
                          subset.row=hvg, correct.all=TRUE)

#####Cluster cells#####
library(bluster)
colLabels(sce) <- clusterCells(sce, use.dimred="corrected") 

#####Visualize the cell clustering#####
sce <- runUMAP(sce, dimred = "corrected") 
gridExtra::grid.arrange(
  plotUMAP(sce, colour_by="label"),
  plotUMAP(sce, colour_by="individual"),
  ncol=2
)

#####Detect differentially expressed marker genes#####
markers <- findMarkers(sce, test.type="wilcox", direction="up", lfc=1)
# Blocking on the individual of origin
