library(scRNAseq)
sce <- MacoskoRetinaData()

#####Perform quality control#####
library(scuttle)
is.mito <- grepl("^MT-", rownames(sce))
qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(x=qcstats, percent_subsets="subsets_Mito_percent")
sce <- sce[, !filtered$discard] # Filter out the ineligible cells

#####Normalize the counts#####
sce <- logNormCounts(sce)

#####Select important features#####
library(scran)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(stats=dec, prop=0.1) # Specify the proportion of genes to report as highly variable genes

#####Build the PCA#####
library(scater)
set.seed(1234)
sce <- runPCA(sce, ncomponents=25, subset_row=hvg)

#####Cluster cells#####
library(bluster)
colLabels(sce) <- clusterCells(sce, use.dimred="PCA",
                               BLUSPARAM=NNGraphParam(cluster.fun="louvain")) 

#####Visualize the cell clustering#####
sce <- runUMAP(sce, dimred = "PCA") # String or integer scalar specifying the existing dimensionality reduction results to use
plotUMAP(sce, colour_by="label")

#####Detect differentially expressed marker genes#####
markers <- findMarkers(sce, test.type="wilcox", direction="up", lfc=1)
# Want to find genes that are upregulated in one group compared to another
# Genes with a log-fold change greater than or equal to 1 will be considered as potential markers
