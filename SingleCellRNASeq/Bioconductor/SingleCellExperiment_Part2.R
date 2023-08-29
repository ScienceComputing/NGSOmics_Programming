library(SingleCellExperiment)
box::use(scran[scr_computeSumFactors = computeSumFactors,
               scr_clusterCells = clusterCells,
               ...])
box::use(scuttle[sct_logNormCounts = logNormCounts, 
                 sct_addPerCellQC = addPerCellQC,
                 sct_addPerFeatureQC = addPerFeatureQC,
                 ...])
box::use(scater[sc_logNormCounts = logNormCounts, 
                sc_runPCA = runPCA,
                sc_runTSNE = runTSNE,
                sc_runUMAP = runUMAP,
                sc_librarySizeFactors = librarySizeFactors,
                ...])
box::use(uwot[uw_umap = umap,...])
box::use(rtracklayer[rt_import = import,...])

#####Dimensionality reduction - reducedDims#####
# Fill the reducedDims slot
# This slot contains a list of numeric matrices of low-reduced representations of the primary data, 
# where the rows represent the columns of the primary data (i.e., cells), and columns represent the dimensions.
sce <- sc_logNormCounts(sce)
sce <- sc_runPCA(sce)
sce <- sc_runUMAP(sce, n_neighbors = 2)
reducedDim(sce, "PCA") |> dim()
assay(sce, "counts") |> ncol()
sce <- sc_runTSNE(sce, perplexity = 0.1) # Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1)
reducedDim(sce, "TSNE") # Access the TSNE lower embedding results
reducedDims(sce) |> View() # Access all lower embedding results using the accessor reducedDims()

# Manually add content to the reducedDims() slot
u <- uw_umap(X = t(logcounts(sce)), n_neighbors = 2)
reducedDim(sce, "UMAP_uwot") <- u
reducedDims(sce) 

#####Alternative experiments#####
# We have data for a distinct set of features but the same set of samples/cells
# -1 to get rid of the first gene length column.
spike_se <- SummarizedExperiment(list(counts=spike.mat[,-1])) # Remove the first gene length column
spike_se

# Store this SummarizedExperiment in the sce object using the altExp() setter
altExp(sce, "spike") <- spike_se

# Retrieve all of the available alternative Experiments with altExps()
altExp(sce, "spike") <- spike_se
altExps(sce)

# Keep the first 2 cells
sub <- sce[,1:2]
altExp(sub, "spike")

#####Size factors#####
# Get or set a numeric vector of per-cell scaling factors used for normalization
# The first approach - deconvolution-based size factors
sce <- scr_computeSumFactors(sce)
summary(sizeFactors(sce))

# The second approach - library size-derived factors
sizeFactors(sce) <- sc_librarySizeFactors(sce)
summary(sizeFactors(sce))

#####Column labels#####
# Get or set a vector or factor of per-cell labels, e.g., groupings assigned by unsupervised clustering, or predicted cell type identities from classification algorithms
colLabels(sce) <- scr_clusterCells(sce, use.dimred="PCA")
table(colLabels(sce))
# scran::findMarkers()) will automatically retrieve the labels if colLabels(x) is available
