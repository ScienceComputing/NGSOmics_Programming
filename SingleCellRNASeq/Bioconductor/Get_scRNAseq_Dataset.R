#####Read the count in csv into R (in-memory representations)#####
# External data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241
library(BiocFileCache)
tools::R_user_dir("BiocFileCache", "cache")
# [1] "/Users/anniliu/Library/Caches/org.R-project.R/R/BiocFileCache"
bfc <- BiocFileCache(cache = getBFCOption("CACHE"), ask=F) # logical(1) Ask before creating, updating, overwriting, or removing cache or local file locations
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85241/suppl/GSE85241%5Fcellsystems%5Fdataset%5F4donors%5Fupdated%2Ecsv%2Egz"

# Make a symbolic link 
# Establish the file path location to load
muraro.fname <- bfcrpath(x = bfc, rnames = url) # character(1) Name of object in file cache
# BFC2 
# "/Users/anniliu/Library/Caches/org.R-project.R/R/BiocFileCache/1b294171ecdb_GSE85241%255Fcellsystems%255Fdataset%255F4donors%255Fupdated%252Ecsv%252Egz" 
local.name <- URLdecode(basename(url))
# [1] "GSE85241_cellsystems_dataset_4donors_updated.csv.gz"
unlink(local.name)
if (.Platform$OS.type=="windows") {
  file.copy(from = muraro.fname, to = local.name)
} else {
  file.symlink(from = muraro.fname, to = local.name)
}

mat <- as.matrix(read.delim(local.name))
dim(mat)
ob.s <- object.size(mat)
print(ob.s , units = "Mb", standard = "legacy")
# 450.1 Mb

# A memory-save method: read the count in sparse format
library(scuttle)
sparse.mat <- readSparseCounts(local.name)
dim(sparse.mat)
ob.s <- object.size(sparse.mat)
print(ob.s , units = "Mb", standard = "legacy")
# 144 Mb

#####Read the count in xls into R (in-memory representations)#####
# External data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61533
bfc <- BiocFileCache(cache = "raw_data", ask=F)
bfc@cache # Cutomize the directory that stores the cache files
# [1] "raw_data"
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE61nnn/GSE61533/suppl/GSE61533%5FHTSEQ%5Fcount%5Fresults%2Exls%2Egz"
wilson.fname <- bfcrpath(x = bfc, 
                         rnames = url)

library(R.utils)
wilson.name2 <- "./derive_data/GSE61533_HTSEQ_count_results.xls"
gunzip(filename=wilson.fname, destname=wilson.name2, remove=F, overwrite=T)

library(readxl)
library(tibble)
all.counts <- read_excel("./derive_data/GSE61533_HTSEQ_count_results.xls")
all.count <- column_to_rownames(.data = all.counts, var = "ID")
dim(all.counts)

#####Read the count from Cellranger output into R (in-memory representations)#####
library(DropletTestFiles)
listTestFiles() |> View()
listTestFiles()@listData[["file.dataset"]] # Find the dataset
cached <- getTestFile(path = "tenx-2.1.0-pbmc4k/1.0.0/filtered.tar.gz")
#                                                                            EH3770 
"/Users/anniliu/Library/Caches/org.R-project.R/R/ExperimentHub/2f833fdd2870_3806" 
fpath <- "tenx-2.1.0-pbmc4k"
untar(tarfile=cached, exdir=fpath) # exdir: the directory to extract files to (the equivalent of tar -C). It will be created if necessary
library(DropletUtils)
sce <- read10xCounts(samples = "tenx-2.1.0-pbmc4k/filtered_gene_bc_matrices/GRCh38") # samples: a character vector containing one or more directory names, each corresponding to a 10X sample. Each directory should contain a matrix file, a gene/feature annotation file, and a barcode annotation file.
sce

cached2 <- getTestFile(path = "tenx-3.1.0-5k_pbmc_protein_v3/1.0.0/filtered.tar.gz")
fpath2 <- "tenx-3.1.0-5k_pbmc_protein_v3"
untar(tarfile=cached2, exdir=fpath2)
# gunzip(filename=paste0(fpath2,"/filtered_feature_bc_matrix/", "barcodes.tsv.gz"), 
#        destname=paste0(fpath2,"/filtered_feature_bc_matrix/", "barcodes.tsv"), remove=F, overwrite=T)
# gunzip(filename=paste0(fpath2,"/filtered_feature_bc_matrix/", "features.tsv.gz"), 
#        destname=paste0(fpath2,"/filtered_feature_bc_matrix/", "features.tsv"), remove=F, overwrite=T)
# gunzip(filename=paste0(fpath2,"/filtered_feature_bc_matrix/", "matrix.mtx.gz"), 
#        destname=paste0(fpath2,"/filtered_feature_bc_matrix/", "matrix.mtx"), remove=F, overwrite=T)

dirA <- "tenx-2.1.0-pbmc4k/filtered_gene_bc_matrices/GRCh38"
dirB <- "tenx-3.1.0-5k_pbmc_protein_v3/filtered_feature_bc_matrix"
sce2 <- read10xCounts(c(dirA, dirB))
# gene information differs between runs

#####Read the count in HDF5 into R (in-memory representations)#####
library(zellkonverter)
demo <- system.file("extdata", "krumsiek11.h5ad", package = "zellkonverter")
# [1] "/Users/anniliu/Library/R/x86_64/4.3/library/zellkonverter/extdata/krumsiek11.h5ad"
sce3 <- readH5AD(demo)
sce3

#####Read the count in loom into R (in-memory representations)#####
library(LoomExperiment)
demo <- system.file("extdata", "L1_DRG_20_example.loom", package = "LoomExperiment")
# [1] "/Users/anniliu/Library/R/x86_64/4.3/library/LoomExperiment/extdata/L1_DRG_20_example.loom"
scle <- import(con=demo, type="SingleCellLoomExperiment")
scle

#####Read the large scRNA-seq data into R (out-of-memory representations)#####
library(TENxBrainData)
sce.brain <- TENxBrainData20k() # The TENxBrainData will return a SingleCellExperiment object containing the full data set, i.e., 1306127 cells. The TENxBrainData20k will return a subset of 20,000 cells from this full data set, as described on the 10X Genomics website. The latter is often useful for quickly testing scripts prior to running them on the full data set.
sce.brain
counts(sce.brain) # This matrix object points to the much larger HDF5 file that actually contains the data
object.size(counts(sce.brain))
file.info(path(counts(sce.brain)))$size

