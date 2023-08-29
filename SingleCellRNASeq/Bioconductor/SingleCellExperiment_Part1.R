library(SingleCellExperiment)
library(scuttle)
library(scran)
box::use(scuttle[sct_logNormCounts = logNormCounts, 
                 sct_addPerCellQC = addPerCellQC,
                 sct_addPerFeatureQC = addPerFeatureQC,
                 ...])
library(scater)
library(uwot)
library(rtracklayer)
box::use(rtracklayer[rt_import = import,...])
library(BiocFileCache)

#####Store primary experimental data#####
# Fill the assays slot
# Download the data
# setwd("~/SingleCellRNASeq/Bioconductor")
bfc <- BiocFileCache("raw_data", ask = F)
calero.counts <- bfcrpath(x = bfc, rnames = "https://www.ebi.ac.uk/biostudies/files/E-MTAB-5522/counts_Calero_20160113.tsv")
# BFC1 
# "raw_data/3f15472e5979_counts_Calero_20160113.tsv" 
mat <- read.delim(calero.counts, header=T, row.names=1, check.names=F)

# Select endogenous genes
spike.mat <- mat[grepl("^ERCC-", rownames(mat)),] 
mat <- mat[grepl("^ENSMUSG", rownames(mat)),] 

# Delete the gene length column
gene.length <- mat[,1]
mat <- as.matrix(mat[,-1]) 
dim(mat)

# Construct a SingleCellExperiment object by providing the counts
sce <- SingleCellExperiment(assays = list(counts = mat))

# Take a peek at various slots in sce
sce

# Access the count data
assay(sce, "counts") # Recommend this general approach
counts(sce)
mat2 <- assay(sce, "counts") 

# Compute a log-transformed normalized expression matrix and store it as another assay
sce <- sct_logNormCounts(sce) # No need to explicitly refer to scuttle::logNormCounts(sce)
sce

# Access the logcount data
assay(sce, "logcounts")
logcounts(sce)
dim(assay(sce, "logcounts"))

# Manually insert the transformed count matrix into the assays
counts_100 <- assay(sce, "counts") + 100
assay(sce, "counts_100") <- counts_100 # Use the  setter function assay() to assign a transformed count matrix to a assays slot
assays(sce)

# Retrieve all the available assays and their names within sce
assays(sce)
assayNames(sce)

# Modify the set of available assays
assays(sce) <- assays(sce)[1:2] # Keep the first two assays
sce

#####Handle cell metadata - colData#####
bfc <- BiocFileCache("raw_data", ask = F)
lun.sdrf <- bfcrpath(x = bfc, rnames = "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5522/E-MTAB-5522.sdrf.txt")
# BFC2 
# "raw_data/3f153d5746eb_E-MTAB-5522.sdrf.txt"
coldata <- read.delim(lun.sdrf, check.names=F) # df
coldata <- coldata[coldata[,"Derived Array Data File"]=="counts_Calero_20160113.tsv",] # Filter out the cells not involved in the count matrix 'mat'
coldata <- DataFrame( # Select interesting columns, and set the library names as the row names
  genotype=coldata[,"Characteristics[genotype]"],
  phenotype=coldata[,"Characteristics[phenotype]"],
  spike_in=coldata[,"Factor Value[spike-in addition]"],
  row.names=coldata[,"Source Name"]
)

coldata

# Populate sce with column data
# Check that the rows of colData refer to the same cells as the columns of count matrix
stopifnot(identical(rownames(coldata), colnames(mat)))

# 1st approach
sce <- SingleCellExperiment(assays = list(counts=mat), colData=coldata)

# 2nd approach
sce <- SingleCellExperiment(list(counts=mat))
colData(sce) <- coldata

# 3rd approach
sce <- SingleCellExperiment(list(counts=mat))
sce$phenotype <- coldata$phenotype
colData(sce)

sce

# Access the column data
colData(sce)

# Extract a single field in the column data
sce$phenotype

# addPerCellQC() function appends a number of quality control metrics to the colData
sce <- sct_addPerCellQC(sce)
colData(sce)

#####Handle feature metadata - rowData/rowRanges#####
# Manually add the gene length for each gene
rowData(sce)$Length <- gene.length

# Access the gene length
rowData(sce)

# Populate sce with rowData
sce <- sct_addPerFeatureQC(sce)
rowData(sce)

# Populate rowRanges slot with genomic coordinates in the form of a GRanges or GRangesList (which stores describes the chromosome, start, and end coordinates of the features (genes, genomic regions))
rowRanges(sce)
# Download the data
# A General Feature Format (GFF) file is a simple tab-delimited text file for describing genomic features. 
# There are several slightly but significantly different GFF file formats. 
# GFF2 file; GFF3 file; GTF file
mm10.gtf <- bfcrpath(bfc, "http://ftp.ensembl.org/pub/release-82/gtf/mus_musculus/Mus_musculus.GRCm38.82.gtf.gz")
gene.data <- rt_import(mm10.gtf)
gene.data <- gene.data[gene.data$type=="gene"]
names(gene.data) <- gene.data$gene_id
mcols(gene.data) # Access the column data using mcols()
colnames(mcols(gene.data)) # Extract the column name
gene.related.cols <- grep("gene_", colnames(mcols(gene.data))) # Extract names of gene-related columns
mcols(gene.data) <- mcols(gene.data)[,gene.related.cols] # Extract gene-related columns
rowRanges(sce) <- gene.data[rownames(sce)]

# Access the first 10 elements in the rowRanges slot
rowRanges(sce)[1:10,]

#####Handle metadata - other data#####
# Populate metadata slot with other data
genes_1 <- c("gene_1", "gene_5")
metadata(sce) <- list(highly_variable_genes = genes_1)
metadata(sce)

genes_2 <- c("gene_6", "gene_9", "gene_20")
metadata(sce)$zero_count_genes <- genes_2
metadata(sce)

#####Subset and combine sce#####
first.10.cells <- sce[,1:10] 
assay(first.10.cells, "counts") |> ncol() # Gene x cell count matrix
colData(first.10.cells) |> nrow() # Cell x phenotype data frame

# Select wild-type cells = slice sce using the colData
wt <- sce[,sce$phenotype=="wild type phenotype"]
assay(wt, "counts") |> ncol() # Gene x cell count matrix
colData(wt) |> nrow() # Cell x phenotype data frame

# Select protein-coding genes = slice sce using the rowData
coding <- sce[rowData(sce)$gene_biotype == "protein_coding",]
assay(coding, "counts") |> nrow() # Gene x cell count matrix
rowData(coding) |> nrow() # Cell x phenotype data frame

# Combine multiple SingleCellExperiment objects
# Assume that all objects involved have the *same row annotation values* and compatible column annotation fields
sce2 <- cbind(sce, sce)
assay(sce2, "counts") |> ncol() 
colData(sce2) |> nrow()

# Assume all objects have the same column annotation values and compatible row annotation fields
sce3 <- rbind(sce, sce)
assay(sce3, "counts") |> nrow() 
rowData(sce3) |> nrow()