#####Set up the environment#####
BiocManager::install("chipseqDBData")
library(chipseqDBData)
nfya.data <- NFYAData()
# NFYAData() will download five NF-YA (nuclear transcription factor Y subunit alpha) ChIP-seq libraries, with two biological replicates for murine terminal neurons, two replicates for embryonic stem cells and one input ontrol. 
# This uses single-end sequencing data from accession GSE25532 of the NCBI Gene Expression Omnibus.
# All mouse datasets were aligned to the mm10 build of the mouse genome. All BAM files are sorted and indexed, with duplicate reads marked with MarkDuplicates. Index files are named by appending .bai onto the BAM file paths.

nfya.data <- nfya.data[1:4,] # Delete the last row - 5  SRR074401 Input <BamFile>
bam.file.path <- nfya.data$Path

cell.type <- sub("NF-YA ([^ ]+) .*", "\\1", nfya.data$Description)
# ([^ ]+): This part of the pattern captures one or more characters that are not spaces and stores them as a group. 
# [^ ] means any character that is not a space, and + means one or more occurrences.
# .*: This part of the pattern matches any remaining characters in the string.
# \\1 is the replacement pattern. It refers to the first captured group in the regular expression.

design <- model.matrix(~factor(cell.type))
design
# (Intercept) factor(cell.type)TN
# 1           1                   0
# 2           1                   0
# 3           1                   1
# 4           1                   1

colnames(design) <- c("intercept", "cell.type")


#####Load BAM data #####
# BMA: https://en.wikipedia.org/wiki/Binary_Alignment_Map
# BAM is the compressed binary representation of SAM (Sequence Alignment Map), a compact and index-able representation of nucleotide sequence alignments.
# The goal of indexing is to retrieve alignments that overlap a specific location quickly without having to go through all of them. 
# Before indexing, BAM must be sorted by reference ID and then leftmost coordinate.
# BiocManager::install("csaw")

library(csaw)
param <- readParam(minq=20)
# readParam() - Specify read loading parameters
# minq - An integer scalar, specifying the minimum mapping quality score for an aligned read.
data <- windowCounts(bam.file.path, ext=110, width=10, param=param, spacing=50)
# windowCounts - Count the number of extended reads overlapping a sliding window at spaced positions across the genome.
# ext	- An integer scalar or a list of two integer scalars/vectors, containing the average length(s) of the sequenced fragments in each library.
# width	- An integer scalar specifying the width of the window.

#####Filter out uninteresting regions#####
binned <- windowCounts(bam.file.path, bin=T, width=10000, param=param)
# bin	- A logical scalar indicating whether binning should be performed.
keep <- filterWindowsGlobal(data=data, background=binned)$filter > log2(5)
# filter - A numeric vector containing the filter statistic for the given type for each row. The definition of this filter statistic will vary across the different methods.
# > - It checks if the value of filter returned from the result of filterWindowsGlobal() is greater than log2(5).
data <- data[keep,]

#####Calculate normalization factors#####
data <- normFactors(object=binned, se.out=data)
# se.out - A SummarizedExperiment object in which normalization factors are to be stored.

#####Identify differential binding windows#####
library(edgeR)
y <- asDGEList(data)
y <- estimateDisp(y, design)
model.fit <- glmQLFit(y, design, robust=T)
# glmQLFit() - Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
# robust - logical, whether to estimate the prior QL dispersion distribution robustly.
results <- glmQLFTest(model.fit)
head(results$table)
# logFC     logCPM        F      PValue
# 1 0.9234028  0.4971102 3.199060 0.074954084
# 2 1.0987247  1.3174063 5.463417 0.020249631
# 3 1.2665121  1.2833911 7.017433 0.008613499
# 4 0.7936013  0.7804788 2.460198 0.118149332
# 5 0.8706493 -0.3806872 2.002955 0.158301700
# 6 0.7999067 -0.3108220 1.722058 0.191074669

#####Correct for multiple testing#####
merged <- mergeResults(data, results$table, tol=1000L)
merged
