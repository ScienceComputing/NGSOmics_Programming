#####Count reads into sliding windows across the genome#####
# Goal: quantify the protein binding intensity across the genome
# Notice that BAM files must be sorted by position and have been indexed

frag.len <- 110
window.width <- 10
param <- readParam(minq=20)
# minq - An integer scalar, specifying the minimum mapping quality score for an aligned read.
data <- windowCounts(bam.files=bam.file.path, ext=frag.len, width=window.width, param=param)
# bam.files	- A character vector containing paths to sorted and indexed BAM files. Alternatively, a list of BamFile objects.
# For character inputs, the index files should have the same prefixes as the BAM files. 
# ext	- An integer scalar or a list of two integer scalars/vectors, containing the average length(s) of the sequenced fragments in each library.
# Notice that an ext list should only be specified for datasets that show large differences in the average fragment sizes between libraries.
# width	- An integer scalar specifying the width of the window.
# param: control the read extraction from the BAM files
class(data) # RangedSummarizedExperiment

head(assay(data)) # Extract the count matrix
# Each row maps to a genomic window
# Each column maps to a library

head(rowRanges(data)) # Extract the genomic coordinates of each window

# Extract the total number of reads in each library [library size]
data$totals

#####Filter out low-quality reads#####
param
# A low mapping score suggests that the sequences may be aligned incorrectly or in a way that is not unique.
# The specific threshold value is contingent upon the score range offered by the aligner such as `subread`.

# When multiple reads align to the same genomic location, they can be identified as potential PCR duplicates. 
# We designate these marked reads in the BAM file to be disregarded by configuring the readParam object with dedup=TRUE. 
# This approach helps minimize the variation stemming from uneven amplification across replicates and prevents spurious duplicate-driven differential binding between different groups.
param.new <- readParam(minq=20, dedup=T)
data.new <- windowCounts(bam.files=bam.file.path, ext=frag.len, width=window.width, 
                         param=param.new)
data.new$totals
# Routine differential binding analyses typically advise against removing duplicates 
# since doing so limits the number of reads at each position, which in turn diminishes the detection power of the differential binding, particularly in regions with high abundance.

#####Estimate the fragment length#####
# Visual approach
max.delay <- 500
# dedup.on <- initialize(param, dedup=T) # Equivalent to param.new; remove the duplicates to reduce the size of the artifact spike such that the fragment length peak will be visible as a separate entity.
# A smaller or absent peak indicates the poor immunoprecipitation efficiency.
ext.est <- correlateReads(bam.files=bam.file.path, max.dist=max.delay, param=param.new)
plot(0:max.delay, ext.est, type="l", ylab="CCF", xlab="Delay (bp)")
# correlateReads() - Computes the auto- or cross-correlation coefficients between read positions across a set of delay intervals.
# The peak's location is employed to approximate the fragment length for extending the reads within the context of windowCounts(). 
# Based on the above plot, we derive an estimated length of approximately 110 base pairs.

# More precise approach
maximizeCcf(ext.est)

#####Estimate the fragment length of narrow histone marks#####
n <- 1000

ac.data <- H3K9acData()
ac.ext <- correlateReads(ac.data$Path[1], max.dist=n, param=param.new)

k27.data <- H3K27me3Data()
k27.ext <- correlateReads(k27.data$Path[1], max.dist=n, param=param.new)

k4.data <- H3K4me3Data()
k4.ext <- correlateReads(k4.data$Path[1], max.dist=n, param=param.new)

plot(
  0:n,
  ac.ext,
  col = "blue",
  ylim = c(0, 0.1),
  xlim = c(0, 1000),
  xlab = "Delay (bp)",
  ylab = "CCF",
  pch = 16,
  type = "l",
  lwd = 2
)
lines(0:n,
      k27.ext,
      col = "red",
      pch = 16,
      lwd = 2)
lines(
  0:n,
  k4.ext,
  col = "forestgreen",
  pch = 16,
  lwd = 2
)
legend(
  "topright",
  col = c("blue", "red", "forestgreen"),
  c("H3K9ac", "H3K27me3", "H3K4me3"),
  pch = 16
)

# As shown above, we can use the cross-correlation plots to estimate the fragment length of narrow histone marks such as histone acetylation and H3K4 methylation.
# However, this type of plot fails to estimate the fragment length for regions of diffuse enrichment where bimodality is not obvious (e.g., H3K27 trimethylation).


