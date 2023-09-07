# For each CpG, we typically have two measurements: a methylated intensity (denoted by M) and an unmethylated intensity (denoted by U).
# These intensity values can be used to determine the proportion of methylation at each CpG locus. 
# Methylation levels are commonly reported as either beta values (β=M/(M+U)) or M-values (Mvalue=log2(M/U)). 
# Beta values are generally preferable for describing the level of methylation at a locus or for graphical presentation because percentage methylation is easily interpretable. 
# M-values are more appropriate for statistical testing, due to their distributional properties.

#####Load methylation specific packages#####
library(methylationArrayAnalysis) 
library(minfi)
library(missMethyl)
library(minfiData)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest) # Illumina manifest is an R object; contains annotation information per CpG probe on the 450k array
library(DMRcate)

#####Load other packages#####
library(knitr)
library(limma)
library(RColorBrewer)
library(Gviz)
library(stringr)

#####Load sample and intensity data#####
# This is a 450k methylation dataset (GSE49667), 
# which contains 10 samples in total: there are 4 different sorted T-cell types (naive, rTreg, act_naive, act_rTreg), collected from 3 different individuals (M28, M29, M30).
data.dir <- system.file("extdata", package = "methylationArrayAnalysis")
list.files(data.dir, recursive = T)

# Return the annotation data
anndata.450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(anndata.450k)

# Import the sample information for the experiment
sample.info <- read.metharray.sheet(data.dir, pattern="SampleSheet.csv") 
sample.info

# Read in the raw data from the IDAT files
rg.set <- read.metharray.exp(targets=sample.info)
# targets: If the targets argument is not NULL it is assumed it has a columned named Basename, and this is assumed to be pointing to the base name of a two color IDAT file, ie. a name that can be made into a real IDAT file by appending either _Red.idat or _Grn.idat.
rg.set # This is an RGChannelSet object that contains raw intensity data from both the red and green colour channels, per sample

sample.info$ID <- paste(sample.info$Sample_Group,sample.info$Sample_Name,sep=".")
sampleNames(rg.set) <- sample.info$ID
rg.set

#####Perform quality control#####
# Estimate a detection p-value for every CpG in every sample
det.p <- detectionP(rg.set) 
head(det.p)
# Calculate detection p-values by comparing the total signal (M+U) for each probe to the background signal level, which is estimated from the negative control probes. 
# Extremely small p-values suggest a strong and reliable signal, whereas larger p-values, such as those greater than 0.01, typically indicate a lower-quality signal

# Examine the mean detection p-values across all samples to identify any samples with low-quality signals
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(det.p), col=pal[factor(sample.info$Sample_Group)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(sample.info$Sample_Group)), fill=pal,
       bg="white")

barplot(colMeans(det.p), col=pal[factor(sample.info$Sample_Group)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(sample.info$Sample_Group)), fill=pal, 
       bg="white")

qcReport(rg.set, sampNames=sample.info$ID, sampGroups=sample.info$Sample_Group, 
         pdf="qcReport.pdf")

# Remove poor quality samples
sample.keep <- colMeans(det.p) < 0.05
rg.set <- rg.set[,sample.keep]
rg.set

# Remove poor quality samples from the sample sheet
dim(sample.info)
sample.info <- sample.info[sample.keep,]
sample.info
dim(sample.info)
colData(rg.set) # Get access to the sample information

# Remove poor quality samples from the detection p-value table
dim(det.p)
det.p <- det.p[,sample.keep]
dim(det.p)

#####Normalization#####
# Goal: weaken the unwanted variation within and between samples
# preprocessFunnorm() is most appropriate for datasets with global methylation differences such as cancer/normal or vastly different tissue types
# preprocessQuantile() is more suited for datasets where you do not expect global differences between your samples, for example a single tissue
norm.set <- preprocessQuantile(rg.set) 
class(norm.set) # GenomicRatioSet

# Create a MethylSet object from the raw data for plotting
raw.set <- preprocessRaw(rg.set) # Convert the Red/Green channel for an Illumina methylation array into methylation signal, without using any normalization
class(raw.set) # MethylSet

# Visualize what the data looks like before and after normalization
par(mfrow=c(1,2))
densityPlot(rg.set, sampGroups=sample.info$Sample_Group,main="Raw", legend=F)
# An **RGChannelSet**, a **MethylSet** or a **matrix**. We either use the getBeta function to get Beta values (for the first two) or we assume the matrix contains Beta values.
legend("topleft", legend = levels(factor(sample.info$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(norm.set), sampGroups=sample.info$Sample_Group,
            main="Normalized", legend=F)
legend("topleft", legend = levels(factor(sample.info$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))


#####Explorative data analysis#####
# Make multi-dimensional scaling (MDS) plots to examine the relationships between the samples in an experiment
# Samples that are more similar to each other cluster together, and samples that are very different are further apart on the plot
# What are the greatest sources of variation are between samples?

# PCA plot with 1st and 2nd PCs
par(mfrow=c(1,2))
plotMDS(getM(norm.set), top=1000, gene.selection="common", 
        col=pal[factor(sample.info$Sample_Group)])
legend("top", 
       legend=levels(factor(sample.info$Sample_Group)), 
       text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(norm.set), top=1000, gene.selection="common",  
        col=pal[factor(sample.info$Sample_Source)]) 
legend("top", 
       legend=levels(factor(sample.info$Sample_Source)), 
       text.col=pal,
       bg="white", cex=0.7)
# This PCA plot shows that the largest source of variation in the methylated intensity data is the difference between individuals (M28, M29, M30) either through the 1st PC or 2nd PC.

# PCA plot with 1st and 3rd PCs (examine higher dimensions to look at other sources of variation)
par(mfrow=c(1,3))
plotMDS(getM(norm.set), top=1000, gene.selection="common", 
        col=pal[factor(sample.info$Sample_Group)], dim=c(1,3))
legend("top", legend=levels(factor(sample.info$Sample_Group)), text.col=pal, 
       cex=0.7, bg="white")

# PCA plot with 2nd and 3rd PCs
plotMDS(getM(norm.set), top=1000, gene.selection="common", 
        col=pal[factor(sample.info$Sample_Group)], dim=c(2,3))
legend("topleft", legend=levels(factor(sample.info$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

# PCA plot with 3rd and 4th PCs
plotMDS(getM(norm.set), top=1000, gene.selection="common", 
        col=pal[factor(sample.info$Sample_Group)], dim=c(3,4))
legend("topright", legend=levels(factor(sample.info$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")
# This PCA plot (3rd PC vs 4th PC) shows that the largest source of variation in the methylated intensity data is the differences between cell types (naive, rTreg, act_naive, act_rTreg).

#####Filter out poor quality CpG probes#####
# Ensure probes are in the same order in the norm.set and det.p objects
det.p.filter <- det.p[match(featureNames(norm.set),rownames(det.p)),] 
dim(det.p)
dim(det.p.filter)

# Remove any CpG probes that have one or more samples showing detection p value >= 0.01
probe.keep <- rowSums(det.p.filter < 0.01) == ncol(norm.set) 
table(probe.keep)
# probe.keep
# FALSE   TRUE 
# 977 484535 

norm.set.filter <- norm.set[probe.keep,]
norm.set.filter

#####Optional section#####
# If the data includes males and females, remove probes on the sex chromosomes
probe.keep <- !(featureNames(norm.set) %in% anndata.450k$Name[anndata.450k$chr %in% c("chrX","chrY")])
table(probe.keep)
norm.set.filter <- norm.set.filter[probe.keep,]

#####Remove all probes affected by SNPs#####
norm.set.filter <- dropLociWithSnps(norm.set.filter)
norm.set.filter

#####Remove cross-reactive probes that map to multiple places in the genome#####
xReactiveProbes <- read.csv(file=paste(data.dir,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=F)
probe.keep <- !(featureNames(norm.set.filter) %in% xReactiveProbes$norm.set.filter)
table(probe.keep)
norm.set.filter <- norm.set.filter[probe.keep,]
norm.set.filter

#####Re-examine the relationship between the samples using MDS plots#####
# PCA plot with 1st and 2nd PCs
par(mfrow=c(1,2))
plotMDS(getM(norm.set.filter), top=1000, gene.selection="common", 
        col=pal[factor(sample.info$Sample_Group)])
legend("top", 
       legend=levels(factor(sample.info$Sample_Group)), 
       text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(norm.set.filter), top=1000, gene.selection="common",  
        col=pal[factor(sample.info$Sample_Source)]) 
legend("top", 
       legend=levels(factor(sample.info$Sample_Source)), 
       text.col=pal,
       bg="white", cex=0.7)
# This PCA plot shows that the largest source of variation in the methylated intensity data is the difference between individuals (M28, M29, M30) only through the 2nd PC.
# In other words, much of the inter-individual variation has been removed in terms of the 1st PC

# PCA plot with 2nd and 3rd PCs (examine higher dimensions to look at other sources of variation)
par(mfrow=c(1,3))
plotMDS(getM(norm.set.filter), top=1000, gene.selection="common", 
        col=pal[factor(sample.info$Sample_Source)], dim=c(1,3))
legend("topleft", legend=levels(factor(sample.info$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")

# PCA plot with 2nd and 3rd PCs
plotMDS(getM(norm.set.filter), top=1000, gene.selection="common", 
        col=pal[factor(sample.info$Sample_Source)], dim=c(2,3))
legend("topleft", legend=levels(factor(sample.info$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")
# This PCA plot shows that the largest source of variation in the methylated intensity data is the difference between individuals (M28, M29, M30) only through the 2nd PC.

# PCA plot with 3rd and 4th PCs
plotMDS(getM(norm.set.filter), top=1000, gene.selection="common", 
        col=pal[factor(sample.info$Sample_Source)], dim=c(3,4))
legend("topright", legend=levels(factor(sample.info$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")

#####Calculate M values and beta values#####
m.val <- getM(norm.set.filter)
head(m.val)

b.val <- getBeta(norm.set.filter)
head(b.val)

#####Visulize the distributions of M values and beta values#####
par(mfrow=c(1,2))
densityPlot(b.val, sampGroups=sample.info$Sample_Group, main="Beta values", 
            legend=F, xlab="Beta values")
legend("top", legend = levels(factor(sample.info$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(m.val, sampGroups=sample.info$Sample_Group, main="M-values", 
            legend=F, xlab="M values")
legend("topleft", legend = levels(factor(sample.info$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))

# Next step, we will perform the probe-wise analysis, where we test each individual CpG probe for differential methylation in terms of the comparisons of interest, and estimate the p-values and moderated t-statistics for each CpG probe.