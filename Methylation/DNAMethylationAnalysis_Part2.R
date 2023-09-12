# We now perform the probe-wise analysis, where we test each individual CpG probe for differential methylation in terms of the comparisons of interest (e.g., different cell types in this case), and estimate the p-values and moderated t-statistics for each CpG probe.

#####Probe-wise differential methylation analysis#####
#####Set up comparisons of interest#####
cell.type <- factor(sample.info$Sample_Group) # Factor of our interest
individual <- factor(sample.info$Sample_Source) # Factor that needs to control for, suggested by prior MDS plots

# !Create a design matrix
design <- model.matrix(~0+cell.type+individual, data=sample.info)
colnames(design) <- c(levels(cell.type),levels(individual)[-1])

# Create a contrast matrix for specific comparisons {limma}
cont.mat <- makeContrasts(naive-rTreg,
                          naive-act_naive,
                          rTreg-act_rTreg,
                          act_naive-act_rTreg,
                          levels=design)
cont.mat

#####Build the linear model#####
fit.m <- lmFit(m.val, design)
fit.m.cont <- contrasts.fit(fit.m, cont.mat)
fit.m.cont <- eBayes(fit.m.cont)

#####View the numbers of differentially methylated CpGs at FDR < 0.05#####
summary(decideTests(fit.m.cont))
# naive - rTreg naive - act_naive rTreg - act_rTreg act_naive - act_rTreg
# Down            1621               398                 0                   553
# NotSig        464326            466732            467351                465898
# Up              1404               221                 0                   900

#####Get the table of results for the first contrast (naive - rTreg)#####
colnames(anndata.450k)
anndata.450k.sub <- anndata.450k[match(rownames(m.val),anndata.450k$Name),
                                 c(1:4,12:19,24:ncol(anndata.450k))]
dmps <- topTable(fit.m.cont, num=Inf, coef=1, genelist=anndata.450k.sub)
head(dmps)

# Save the results
write.table(dmps, file="dmps.csv", sep=",", row.names=F)

#####Plot sample-wise beta values for the top differentially methylated CpG sites#####
par(mfrow=c(2,2))
sapply(rownames(dmps)[1:4], function(cpg){
  plotCpg(b.val, cpg=cpg, pheno=sample.info$Sample_Group, ylab = "Beta values")
})
# cpg: A character vector of the CpG position identifiers to be plotted.

#####Differential methylation analysis of regions#####
# Interest: infer the proximal CpGs concordantly differentially methylated, that is, differentially methylated regions
m.val.annotate <- cpg.annotate(object = m.val, datatype = "array", what = "M", 
                               analysis.type = "differential", design = design, 
                               contrasts = T, cont.matrix = cont.mat, 
                               coef = "naive - rTreg", arraytype = "450K")
# cpg.annotate {DMRcate}
str(m.val.annotate)

# Compute a kernel estimate against a null comparison to identify significantly differentially (or variable) methylated regions
dmrs <- dmrcate(m.val.annotate, lambda=1000, C=2)
# C: scaling factor for bandwidth. Gaussian kernel is calculated where lambda/C = sigma. 
# Empirical testing has demonstrated that for both Illumina and bisulfite sequencing data, achieving near-optimal prediction of sequencing-derived DMRs is possible when the value of lambda is set to 1000 and the parameter C is approximately 2. In other words, a Gaussian kernel with a standard deviation of 500 base pairs provides an effective approach.
# Cannot be < 0.2.
# lambda: Gaussian kernel bandwidth for smoothed-function estimation. Also informs DMR bookend definition; gaps >= lambda between significant CpG sites will be in separate DMRs. Support is truncated at 5*lambda. Default is 1000 nucleotides. 
# The values of lambda and C should be chosen with care. For array data, we currently recommend that half a kilobase represent 1 standard deviation of support (lambda=1000 and C=2). 
# If lambda is too small or C too large then the kernel estimator will not have enough support to significantly differentiate the weighted estimate from the null distribution. If lambda is too large then dmrcate will report very long DMRs spanning multiple gene loci, and the large amount of support will likely give Type I errors. If you are concerned about Type I errors we highly recommend using the default value of pcutoff, although this will return no DMRs if no DM CpGs are returned by limma/DSS either.
results.ranges <- extractRanges(dmrs)
results.ranges

#####Set up the grouping variables and colours#####
groups <- pal[1:length(unique(sample.info$Sample_Group))]
names(groups) <- levels(factor(sample.info$Sample_Group))
cols <- groups[as.character(factor(sample.info$Sample_Group))]

#####Draw the plot for the top DMR#####
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, dmr = 2, CpGs = b.val, phen.col = cols, 
         what = "Beta", arraytype = "450K", genome = "hg19")
# dmr: index of ranges (one integer only) indicating which DMR to be plotted
# This plot shows (1) the location of the differentially methylated region in the genome, 
# (2) the position of any genes that are nearby, 
# (3) the base pair positions of the CpG probes, 
# (4) the methylation levels of the individual samples as a heatmap and the mean methylation levels for the various sample groups in the experiment.

