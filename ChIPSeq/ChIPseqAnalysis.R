library(chipseq) 
library(GenomicFeatures) 
library(lattice)

#####Load data#####
data(cstest)
cstest
class(cstest) # CompressedGRangesList
cstest$ctcf |> class() # GenomicRanges

#####Estimate fragment length#####
frag.len <- estimate.mean.fraglen(cstest$ctcf, method="correlation")
frag.len[!is.na(frag.len)] 
# chr10 chr11 chr12 
# 340   340   340 

#####Extend all reads for one lane#####
ctcf.long <- resize(cstest$ctcf, width = 200)
ctcf.long

#####Report coverage for one lane#####
cov.ctcf <- coverage(ctcf.long)
cov.ctcf

#####Report island for one lane#####
islands <- slice(cov.ctcf, lower = 1)
islands

#####Report number of reads in the island for one lane#####
viewSums(islands)
nread.tab <- table(viewSums(islands) / 200)
nread.tab

#####Report maximum coverage depth within the island for one lane#####
viewMaxs(islands)
depth.tab <- table(viewMaxs(islands))
depth.tab

#####Report number of reads in the island for multiple lanes#####
islandReadSummary <- function(x) {
  r <- resize(x, 200)
  s <- slice(coverage(r), lower = 1)
  tab <- table(viewSums(s) / 200)
  df <- DataFrame(tab)
  colnames(df) <- c("chromosome", "nread", "count")
  df$nread <- as.integer(df$nread)
  df
}
# islandReadSummary(cstest$ctcf) one lane

#####Report number of reads in the island for multiple lanes#####
nread.islands <- DataFrameList(lapply(cstest, islandReadSummary)) 
nread.islands$ctcf
nread.islands$gfp
nread.islands <- stack(nread.islands, "sample") # Stack nreads of ctcf and nreads of gfp into one big data frame
nread.islands 

xyplot(log(count) ~  nread | sample, as.data.frame(nread.islands),
       subset = (chromosome == "chr10" & nread <= 40), 
       layout = c(1, 2), pch = 16, type = c("p", "g"))

xyplot(log(count) ~ nread | sample, data = as.data.frame(nread.islands), 
       subset = (chromosome == "chr10" & nread <= 40), 
       layout = c(1, 2), pch = 16, type = c("p", "g"), 
       panel = function(x, y, ...) {
         panel.lmline(x[1:2], y[1:2], col = "black")
         panel.xyplot(x, y, ...)
       })

#####Report maximum coverage depth within the island for multiple lanes#####
islandDepthSummary <- function(x) 
{
  r <- resize(x, 200) 
  s <- slice(coverage(r), lower = 1) 
  tab <- table(viewMaxs(s) / 200) 
  df <- DataFrame(tab) 
  colnames(df) <- c("chromosome", "depth", "count")
  df$depth <- as.integer(df$depth) 
  df
} 
depth.islands <- DataFrameList(lapply(cstest, islandDepthSummary))
depth.islands <- stack(depth.islands, "sample")

#####Visualize maximum coverage depth within the island for multiple lanes#####
xyplot(
  log(count) ~ depth | sample,
  as.data.frame(depth.islands),
  subset = (chromosome == "chr10" & depth <= 20),
  layout = c(1, 2),
  pch = 16,
  type = c("p", "g"),
  panel = function(x, y, ...) {
    lambda <- 2 * exp(y[2]) / exp(y[1])
    null.est <- function(xx) {
      xx * log(lambda) - lambda - lgamma(xx + 1)
    }
    log.N.hat <- null.est(1) - y[1]
    panel.lines(1:10,-log.N.hat + null.est(1:10), col = "black")
    panel.xyplot(x, y, ...)
  }
)

#####Visualize maximum coverage depth within the island for one lane#####
islandDepthPlot(cov.ctcf)

#####Estimate the peak cutoff for a specific FDR#####
peakCutoff(cov.ctcf, fdr = 0.0001) 
# [1] 6.959837

#####Obtain putative binding sites#####
peaks.ctcf <- slice(cov.ctcf, lower = 8) 
peaks.ctcf 

#####Summarize peaks#####
peaks <- peakSummary(peaks.ctcf) 

#####Compute strand-specific coverage#####
peak.depths <- viewMaxs(peaks.ctcf)
cov.pos <- coverage(ctcf.long[strand(ctcf.long) == "+"]) 
cov.neg <- coverage(ctcf.long[strand(ctcf.long) == "-"])
peaks.pos <- Views(cov.pos, ranges(peaks.ctcf)) 
peaks.neg <- Views(cov.neg, ranges(peaks.ctcf))

wpeaks <- tail(order(peak.depths$chr10), 4)
wpeaks
# [1]  971  989 1079  922
coverageplot(peaks.pos$chr10[wpeaks[1]], peaks.neg$chr10[wpeaks[1]])

#####Visualize four highest peaks on chromosome 10#####
lapply(1:4, function (i) {
  coverageplot(peaks.pos$chr10[wpeaks[i]], peaks.neg$chr10[wpeaks[i]])
})

#####!Differential peaks#####
cov.gfp <- coverage(resize(cstest$gfp, 200)) # Extend all reads for one lane and then report coverage
peakCutoff(cov.gfp, fdr = 0.0001) # # Estimate the peak cutoff for a specific FDR; [1] 8.846515
peaks.gfp <- slice(cov.gfp, lower = 9) # Obtain putative binding sites
peakSummary <- diffPeakSummary(peaks.gfp, peaks.ctcf) # Produce summary statistics for differentially expressed peaks

xyplot(
  asinh(sums2) ~ asinh(sums1) | seqnames,
  data = as.data.frame(peakSummary),
  panel = function(x, y, ...) {
    panel.xyplot(x, y, ...)
    panel.abline(median(y - x), 1)
  },
  type = c("p", "g"),
  alpha = 0.5,
  aspect = "iso"
)

mcols(peakSummary) <-
  within(mcols(peakSummary),
         {
           diffs <- asinh(sums2) - asinh(sums1)
           resids <- (diffs - median(diffs)) / mad(diffs)
           up <- resids > 2
           down <- resids < -2
           change <- ifelse(up, "up", ifelse(down, "down", "flat"))
         })

mcols(peakSummary) # mcols accessor (mcols stands for metadata columns) 

#####Place peaks in genomic context#####
library(TxDb.Mmusculus.UCSC.mm9.knownGene) 
gregions <- transcripts(TxDb.Mmusculus.UCSC.mm9.knownGene) 
gregions 

promoters <- flank(gregions, 1000, both = TRUE)  # Estimate the promoter for each transcript
peakSummary$inPromoter <- peakSummary %over% promoters
xtabs(~ inPromoter + change, peakSummary) # Count the peaks that fall into a promoter
#            change
# inPromoter down flat
# FALSE   16 5152
# TRUE     2  624

upstream <- flank(gregions, 20000) 
peakSummary$inUpstream <- peakSummary %over% upstream # Count the peaks that fall into the upstream
xtabs(~ inUpstream + change, peakSummary)
#            change
# inUpstream down flat
# FALSE    2 3291
# TRUE    16 2485

peakSummary$inGene <- peakSummary %over% gregions
xtabs(~ inGene + change, peakSummary)
#         change
# inGene  down flat
# FALSE    1 2716
# TRUE    17 3060

sumtab <- # Put all count results together
  as.data.frame(rbind(total = xtabs(~ change, peakSummary),
                      promoter = xtabs(~ change, 
                                       subset(peakSummary, inPromoter)),
                      upstream = xtabs(~ change, 
                                       subset(peakSummary, inUpstream)),
                      gene = xtabs(~ change, subset(peakSummary, inGene))))
sumtab
#            down flat
# total      18 5776
# promoter    2  624
# upstream   16 2485
# gene       17 3060

#####Visualize peaks in genomic context#####
library(rtracklayer) 
session <- browserSession() 
genome(session) <- "mm9" 
session$gfpCov <- cov.gfp
session$gfpPeaks <- peaks.gfp 
session$ctcfCov <- cov.ctcf
session$ctcfPeaks <- peaks.ctcf 

# View the tallest peak on chr10 in the CTCF data
peak.ord <- order(unlist(peak.depths), decreasing=TRUE) 
peak.sort <- as(peaks.ctcf, "GRanges")[peak.ord] 
view <- browserView(session, peak.sort[1], full = c("gfpCov", "ctcfCov")) 

# Display a view for the top 5 tallest peaks
views <- browserView(session, head(peak.sort, 5), full = c("gfpCov", "ctcfCov")) 
