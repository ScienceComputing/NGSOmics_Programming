#####Make custom visuals of methylation data#####
#####Set up the genomic region#####
gen <- "hg19" # Indicate which genome is being used
dmr.index <- 1 # The index of the DMR that will be plot

# Extract chromosome number and location from DMR results 
chrom <- as.character(seqnames(results.ranges[dmr.index]))
start <- as.numeric(start(results.ranges[dmr.index]))
end <- as.numeric(end(results.ranges[dmr.index]))

# Add 25% extra space to plot
minbase <- start - (0.25 * (end - start))
maxbase <- end + (0.25 * (end - start))

#####Add genomic annotations of interest - locations of CpG islands#####
island.HMM <- read.csv(paste0(data.dir,
                              "/model-based-cpg-islands-hg19-chr17.txt"),
                       sep="\t", stringsAsFactors=F, header=F)
head(island.HMM)

island.data <- GRanges(seqnames=Rle(island.HMM[,1]), 
                       ranges=IRanges(start=island.HMM[,2], end=island.HMM[,3]),
                       strand=Rle(strand(rep("*",nrow(island.HMM)))))
island.data

#####Add genomic annotations of interest - DNAseI hypersensitive sites#####
dnase <- read.csv(paste0(data.dir,"/wgEncodeRegDnaseClusteredV3chr17.bed"),
                  sep="\t",stringsAsFactors=F,header=F)
head(dnase)

dnase.data <- GRanges(seqnames=dnase[,1],
                      ranges=IRanges(start=dnase[,2], end=dnase[,3]),
                      strand=Rle(rep("*",nrow(dnase))),
                      data=dnase[,5])
dnase.data

#####Set up the ideogram, genome and RefSeq tracks#####
i.track <- IdeogramTrack(genome = gen, chromosome = chrom, name="")
# Ideogram represents the schematic display of a chromosome
g.track <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
# Genomeaxis represents a customizable genomic axis
r.track <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq", 
                     from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                     rstarts="exonStarts", rends="exonEnds", gene="name", 
                     symbol="name", transcript="name", strand="strand", 
                     fill="darkblue",stacking="squish", name="RefSeq", 
                     showId=T, geneSymbol=T)
# This returns an annotation track object 
# track: Character, the name of the track to fetch from UCSC
# genome: Character, a valid USCS genome identifier for which to fetch the data
# chromosome: Character, a valid USCS character identifier for which to fetch the data
# from, to: A range of genomic locations for which to fetch data

#####Order the methylation data by chromosome and base position#####
anndata.450k.ord <- anndata.450k.sub[order(anndata.450k.sub$chr,anndata.450k.sub$pos),]
head(anndata.450k.ord)

b.val.ord <- b.val[match(anndata.450k.ord$Name,rownames(b.val)),]
head(b.val.ord)

#####Create the data tracks#####
# Create genomic ranges object from methylation data
cpg.data <- GRanges(seqnames=Rle(anndata.450k.ord$chr),
                    ranges=IRanges(start=anndata.450k.ord$pos, end=anndata.450k.ord$pos),
                    strand=Rle(rep("*",nrow(anndata.450k.ord))),
                    betas=b.val.ord)

# Extract the data on CpGs in DMR
cpg.data <- subsetByOverlaps(cpg.data, results.ranges[dmr.index])

# Build the CpG island track
island.track <- AnnotationTrack(range=island.data,genome=gen,name="CpG Is.", 
                                chromosome=chrom,fill="darkgreen")

# Build the methylation data track
meth.track <- DataTrack(range=cpg.data, groups=sample.info$Sample_Group,genome = gen,
                        chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                        type=c("a","p"), name="DNA Meth.\n(beta value)",
                        background.panel="white", legend=T, cex.title=0.8,
                        cex.axis=0.8, cex.legend=0.8)

# Build the DNaseI hypersensitive site data track
dnase.track <- DataTrack(range=dnase.data, genome=gen, name="DNAseI", 
                         type="gradient", chromosome=chrom)

# Build the DMR position data track
dmr.track <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                             chromosome=chrom,fill="darkred")

# Combine all data tracks in one plot
track.list <- list(i.track, g.track, meth.track, dmr.track, island.track, dnase.track,
                   r.track)
sizes <- c(2,2,5,2,2,2,3) # Set up the relative sizes of the tracks
plotTracks(track.list, from=minbase, to=maxbase, showTitle=T, 
           add53=T, # 5' -> 3'
           add35=T, # 3' -> 5'
           grid=T, lty.grid=3, sizes = sizes, length(track.list))
# This plot shows one of the DMRs identified by the DMRcate analysis