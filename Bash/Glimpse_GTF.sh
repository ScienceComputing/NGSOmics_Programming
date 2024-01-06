# https://hgdownload.soe.ucsc.edu/downloads.html#human
grep ENST00000473358 hg38.knownGene.gtf | head -n 6
grep ENST00000473358 hg38.knownGene.gtf | wc -l # Count the number of transcript isoforms associated with the gene ENST00000473358
