# https://hgdownload.soe.ucsc.edu/downloads.html#human
grep ENST00000473358 hg38.knownGene.gtf | head -n 6
grep ENST00000473358 hg38.knownGene.gtf | wc -l # Count the number of transcript isoforms associated with the gene ENST00000473358

# View the 1st and 9th columns of a GTF file
cut -f 1,9 hg38.knownGene.gtf | head -n 6
cut -f 1,9 hg38.knownGene.gtf | less

# View the contents in the 9th column seperated by ';'
cut -f 9 hg38.knownGene.gtf | cut -f 1-4 -d ';' | head -n 3
