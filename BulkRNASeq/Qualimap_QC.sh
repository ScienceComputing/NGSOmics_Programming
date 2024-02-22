qualimap rnaseq \
-outdir /path/bulk_RNAseq/results/qualimap/Human_Infectious_Disease \
-a proportional \
-bam /path/bulk_RNAseqresults/STAR/Human_Infectious_Disease_Aligned.sortedByCoord.out.bam \
-p strand-specific-reverse \
-gtf /path/bulk_RNAseq/reference_data/homo_sapiens.GRCh38.92.gtf \
--java-mem-size=20G

# -outdir: provide the output directory for html report
# -a: provide the option for the counting algorithm - uniquely-mapped-reads(default) or proportional (each multi-mapped read is weighted according to the number of mapped locations)
# -bam: provide the directory that stores the BAM file
# -p: provide the option for the sequencing library protocol: 1) strand-specific-forward, 2) strand-specific-reverse, 3) non-strand-specific (default)
# -gtf: provide the directory that stores the GTF used in previous alignment
# --java-mem-size=: set Java memory
