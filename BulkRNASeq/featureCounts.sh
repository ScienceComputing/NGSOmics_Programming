# Summarize paired-end reads and count fragments 
featureCounts -T 20 \
-p --countReadPairs
-t exon \
-g gene_id \
-a annotation.gtf \
-o human_featureCounts_output.txt \
Human_Infectious_Disease_Aligned.sortedByCoord.out.bam

# T: use 20 threads
