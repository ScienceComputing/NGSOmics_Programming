featureCounts -T 20 \
-t exon \
-g gene_id \
-a annotation.gtf \
-o human_featureCounts_output.txt \
sorted_human_alignment.bam
