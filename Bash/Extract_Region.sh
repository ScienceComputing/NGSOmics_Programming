#!/bin/bash

# Reference: 
# https://www.biostars.org/p/46327/
# https://genome.ucsc.edu/cgi-bin/hgGeneGraph?gene=BRCA1&1=OK&supportLevel=text&geneCount=25&1=OK&geneCount=25
# This script helps efficiently extract reads that align to a specific genomic region of interest (e.g., a gene or a regulatory element) from numerous alignment files in BAM format

N=6 # Adjust concurrency based on system resources

for i in aligned_reads*.bam; do
  samtools view "$i" chr17:41196311-41277500 -b > output_"$i" &&  
  (( ++count % N == 0)) && wait
done

