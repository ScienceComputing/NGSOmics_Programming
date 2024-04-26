#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --output=featureCounts.out
#SBATCH --partition=hpc_1,hpc_2,hpc_3


ANALYSIS_DIR=$1

export PATH=/scratch/analysis_resources/tools/Subread/subread-2.0.3-Linux-x86_64/bin:$PATH

featureCounts -T 20 -p -g gene_name \ 
  -a /scratch/reference_genome/Homo_sapiens.GRCh38.103.gtf \
  -o ${ANALYSIS_DIR}/COUNT/featurecounts.txt \
  ${ANALYSIS_DIR}/BAM/*.bam
