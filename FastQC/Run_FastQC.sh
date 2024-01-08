# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt

# Locate to the desired directory
cd scRNA/data

# View the first 6 entries
zcat HumanBrain_S1_L001_R1_001.fastq.gz | head -6 

# View the number of lines in the FASTQ file
zcat HumanBrain_S1_L001_R1_001.fastq.gz | wc -l

# Before trimming: perform quality control analysis on particular or all FASTQ files stored in compressed (gzipped) format
zcat HumanBrain_S1_L001_R1_001.fastq.gz | fastqc stdin --outdir=../report/

zcat *fastq.gz | fastqc stdin --outdir=../report/

#!/bin/bash
cd ~/single_cell_data/raw_fastq
for filename in *.fastq.gz; do
    echo "Processing $filename"
    zcat $filename | fastqc stdin --outdir=./report/
done

# After trimming: do the FastQC again
cat *trimmed.fastq | fastqc stdin --outdir=../report/

#!/bin/bash
cd ~/single_cell_data/trim_fastq
for filename in *trimmed.fastq; do
    echo "Processing trimmed $filename"
    cat $filename | fastqc stdin --outdir=./report/
done
