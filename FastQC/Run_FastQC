# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt

# Perform quality control analysis on FASTQ files stored in compressed (gzipped) format. 
cd scRNA/data
zcat HumanBrain_S1_L001_R1_001.fastq.gz | fastqc stdin --outdir=../report/
zcat *fastq.gz | fastqc stdin --outdir=../report/
