# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt

# Perform quality control analysis on particular or all FASTQ files stored in compressed (gzipped) format. 
cd scRNA/data
# View the first 6 entries
zcat HumanBrain_S1_L001_R1_001.fastq.gz | head -6 
# View the number of lines in the FASTQ file
zcat HumanBrain_S1_L001_R1_001.fastq.gz | wc -l

zcat HumanBrain_S1_L001_R1_001.fastq.gz | fastqc stdin --outdir=../report/

zcat *fastq.gz | fastqc stdin --outdir=../report/
