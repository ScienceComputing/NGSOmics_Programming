# Reference: https://github.com/ScienceComputing/NGSOmics_Programming/tree/main/SingleCellRNASeq/kb-python
# https://phoenixnap.com/kb/wget-command-with-examples
# http://aria2.github.io/manual/en/html/aria2c.html#synopsis

vim target_fastq.txt
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample1_L005_R1.fastq.gz
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample1_L005_R2.fastq.gz
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample2_L005_R1.fastq.gz
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample2_L005_R2.fastq.gz
wget -q -i target_fastq.txt
# -i input file; -q: suppress output

# A faster downloading of large files
brew install aria2
aria2c -x 16 -s 16 https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample1_L005_R1.fastq.gz
aria2c -x 16 -s 16 https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample1_L005_R2.fastq.gz
aria2c -x 16 -s 16 https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample2_L005_R1.fastq.gz
aria2c -x 16 -s 16 https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample2_L005_R2.fastq.gz

# Quantify the read counts from paired-end FASTQ files

