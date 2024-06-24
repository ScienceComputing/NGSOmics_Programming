# Reference: 
# https://github.com/ScienceComputing/NGSOmics_Programming/tree/main/SingleCellRNASeq/kb-python
# https://phoenixnap.com/kb/wget-command-with-examples
# http://aria2.github.io/manual/en/html/aria2c.html#synopsis
# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6108

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
# Single cell suspensions were loaded onto 10X Genomics Single Cell 3' Chips along with the reverse transcription master mix as per the manufacturer's protocol for the Chromium Single Cell 3' v2 Library (10X Genomics; PN-120233), to generate single cell gel beads in emulsion.
kb ref -d human -i human_index.idx -g human_t2g.txt --verbose --kallisto /your_path/to/kallisto/build/src/kallisto
kb count -i human_index.idx -g human_t2g.txt -x 10xv2 -o sample1 --h5ad -t 8 --filter bustools \
iPSC_RGCscRNAseq_Sample1_L005_R1.fastq.gz iPSC_RGCscRNAseq_Sample1_L005_R2.fastq.gz
