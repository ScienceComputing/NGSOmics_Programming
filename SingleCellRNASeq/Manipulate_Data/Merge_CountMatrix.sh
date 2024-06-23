vim target_fastq.txt
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample1_L005_R1.fastq.gz
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample1_L005_R2.fastq.gz
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample2_L005_R1.fastq.gz
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6108/iPSC_RGCscRNAseq_Sample2_L005_R2.fastq.gz
wget -q -i target_fastq.txt
# -i input file; -q: suppress output


