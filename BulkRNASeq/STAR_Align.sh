# https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
# http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
# The workflow of STAR: 1) build a reference genome index; 2) map reads to the reference genome
# STAR’s default parameters are optimized for mammalian genomes (e.g., human genome)

wget https://github.com/alexdobin/STAR/archive/2.7.11a.tar.gz
tar -xzf 2.7.11a.tar.gz
cd STAR-2.7.11a
alias STAR="/Users/your_name/STAR-2.7.11a/bin/MacOSX_x86_64/STAR"
# echo "export PATH=/Users/your_name/STAR-2.7.11a/bin/MacOSX_x86_64:$PATH" > .bashrc
# source .bashrc
which STAR

# or STAR=/Users/your_name/STAR-2.7.11a/bin/MacOSX_x86_64/STAR

# Step 1: build a reference genome index
mkdir -p /Users/your_name/bulk_RNAseq/data/hg38_chr1_index

STAR --runThreadN 20 \
--runMode genomeGenerate \ 
--genomeDir /Users/your_name/bulk_RNAseq/data/hg38_chr1_index \
--genomeFastaFiles /Users/your_name/bulk_RNAseq/data/reference/homo_sapiens.GRCh38.dna.chromosome.1.fa \
--sjdbGTFfile /Users/your_name/bulk_RNAseq/data/reference/homo_sapiens.GRCh38.92.gtf \
--sjdbOverhang 99

# Alternatively
# ${STAR} --runThreadN 20 \
# --runMode genomeGenerate \ 
# --genomeDir /Users/your_name/bulk_RNAseq/data/hg38_chr1_index \
# --genomeFastaFiles /Users/your_name/bulk_RNAseq/data/reference/homo_sapiens.GRCh38.dna.chromosome.1.fa \
# --sjdbGTFfile /Users/your_name/bulk_RNAseq/data/reference/homo_sapiens.GRCh38.92.gtf \
# --sjdbOverhang 99

# genomeDir: the path to store the genome indices
# genomeFastaFiles: the path to one or more reference genomes in FASTA format
# sjdbGTFfile: the path to the annotations of the reference genome in GTF format
# sjdbOverhang: the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. 
# Ideally, this length should equal to the max(ReadLength)-1, where ReadLength is the length of the reads. e.g., for Illumina 2x100b paired-end reads, the ideal value is 100-1=99. 
# Often, 100 works quite well.

# Step 2: map reads to the reference genome
STAR --genomeDir /Users/your_name/bulk_RNAseq/data/hg38_chr1_index \
--runThreadN 20 \
--runMode alignReads
--readFilesIn /Users/your_name/bulk_RNAseq/data/fastq/s1read1.fq.gz,s2read1.fq.gz s1read2.fq.gz,s2read2.fq.gz \ # !
--readFilesCommand zcat
--outFileNamePrefix /Users/your_name/bulk_RNAseq/results/STAR/Human_Infectious_Disease_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard 

# use --readFilesCommand bunzip2 -c for bzip2 compressed files



