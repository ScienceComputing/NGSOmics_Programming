# https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
# The workflow of STAR: 1) build a reference genome index; 2) map reads to the reference genome

# Download latest STAR source
wget https://github.com/alexdobin/STAR/archive/2.7.11a.tar.gz
tar -xzf 2.7.11a.tar.gz
cd STAR-2.7.11a

# Compile under Mac OS X
brew install gcc
cd source
make STARforMacStatic CXX=/usr/local/Cellar/gcc/8.2.0/bin/g++-8
cp STAR /usr/local/bin

# Build a reference genome index
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir chr1_hg38_index \
--genomeFastaFiles /n/groups/hbctraining/intro_rnaseq_hpc/reference_data_ensembl38/Homo_sapiens.GRCh38.dna.chromosome.1.fa \
--sjdbGTFfile /n/groups/hbctraining/intro_rnaseq_hpc/reference_data_ensembl38/Homo_sapiens.GRCh38.92.gtf \
--sjdbOverhang 99
