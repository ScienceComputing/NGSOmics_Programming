# Set up the environment
work_dir = "/Users/your_name/SingleCell"
mkdir data result
cd $work_dir/result/

# Download the FASTQ data and prebuilt human reference transcriptome
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar -xvf pbmc_1k_v3_fastqs.tar
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -zxvf refdata-gex-GRCh38-2020-A.tar.gz

# Align reads and generate Feature Barcode matrices
cellranger count --id=pbmc_1k_count \
   --fastqs=../data/pbmc_1k_v3_fastqs \
   --sample=pbmc_1k_v3 \
   --transcriptome=../data/refdata-gex-GRCh38-2020-A

# List the output
ls -1 pbmc_1k_count/outs # Display one entry per line
ls -Flha pbmc_1k_count/outs # Thoroughly examine the output
du -sh
