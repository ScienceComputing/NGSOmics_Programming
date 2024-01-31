#! /bin/bash -l
#SBATCH --partition=panda
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --job-name=download_fastq
#SBATCH --time=24:00:00 # HH/MM/SS
#SBATCH --mem=32G # units: K,M,G,T
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_account@med.cornell.edu

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRRXXX/XXX/SRRXXXXXXXX/SRRXXXXXXXX.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRRXXX/XXX/SRRXXXXXXXX/SRRXXXXXXXX.fastq.gz
