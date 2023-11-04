#!/bin/bash

#SBATCH --partition=panda 
#SBATCH --nodes=5
#SBATCH --ntasks=1
#SBATCH --job-name=cellranger_count
#SBATCH --output=cellranger_count-%j.out
#SBATCH --error=cellranger_count-%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=48:00:00

# Load software modules
module load cellranger/7.0.0

# Define the paths to your input data and reference
fastqs_path=../data/pbmc_1k_v3_fastqs
reference_path=../data/refdata-gex-GRCh38-2020-A
id_name=pbmc_1k_run_count
sample_name=pbmc_1k_v3

# Run Cell Ranger count
cellranger count --id=${id_name} \
                 --transcriptome=${reference_path} \
                 --fastqs=${fastqs_path} \
                 --sample=${sample_name} \
