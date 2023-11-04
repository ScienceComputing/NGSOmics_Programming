#!/bin/bash

#SBATCH --partition=panda 
#SBATCH --nodes=5
#SBATCH --ntasks=1
#SBATCH --job-name=scRNAseq_analysis
#SBATCH --output=scRNA-%j.out
#SBATCH --error=scRNA-%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

# Load software modules
module load R/4.0.2

# Run the R script
Rscript scRNAseq_analysis.R
