#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --output=multiqc.out
#SBATCH --partition=hpc_1,hpc_2,hpc_3

export PATH="/scratch/lab/analysis_resources/tools:$PATH"
multiqc -o $1 $1
