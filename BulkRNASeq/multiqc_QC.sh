target_dir=/path/to/your_folder
multiqc -o $target_dir $target_dir

# In HPC, multiqc_AL.sh
#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --output=multiqc.out
#SBATCH --partition=hpc_1,hpc_2,hpc_3

export PATH="/scratch/lab/analysis_resources/tools:$PATH"
multiqc -o $1 $1
