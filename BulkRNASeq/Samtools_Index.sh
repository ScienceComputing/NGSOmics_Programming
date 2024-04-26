ANALYSIS_DIR=/RNAseq/analysis/Project_IRF1
ls $ANALYSIS_DIR/BAM
cd $ANALYSIS_DIR/BAM

# Create a folder containing both BAM and BAI files
mkdir -p $ANALYSIS_DIR/BAMBAI
mkdir -p $ANALYSIS_DIR/BAMBAI_2

ls | grep '^XYZI1.*\.bam$' | xargs -I {} cp {} $ANALYSIS_DIR/demo_BAMBAI
ls | grep '^CTL1.*\.bam$' | xargs -I {} cp {} $ANALYSIS_DIR/demo_BAMBAI_2

# cd $ANALYSIS_DIR/BAMBAI
# ls | grep '^CTL1.*\.bam$' | xargs rm 
# rm CTL*

ls $ANALYSIS_DIR/BAMBAI

ls $ANALYSIS_DIR/BAMBAI_2

# Run samtools index in the HPC environment
# Create the script
vim $HOME/pipeline_AL/SamIndex_AL.sh

#!/bin/bash
#SBATCH --job-name=demo_SamIndex
#SBATCH --output=demo_SamIndex.out
#SBATCH --partition=hpc,casa,bigmem

samtools index -M $1/*.bam

cd $ANALYSIS_DIR/BAMBAI
sbatch $HOME/pipeline_AL/SamIndex_AL.sh $ANALYSIS_DIR/BAMBAI

cd $ANALYSIS_DIR/BAMBAI_2
sbatch $HOME/pipeline_AL/SamIndex_AL.sh $ANALYSIS_DIR/BAMBAI_2

squeue -u your_username -l # Show the info of the current job
sacct | tail -n 3 # Show the info of the past job

# Transfer the BAM and BAI files from HPC to local machines
# Open the local terminal

# Transfer the whole directory
rsync -avP your_username@login03-hpc.XXX.edu://$ANALYSIS_DIR/BAMBAI/'*' ./Documents/Project_IRF1/data/raw/SamIndex/ 
rsync -avP your_username@login03-hpc.XXX.edu://$ANALYSIS_DIR/BAMBAI_2/'*' ./Documents/Project_IRF1/data/raw/SamIndex/ 
