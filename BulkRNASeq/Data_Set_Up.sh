#!/bin/bash
#SBATCH --job-name=data_set_up
#SBATCH --output=data_set_up.out
#SBATCH --partition=hpc_1,hpc_2,hpc_3

PROJECT_DIR=$1

DATA_DIR=$2

ANALYSIS_DIR=$3

if [ -d "$DATA_DIR" ]; then # Enclose the $DATA_DIR in double-quotes to handle directories having white/black spaces in their names
  echo "The real data has been in the specified data directory; no need to create the symbolic links projecting from project directory to data directory."
  mkdir -p  $ANALYSIS_DIR/FASTQ/
  ln -s $DATA_DIR/*.fastq.gz  $ANALYSIS_DIR/FASTQ/
else
  mkdir -p $DATA_DIR
  cd $PROJECT_DIR
  cp *.fastq.gz $DATA_DIR
  mkdir $ANALYSIS_DIR
  cd $ANALYSIS_DIR
  mkdir -p FASTQ
  ln -s $DATA_DIR/*.fastq.gz FASTQ/ 
fi
