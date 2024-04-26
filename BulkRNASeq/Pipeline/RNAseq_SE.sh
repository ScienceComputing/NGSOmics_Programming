#!/bin/bash
#SBATCH --job-name=RNAseq_SE
#SBATCH --output=RNAseq_SE.out
#SBATCH --partition=hpc_1,hpc_2,hpc_3

ANALYSIS_DIR=$1

cd ${ANALYSIS_DIR}/FASTQ
ls *_R1* > R1_files_list.txt                           

awk '{gsub(/_R1_001.fastq.gz$/, ""); print > "../samples_list.txt"}' R1_files_list.txt

rm R1_files_list.txt

WORK_DIR=$HOME/pipeline_AL/
LOG_DIR=${ANALYSIS_DIR}/logs_fastqc_STAR/

mkdir -p ${LOG_DIR}

cd ${WORK_DIR}

OUT_1=$(sbatch --array=1-$(wc -l < ${ANALYSIS_DIR}/samples_list.txt) -o ${LOG_DIR}/fastqc_STAR_%a.out ${WORK_DIR}/fastqc_STAR_AL.sh ${ANALYSIS_DIR})

OUT_2=$(sbatch --dependency=afterok:${OUT_1##* } ${WORK_DIR}/featureCounts_AL.sh ${ANALYSIS_DIR})

OUT_3=$(sbatch --dependency=afterok:${OUT_2##* } ${WORK_DIR}/multiqc_AL.sh ${ANALYSIS_DIR})
