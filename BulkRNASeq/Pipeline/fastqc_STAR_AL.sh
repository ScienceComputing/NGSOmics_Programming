#! /usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --job-name=demo_fastqc_STAR
#SBATCH --partition=hpc_1,hpc_2,hpc_3

ANALYSIS_DIR=$1

cd ${ANALYSIS_DIR}

mkdir -p QC BAM COUNT 

FILE_NAME=${ANALYSIS_DIR}/samples_list.txt
INPUT_LISTS=${FILE_NAME}

readarray -t samples < ${INPUT_LISTS} 
CURRENT_INPUT=${samples[$((SLURM_ARRAY_TASK_ID-1))]}

echo processing : ${CURRENT_INPUT}

if [ -s BAM/${CURRENT_INPUT}_starAligned.sortedByCoord.out.bam ]; then
       echo "The file BAM/${CURRENT_INPUT}_starAligned.sortedByCoord.out.bam exists with > 0 size."
else

    export PATH=/scratch/analysis_resources/tools/FASTQC/0.11.7/binary/FastQC:$PATH

    fastqc FASTQ/${CURRENT_INPUT}* -t 12 -o demo_QC 

    IDX="/scratch/analysis_resources/tools/STAR/STAR-2.7.10a/genome"

    export PATH=/scratch/analysis_resources/tools/STAR/STAR-2.7.10a/bin/Linux_x86_64_static:$PATH

    STAR    --runThreadN 12 \
                --genomeDir ${IDX}\
                --readFilesIn FASTQ/${CURRENT_INPUT}_R1_001.fastq.gz FASTQ/${CURRENT_INPUT}_R2_001.fastq.gz \
                --readFilesCommand zcat \
                --outFileNamePrefix BAM/${CURRENT_INPUT}_star \
                --outSAMtype BAM SortedByCoordinate \
                --outReadsUnmapped Fastx \
                --limitBAMsortRAM 30000000000 \
                --twopassMode Basic \
                --outSAMunmapped Within \
                --outFilterMultimapNmax 10 \
                --genomeSAindexNbases 10
fi
