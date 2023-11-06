#!/bin/bash
trimmed_fastq_dir="trimmed_fastq"
bowtie_ref_dir="reference/GRCm38"
alignment_dir="aligned_bam"

mkdir -p "$alignment_dir"

for fq_gz in "$trimmed_fastq_dir"/*.fastq.gz; do
  filehead=$(basename "$fq_gz" .fastq.gz)
  echo "Processing $filehead"
  zcat "$fq_gz" |
  bowtie "$bowtie_ref_dir" --phred33 --threads 4 -l 15 -q -S - |
  samtools view -bSF 4 -@ 4 > "$alignment_dir/${filehead}_aligned.bam"
done
