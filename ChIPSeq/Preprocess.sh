# Organize the directories
mkdir -p raw_fastq qc_reports trimmed_fastq aligned_bam sorted_bam peak_calls peak_annotations tag_directories profiles motifs

# QC
fastqc raw_fastq/*.fastq.gz -o qc_reports/

# Trim the sequences
cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGA -o trimmed_fastq/sample_trim.fastq.gz raw_fastq/sample.fastq.gz
