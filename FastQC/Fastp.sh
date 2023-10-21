# Fastp is an 'all-in-one' preprocessing tool designed for use with fastq files that seamlessly incorporates various quality profiling aspects for both pre- and post-filtering data. These include quality curves, base composition, KMER analysis, Q20/Q30 evaluation, GC ratio assessment, duplication detection, adapter content, and more.

cd report
mkdir fastp_results
trim_galore --nextera -o trimmed_results ../data/file_1.fastq ../data/file_2.fastq

cd report
fastp -i ../data/file_1.fastq -I ../data/file_2.fastq \
      -o fastp_results/file_1.fastp.fastq -O fastp_results/file_1.fastp.fastq \
      --length_required 20 --average_qual 20 --detect_adapter_for_pe --correction \
      -h fastp_results/file_1.html -j fastp_results/file_2.json
