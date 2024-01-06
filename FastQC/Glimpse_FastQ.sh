# Print the first 3 sequecing reads
gzip -cd SRR_number.fastq.gz | head -n 12

# View the whole records of N bases in one/multiple FASTQ files
grep -B 1 -A 2 N SRR_number.fq
grep -B 1 -A 2 N SRR*

# View the whole records of N bases in the compressed FASTQ file
gzip -cd SRR_number.fastq.gz | grep -B 1 -A 2 N



