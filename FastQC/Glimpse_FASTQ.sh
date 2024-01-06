# Print the first 3 sequecing reads
gzip -cd SRR_number.fastq.gz | head -n 12

# View N bases
grep -B 1 -A 2 N SRR_number.fq



