# Print the first 3 sequecing reads
gzip -cd SRR_number.fastq.gz | head -n 12

# View the whole records of N bases in one/multiple FASTQ files
grep -B 1 -A 2 N SRR_number.fastq # -B 1: 1 line before the matched line; -A 2: 2 lines after the matche line
grep -B 1 -A 2 -n N SRR_number.fastq # -n: show the line number;  `line_number :` -> lines with the match
grep -B 1 -A 2 -n N SRR_number.fastq | less #!
grep -B 1 -A 2 -n N SRR_number.fastq | head -n 12
grep -B 1 -A 2 N SRR*
grep -B 1 -A 2 -n N SRR*
grep -B 1 -A 2 -n N SRR* | less
grep -B 1 -A 2 -n N SRR* | head -n 12
grep -B 1 -A 2 -n N SRR_number.fastq > low_quality_read.txt # > overwrite
grep -B 1 -A 2 -n N SRR* >> low_quality_read.txt # >> append
cat low_quality_read.txt

# ! View the whole records of N bases in the compressed FASTQ file
gzip -cd SRR_number.fastq.gz | grep -B 1 -A 2 N
gzip -cd SRR_number.fastq.gz | grep -B 1 -A 2 N | less #!
gzip -cd *.fastq.gz | grep -B 1 -A 2 N
gzip -cd *.fastq.gz | grep -B 1 -A 2 -n N
gzip -cd *.fastq.gz | grep -B 1 -A 2 -n N | less 
gzip -cd *.fastq.gz | grep -B 1 -A 2 -n N | head -n 12
gzip -cd SRR_number.fastq.gz | grep -B 1 -A 2 -n N > low_quality_read.txt
gzip -cd *.fastq.gz | grep -B 1 -A 2 -n N >> low_quality_read.txt

# Remove the group separator (--), when viewing the whole records of N bases in the compressed FASTQ file
brew install grep
gzip -cd SRR_number.fastq.gz | ggrep -B 1 -A 2 --no-group-separator N
gzip -cd *.fastq.gz | ggrep -B 1 -A 2 --no-group-separator N
gzip -cd SRR_number.fastq.gz | ggrep -B 1 -A 2 --no-group-separator N > low_quality_read.txt
gzip -cd *.fastq.gz | ggrep -B 1 -A 2 --no-group-separator N >> low_quality_read.txt

# Count the number of reads detected with a N base in one/multiple FASTQ/compressed FASTQ files
grep N SRR_number.fastq | wc -l
grep N *.fastq | wc -l
gzip -cd SRR_number.fastq.gz | grep N | wc -l
gzip -cd *.fastq.gz | grep N | wc -l


