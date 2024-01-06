wget -c \
ftp://ftp.sra.ebi.ac.uk/XXXX.fastq.gz \
ftp://ftp.sra.ebi.ac.uk/YYYY.fastq.gz \ 
ftp://ftp.sra.ebi.ac.uk/ZZZZ.fastq.gz \
ftp://ftp.sra.ebi.ac.uk/AAAA.fastq.gz


# Inside download_pair_fastq.sh
#!/bin/bash

read -p "Enter SRR IDs separated by spaces: " srr_ids
IFS=' ' read -r -a srr_array <<< "$srr_ids"

for srr_id in "${srr_array[@]}"; do
  echo "Downloading paired-end FASTQ files for SRR ID $srr_id..."

  wget -c \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRRXXX/${srr_id}_1.fastq.gz \ # Customize
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRRXXX/${srr_id}_2.fastq.gz

  echo "Complete downloading for SRR ID $srr_id!"
done
