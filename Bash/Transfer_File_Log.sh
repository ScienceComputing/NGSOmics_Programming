#!/bin/bash

# Reference: https://linux.die.net/man/1/rsync

echo "Please enter your target foldername:"
read -r foldername

cd ~

if [ -d "$foldername" ]; then
    echo "$foldername exists."
else
    echo "$foldername does not exist. Creating it now."
    mkdir "$foldername"
fi

log_file="rsync.newlog.txt"

# Transfer files with logging and remove source files
rsync --remove-source-files --log-file="$foldername"/"$log_file" -avH --stats \
    *genome_bam.bam \
    *bc_matrix_h5.h5 \
    *summary.html \
    "$foldername"

echo "Transfer data to $foldername"
cd ..
