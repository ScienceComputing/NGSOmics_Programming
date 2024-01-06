#!/bin/bash

# Reference: http://www.htslib.org/doc/samtools-view.html
# -S Ignored for compatibility with previous samtools versions. 
# Previously this option was required if input was in SAM format, but now the correct format is automatically detected by examining the first few characters of input.

sam_files_count=$(find . -maxdepth 1 -name "*.sam" | wc -l)
echo "A total ${sam_files_count} SAM files are detected in the current directory."

for sam_file in *.sam; do
    # Check if the file is a regular file
    if [[ -f "$sam_file" ]]; then
        echo -n "Converting ${sam_file}... "
        bam_file="${sam_file%.sam}.bam" # Remove .sam in the value contained in sam_file
        # Convert SAM to BAM, sort and index the BAM file
        if samtools view -bS "$sam_file" | samtools sort -o "$bam_file" - && samtools index "$bam_file"; then
            echo "Done!"
        else
            echo "Error occurred while converting ${sam_file}"
        fi
    else
        echo "No sam file is processed."
    fi
done
