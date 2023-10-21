# PAGE 14: https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf
# https://www.youtube.com/embed/gz3H_pcLe7c?autoplay=1&rel=0

bcl2fastq --use-bases-mask=Y26,I8,Y98 \
            --create-fastq-for-index-reads \
            --minimum-trimmed-read-length=8 \
            --mask-short-adapter-reads=8 \
            --ignore-missing-positions \
            --ignore-missing-controls \
            --ignore-missing-filter \
            --ignore-missing-bcls \
            -r 6 -w 6 \
            -R ../flow_cell/data \
            --output-dir=../fastq \
            --interop-dir=../op \
            --sample-sheet=../sample_sheet.csv

# --use-bases-mask: specify the number of cycles to use in each read and how to use them
# y = read, i = index, n = ignore
# If --use-bases-mask is not specified, the default values will be taken from the RunInfo.xml file 
# Usage for all cycles: --use-bases-mask=Y*,I*,Y*
# Y26,I8,Y98: consider 26 bases for the first read, 8 bases for the index read, and 98 bases for the second read

# -r 6 -w 6: both reading BCL files and writing FASTQ files will use 6 threads
# -R ${FLOWCELL_DIR}: specify the directory that contains a flow cell's data folder
# --output-dir=${OUTPUT_DIR}: specify the directory that we want to output the converted FASTQs to
# --interop-dir=${INTEROP_DIR}: specify the directory where InterOp files are located. InterOp files contain various run metrics.
# --sample-sheet=${SAMPLE_SHEET_PATH}: specify the path to the sample sheet CSV we create
