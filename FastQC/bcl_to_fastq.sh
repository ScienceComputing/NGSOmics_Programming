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
# y = read, i = index, n = ignore; the default values are assumed from the RunInfo.xml file 
# Usage for all cycles: --use-bases-mask=Y*,I*,Y*
# -R ${FLOWCELL_DIR}: specify the directory that contains a flow cell's Data folder
# --output-dir=${OUTPUT_DIR}: specify the directory that we want to output the converted FASTQs to
# --interop-dir=${INTEROP_DIR}: specify the directory where InterOp files are located. InterOp files contain various run metrics.
# --sample-sheet=${SAMPLE_SHEET_PATH}: specify the path to the sample sheet CSV we create
