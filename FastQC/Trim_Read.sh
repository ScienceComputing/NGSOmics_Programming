# Use FastQC to check read quality again after trimming; FastQC should show that our reads have nice quality on the 'Per Base Sequence Quality' and 'Adaptor Content' plots 
# https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md

cd report
mkdir trimmed_results
trim_galore --illumina --paired -o trimmed_results ../data/file_1_trimmed.fastq ../data/file_2_trimmed.fastq
