# Use FastQC to check read quality again after trimming; FastQC should show that our reads have nice quality on the 'Per Base Sequence Quality' and 'Adaptor Content' plots 

cd report
mkdir trimmed_results
trim_galore --nextera -o trimmed_results ../data/file_1.fastq ../data/file_2.fastq
