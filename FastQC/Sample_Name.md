## illumina Sequencing Data File Name
General expression: SampleName_SampleNumber_LaneNumber_ReadNumer_001.fastq.gz
- SampleName_S1_L001_R1_001.fastq.gz
- SampleName_S1_L001_R2_001.fastq.gz
- SampleName_S1_L002_R1_001.fastq.gz
- SampleName_S1_L002_R2_001.fastq.gz

```
sed -n '3p' sample_list.txt # Print the third sample name; -n: suppresse line numbering
```
