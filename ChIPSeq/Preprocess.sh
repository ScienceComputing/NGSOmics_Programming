# Further reading:
# https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#how-is-bowtie-2-different-from-bowtie-1
# http://www.htslib.org/doc/samtools-view.html

# Organize the directories
mkdir -p raw_fastq qc_reports trimmed_fastq aligned_bam sorted_bam peak_calls peak_annotations tag_directories profiles motifs

# QC
fastqc raw_fastq/*.fastq.gz -o qc_reports/

# Trim the sequences
cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGA -o trimmed_fastq/sample_trim.fastq.gz raw_fastq/sample.fastq.gz

# Align and convert SAM to BAM with Samtools
bowtie reference/GRCm38/index --phred64-quals --threads 4 -l 15 -q -S - | samtools view -bSF 4 -@ 4 > aligned_bam/sample.bam

# Sort and index BAM files
samtools sort aligned_bam/sample.bam -o sorted_bam/sample_sorted.bam
samtools -@ 4 index sorted_bam/sample_sorted.bam

# Peak calling with MACS2
macs2 callpeak -t sorted_bam/sample_sorted.bam \ 
  -c sorted_bam/input_sorted.bam \ 
  -f BAM -g mm -n sample \
  --outdir peak_calls

# Annotate the peak
annotatePeaks.pl peak_calls/sample_peaks.narrowPeak mm10 > peak_annotations/sample_annot.txt

# Make tag directories for Homer (after adding "chr" prefix to chromosome names)
samtools view -h aligned_bam/sample.bam | \
sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2'|awk -F ' ' '$7=($7=="=" || $7=="*"?$7:sprintf("chr%s",$7))' |tr " " "\t" | \
samtools reheader - aligned_bam/sample.bam > aligned_bam/sample_chr.bam
makeTagDirectory tag_directories/sample_td aligned_bam/sample_chr.bam

# Generate profiles around TSS and gene bodies
makeMetaGeneProfile.pl rna mm10 -d tag_directories/sample_td/ > profiles/sample_MetaGene.txt
annotatePeaks.pl tss mm10 -hist 50 -ghist -d tag_directories/sample_td/ > profiles/sample_TSS.txt

# Motif calling with Homer
findMotifsGenome.pl peak_calls/sample_peaks.narrowPeak mm10 motifs/sample -preparsedDir motifs/preparsed -size given -mask
