
# Build the index 
curl ftp://ftp.ensemblgenomes.org/pub/plants/release-28/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.28.cdna.all.fa.gz -o athal.fa.gz
# Build the transcriptome index

# Obtain the sequencing data
mkdir raw_data
cd raw_data
for i in `seq 1 20`; 
do 
  mkdir ABC0290${i}; 
  cd ABC0290${i}; 
  wget html_path/ABC0290${i}/ABC0290${i}_1.fastq.gz; 
  wget html_path/ABC0290${i}/ABC0290${i}_2.fastq.gz; 
  cd ..; 
done
cd .. 

# Quantify the reads
for file_path in raw_data/ABC0290{1..20}; # This loop iterates over a range of sample numbers from 1 to 20
do
    samp=`basename ${file_path}` # Extracts the sample name from the file path
    echo "Processing sample ${samp}" # Prints a message indicating which sample is being processed

    # Runs the Salmon tool to quantify gene expression for the current sample
    salmon quant -i transcriptome_index -l A \
         -1 ${fn}/${samp}_1.fastq.gz \
         -2 ${fn}/${samp}_2.fastq.gz \
         -p 8 --validateMappings -o quants/${samp}_quant
done
