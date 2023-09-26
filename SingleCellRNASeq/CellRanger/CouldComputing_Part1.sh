# Download the 10x Genomics Cloud command line interface and unpack it
curl -f -o txg-macos-v1.3.1.zip "https://cf.10xgenomics.com/cloud-cli/v1.3.1/txg-macos-v1.3.1.zip"
tar -zxvf txg-macos-v1.3.1.tar.gz
cd /Users/your_name/txg-macos-v1.3.1

# Check if we can run the executable txg
./txg --version
# txg version v1.3.1

# Add the executable txg into our $PATH
touch ./bashrc
echo "export PATH=$PATH:/Users/your_name/txg-macos-v1.3.1" > .bashrc
cat .bashrc
source .bashrc
where txg

# Upload all the FASTQ files in one folder
# Notice that all FASTQ files must follow the Illumina naming convention (e.g., samplename_S1_L001_R1_001.fastq.gz)
# More details on the naming convention can be viewed here: https://help.basespace.illumina.com/files-used-by-basespace/fastq-files
txg fastqs upload --project-id hHCuUsZhQJeSSSD90699s6gg scRNA/data/ 

# Upload one FASTQ file
txg fastqs upload --project-id hHCuUsZhQJeSSSD90699s6gg scRNA/data/Human_S1_L001_R1_001.fastq.gz

# Upload all the FASTQ files of which the names start with Neuron
txg fastqs upload --project-id hHCuUsZhQJeSSSD90699s6gg scRNA/data/Neuron*

# Now, we'are going to fetch the information and upload the snRNA-seq data (SRR17380406) from one human brain sample
# This data comes from the recent Nature paper: molecular features driving cellular complexity of human brain evolution
ffq SRR17380406
ffq --ftp SRR17380406 | jq -r '.[] | .url' | xargs curl -O
# Or we can use the NCBI SRA Toolkit to download the target FASTQ.gz data
prefetch SRR17380406 --max-size u
fasterq-dump SRR17380406 -e 20 # Here we use 20 threads to speed up the downloading of this large-size file

# We put this fastq.gz file onto the 10x Genomics clouding analytics platform
txg fastqs upload --project-id hHCuUsZhQJeSSSD90699s6gg scRNA/data/HumanBrain_S1_L001_R1_001.fastq.gz
