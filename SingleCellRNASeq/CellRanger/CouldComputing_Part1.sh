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
