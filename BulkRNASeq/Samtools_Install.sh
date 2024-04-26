cd /scratch/softwares
wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2
# https://github.com/samtools/samtools/releases/tag/1.19.2

# Extract the files
tar -xvjf samtools-1.19.2.tar.bz2 

# Compile samtools
cd samtools-1.19.2
./configure
make

# Install Samtools in a specific directory
make prefix=/scratch/softwares/samtools install

vim ~/.bashrc
SAMTOOLS_PATH=/scratch/softwares/samtools/bin/
export PATH=$SAMTOOLS_PATH:$PATH

exit
ssh -l your_username login03-hpc.XXX.edu

# Verify installation 
which samtools
samtools —version

