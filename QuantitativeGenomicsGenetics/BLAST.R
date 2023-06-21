#################################################################
##          Basic Local Alignment Search Tool (BLAST)          ##
#################################################################
##------Attach libraries------
library(rBLAST)
library(Biostrings)

##------Downlaod the 16S Microbial database from NCBI------
download.file(url = "https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz",
              destfile = "data/16S_ribosomal_RNA.tar.gz", mode = "wb") # wb: binary files
untar(tarfile = "data/16S_ribosomal_RNA.tar.gz", exdir = "data/16S_rRNA_DB") # The directory to extract files to (the equivalent of tar -C). It will be created if necessary.

##------Load sequences------
seq <- readRNAStringSet(filepath = system.file("examples/RNA_example.fasta", package = "rBLAST"))
seq

##------Load a BLAST database------
# Set the environmental variable so that R can find the blast+ executable
new.path <- paste(Sys.getenv("PATH"), "/usr/local/ncbi/blast/bin", sep = .Platform$path.sep)
Sys.setenv(PATH = new.path)
blast.db <- blast(db = "data/16S_rRNA_DB/16S_ribosomal_RNA", type = "blastn")
blast.db

##------Align a sequence using BLAST with a 99% percent identity or higher------
predict(object = blast.db, newdata = seq[1, ],
        BLAST_args = "-perc_identity 99")
