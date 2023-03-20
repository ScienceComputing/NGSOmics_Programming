# Goal: split a pair of two consecutive alleles stored in a column to two separate columns; repeat this splitting for each original column 
# so that the column number of the output genotype dataset should be double the original genotype dataset

# Load genotype data into R
geno_import <- read.delim(file = "QG23 - hw4_genotypes.txt",
                          header = F,
                          sep = " ")

# Split genotype table to get even numbered and odd numbered rows
evenNumbRowsGenotypes <- geno_import[seq(2, nrow(geno_import), 2), ]
oddNumbRowsGenotypes <- geno_import[seq(1, nrow(geno_import), 2), ]
# oddNumbRowsGenotypes[1:10, 1:6] 

# Create a function to create alleles (2 bases) from genotypes
transformGeno <- function(x, y) {
  genoOut <- paste0(x, y)
}

# Make alleles per loci for each individual
mGenotype <- mapply(transformGeno, oddNumbRowsGenotypes, evenNumbRowsGenotypes)

# Rename columns as "genotypes"
colnames(mGenotype) <- paste0("genotype", seq(ncol(mGenotype))) 
# mGenotype[1:10,1:6] 

# Separate two alleles in a cell into two cells (each of two cells carries one allele) per column
analysisReady_genotype_df <- mGenotype %>% 
  as.data.frame() %>% 
  mutate(across(everything(), ~ str_split(.x, "", simplify = T))) %>% 
  as.matrix() %>% 
  as.data.frame()
