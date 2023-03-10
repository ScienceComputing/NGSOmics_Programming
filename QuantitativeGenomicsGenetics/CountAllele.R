# Load the data
geno_import <- read.csv("./genotype_data.csv", 
                        header = TRUE, 
                        stringsAsFactors = F,
                        row.names = 1)

for(i in 1:ncol(geno_import)) {
  geno_import[, i] <- as.character(geno_import[, i])
  geno_import[geno_import[, i] == "TRUE", i] = "T" 
}
head(geno_import[, 253:254])


# Grab the genotypes at our first two nucleotide positions
geno_in <- geno_import[, 1:2]
head(geno_in)


# The frequency of each allele at the first two nucleotide positions
geno_count <- table(c(geno_in[, 1], geno_in[, 2]))
geno_count


# Find the least frequent allele (minor allele)
minor_allele <- names(geno_count[geno_count == min(geno_count)])
minor_allele


# Count the nucleotide positions showing the minor allele across samples = minor allele frequency
xa_result <- (geno_in[,1] == minor_allele) + (geno_in[,2] == minor_allele)
xa_result


# Construct a function that calculates the frequency of minor alleles across samples at the two nucleotide positions
xa_converter <- function(geno_in){ 
  # Estimate the frequency of each allele of a genotype spanning the two nucleotide positions
  geno_count <- table(c(geno_in[, 1],geno_in[, 2]))
  
  # Find the name of the least frequent allele (minor allele)
  minor_allele <- names(geno_count[geno_count == min(geno_count)])[1]
  
  # Estimate the frequency of the minor allele of a genotype spanning the two nucleotide positions per sample
  xa_result <- (geno_in[, 1] == minor_allele) + (geno_in[, 2] == minor_allele) # Return an integer vector taking the value 0, 1, or 2 of the sample size length 
  return(xa_result)
}
