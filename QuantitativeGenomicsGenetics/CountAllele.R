# Load the data
geno_import <- read.csv("./genotype_data.csv", 
                        header = T, 
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


# Apply the function `xa_converter()` to each genotype
xa_matrix <- matrix(NA, nrow = nrow(geno_import), ncol = ncol(geno_import)/2)
  # What are the rows and columns? row: sample; column: nucleotide position
  # Why divide by 2? a genotype spanning two nucleotide positions: A1A2

for (i in 1:(ncol(geno_import)/2)){
    xa_matrix[, i] <- xa_converter(geno_import[, c(2*i - 1, 2*i)]) # Grab pairs of columns; (1,2), (3,4)...(399,400); Each column of xa_matrix is an integer vector taking the value 0, 1, or 2 of the sample size length 
}
xa_matrix <- xa_matrix - 1 # Center on 0
xa_matrix[1:10, 1:10]


# Write an expanded function that handle diverse situations in the class and dim of the imported data, the frequencies of minor and major alleles, the threshold of the minor allele frequency, and dummary variables in the additive genetic model
xa_converter <- function(geno_in, maf_limit = 0.01, return_minor_allele = F){
  # Combine a pair of two nucleotide positions to each genotype which is a vector of double letter elements: AA, AG, GG, ...
  if(class(geno_in) == "matrix" | class(geno_in) == "data.frame"){
    if(ncol(geno_in) > 1){
      geno <- paste0(geno_in[, 1], geno_in[, 2])
    } else {
      geno <- as.vector(geno_in)
    }
  } else {
    geno <- geno_in
  }
  
  return_coding <- rep(NA, length(geno))
  
  # Determine the minor and major alleles
  geno_undone <- unlist(strsplit(paste0(geno_in[, 1], geno_in[, 2]), ''))
  geno_count <- table(geno_undone)
  maf <- min(geno_count)[1]/length(geno) # maf: minor allele frequency
  if(length(geno_count) == 1 | maf < maf_limit){
   return(return_coding) 
  } else if(geno_count[1] == geno_count[2]){
    minor_allele <- sample(names(geno_count), 1)
    major_allele <- names(geno_count)[names(geno_count) != minor_allele]
  } else {
    minor_allele <- names(geno_count[geno_count == min(geno_count)])
    major_allele <- names(geno_count[geno_count == max(geno_count)])
  }

  # Categorize genotypes to the "dummy variables" centered on 0 (for additive genetic model)
  # The additive genetic model assumes that the effects of two alleles at a given locus (location on a chromosome) are additive, meaning that the more copies of the "good" allele a person has, the better their trait value will be.
  return_coding[geno == paste0(minor_allele, minor_allele)] <- 1
  return_coding[geno == paste0(minor_allele, major_allele) | geno == paste0(major_allele, minor_allele) ] <- 0
  return_coding[geno == paste0(major_allele, major_allele)] <- -1
  
  if(return_minor_allele){
    return(list(return_coding, minor_allele))
  } else {
    return(return_coding)
  }
}

xa_matrix <- matrix(NA, nrow = nrow(geno_import), ncol = ncol(geno_import)/2)

for (i in 1:(ncol(geno_import)/2)){
    xa_matrix[, i] <- xa_converter(geno_import[,c(2*i - 1, 2*i)]) # Grab pairs of columns; (1,2), (3,4)...(399,400); Each column of xa_matrix is the returned vector `return_coding`, which is an integer vector taking the value 1, 0, or -1 of the sample size length   
}

xa_matrix[1:10, 1:10]

