# Goal: covert the ATCGs into dummy numbers under either additive or dominant genetic model

# Input genotype file: each column represents a specific SNP (column 1 = column 1 = genotype 1, column 2 = genotype 2)
# Each consecutive pair of rows represent all of the genotype states for an individual for the entire set of SNPs
# For example, rows 1 and 2 = all of individual 1’s genotypes, rows 3 and 4 = all individual 2’s genotypes

geno.recode <- function(geno.import, maf.lim, error.val = 3) {
  
  # Inputs: geno.col is a vector of genotype where each consecutive pair of two elements in this vector represents one particular genotype for an individual), num.sample is the number of individuals, maf.lim is the threshold of minor allele frequency (MAF)
  xa.convert <- function(geno.col, num.sample, maf.lim) {
    
    # Frequencies of alleles at this locus
    geno.count <- table(geno.col)
    
    # If the MAF is less than the pre-set threshold OR this locus has less than 2 alleles
    if (min(geno.count) / length(geno.col) <= maf.lim | length(geno.count) < 2) {
      
      # Return a vector of 3's, to be filtered out later
      return(rep(error.val, num.sample))
      
    }
    
    # Otherwise, figure out which is the minor allele, and grab the first element in case we end up with a tie
    minor.allele <- names(geno.count[geno.count == min(geno.count)])[1]
    
    # geno.col[seq(1, num.sample*2 - 1, 2)] refers to the first allele in a genotype across individuals
    # geno.col[seq(2, num.sample*2, 2)] refers to the second allele in a genotype across individuals
    # The conditional gives a Boolean (T/F) but the vector operation (+) converts it to (1/0)
    # For example, an individual with TT where T is the minor allele would have TRUE + TRUE = 2; an individual with AT where T is the minor allele would have FALSE + TRUE = 1
    xa <- (geno.col[seq(1, num.sample*2 - 1, 2)] == minor_allele) + (geno.col[seq(2, num.sample*2, 2)] == minor.allele)
    # xa originally takes the value of 0, 1, 2; now we rescale it so that it takes the value of -1, 0, or 1
    xa <- xa - 1
    
    return(xa)
  }
  
  # Apply the xa.convert function to each column of geno.import, with num.sample = nrow(geno.import)/2, and maf.limit = 0.05; this returns a matrix of xa dummy variable coding for all genotypes and the row number of this matrix is the number of individuals
  xa.mat  <- apply(X = geno.import, MARGIN = 2, FUN = xa.convert, num.sample = nrow(geno.import)/2, maf.lim = 0.05)
  
  # Filter out any columns that show the error value 3
  xa.mat <- xa.mat[, xa.mat[1, ] != error.val] # A complicated way: xa.mat <- xa.mat[, apply(xa.mat == 3, 2, sum) == 0]
  
  # Rescale xa.mat so that it takes the value of -1 or 1
  xd.mat <- 1 - 2 * abs(xa.mat)
  
  # Return a list containing both the xa and xd dummy variable coding for our genotype data
  return(list(xa.mat, xd.mat))
}

# Additional part: test the function xa.convert
# Example 1
# geno.import <- data.frame(genotype1 = c("T", "A", "A", "A"), genotype2 = c("T", "A", "A", "A"), genotype3 = c("T", "T", "T", "T"))
# xa.mat  <- apply(X = geno.import, MARGIN = 2, FUN = xa.convert, num.sample = nrow(geno.import)/2, maf.lim = 0.05)
# xa.mat
#      genotype1 genotype2 genotype3
# [1,]         0         0         3
# [2,]        -1        -1         3
# xa.mat <- xa.mat[, apply(xa.mat == 3, 2, sum) == 0]
# xa.mat
#      genotype1 genotype2
# [1,]         0         0
# [2,]        -1        -1
# Example 2
# geno.import <- data.frame(genotype1 = sample(c("A", "C", "T", "G"), size = 100, replace = T), genotype2 = sample(c("A", "C", "T", "G"), size = 100, replace = T), genotype3 = c("A", rep("C", times = 99)))
# xa.mat  <- apply(X = geno.import, MARGIN = 2, FUN = xa.convert, num.sample = nrow(geno.import)/2, maf.lim = 0.05)
# xa.mat <- xa.mat[, apply(xa.mat == 3, 2, sum) == 0]

xa <- geno.recode(geno.import = geno.data, maf.lim = 0.05, error.val = 3)[[1]]
xd <- geno.recode(geno.import = geno.data, maf.lim = 0.05, error.val = 3)[[2]]
