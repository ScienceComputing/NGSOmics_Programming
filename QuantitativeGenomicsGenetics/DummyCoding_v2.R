# Goal: covert the ATCGs into dummy numbers under either additive or dominant genetic model

# Input genotype file: each column represents a specific SNP (column 1 = allele 1 of the genotype 1, column 2 = allele 2 of the genotype 1).
# Each consecutive pair of columns represent two alleles of a genotype across individuals.
# Each row represents all of individual 1’s SNPs.

# Takes in the genotype data and is applied two columns at a time (one locus at a time)   
genotype_coder <- function(geno_import, maf_limit, error_value = 3){
  
  # This line stacks the two columns of the current SNP into one double long vector
  geno_input = mapply(FUN = c, geno_import[,seq(1,ncol(geno_import),2)], geno_import[,seq(2,ncol(geno_import),2)])
  # 把(1,2)列合并成一个长列, (3,4)合并成一个长列, (5,6)合并成一个长列。。。最终效果：geno_input的行数是geno_import的2倍，geno_input的列数是geno_important的1/2.
  
  # Inputs: geno_col is the vector of alleles from the last line, numSamples=number of rows, maf_limit
  xa_converter <- function(geno_col, numSamples, maf_limit){
    # What alleles are present at this locus?
    geno_count <- table(geno_col) 
    # If the MAF is less than the pre-set threshold OR this locus has less than 2 alleles (e.g., this site is mono-allelic (has only 1 allele))
    if(min(geno_count)/length(geno_col) <= maf_limit | length(geno_count) < 2){
      # Return a vector of 3's, to be filtered out later.
      return(rep(error_value, numSamples)) 
    }
    # Otherwise find our minor allele
    minor_allele <- names(geno_count[geno_count == min(geno_count)])
      # geno_col[1:numSamples] refers to the 1st alleles of all genotypes
      # geno_col[(numSamples+1):length(geno_col)] refers to the 2nd alleles of all genotype.
      # the conditional results in a boolean (T/F) but the vector operation (+) converts it to (1/0)
      # e.g., an individual with TT where T is the minor allele would have True + True = 2
    xa <- (geno_col[1:numSamples]==minor_allele) + (geno_col[(numSamples+1):length(geno_col)]==minor_allele)
    # We want our dummy variable coding centered on 1 so take 0,1,2 and minus 1
    xa <- xa-1
    return(xa)
  }
  
  # Take our input, apply the xa_converter function to the input by column (2), with numSamples=nrow(data), and maf_limit=0.05
  xa_mat  <- apply(geno_input, 2, xa_converter, numSamples = nrow(geno_import), maf_limit = 0.05)
  # Filter out any columns that have the set error value 3. 
  xa_mat <- xa_mat[,xa_mat[1,]!=error_value]
  # To get our Xd dominance dummy variable coding use the algebra below.
  xd_mat <- 1 - 2*abs(xa_mat)
  
  # Return both the Xa and Xd dummy variable coding for our input matrix/data
  return(list(xa_mat,xd_mat))
}
