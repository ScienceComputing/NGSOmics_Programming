# First approach:
# We can specify the column classes as characters...
geno_import <- read.csv("./genotype_data.csv", 
                        header = T, 
                        stringsAsFactors = F,
                        row.names = 1, 
                        colClasses = "character") # ...Here!

# To avoid T being interpreted as TRUE, we can project the characters to lower case
geno_import <- apply(geno_import, 2, tolower)
head(geno_import[, 253:254])


# Second approach:
# We could change the data type/value by column
geno_import <- read.csv("./genotype_data.csv", 
                        header = TRUE, 
                        stringsAsFactors = F,
                        row.names = 1)

for(i in 1:ncol(geno_import)) {
  geno_import[, i] <- as.character(geno_import[, i])
  geno_import[geno_import[, i] == "TRUE", i] = "T" 
  # ALiu: geno_import[,i] -> march through each element in the ith column; return a logical vector which has the length of the ith column
}
head(geno_import[, 253:254])
