# Load genotypes data
geno.data <- read.table("genotypes.txt", 
                        sep = ",", header = F)
                        
# Convert SNP coded as 0, 1, 2 to xa and xd
## Calculate xa and xd
xa <- geno.data - 1
xa.mat <- data.matrix(xa)
xd.mat <- 1 - 2*abs(xa.mat)
# xa.mat[1:6, 1:6]
# table(xa.mat[, 1])
# xd.mat[1:6, 1:6]
# table(xd.mat[, 1])

## Save xa.mat and xd.mat for fast loading next time
save(xa.mat, xd.mat, file = "xa.xd.RData")
