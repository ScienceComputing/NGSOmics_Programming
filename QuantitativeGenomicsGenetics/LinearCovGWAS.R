# Attach required libraries
library(MASS)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(data.table)

# Read-in genotype data
gen_import <- read.csv(file = './genotype_data.csv', 
                       header = TRUE, 
                       stringsAsFactors = FALSE,
                       row.names = 1, colClasses = "character")

# Read-in phenotype data
sim_pheno_mx <- read.csv(file = './phenotype_data.csv', 
                         header = TRUE, row.names = 1)


# Read-in covariate data
xc_mat = read.csv(file = 'covar_data.csv', row.names = 1)


# Convert genotypes to Xa & Xd
genotype_coder <- function(gen_import, maf_limit){
  gen_input <- mapply(c, gen_import[, seq(1, ncol(gen_import), 2)], gen_import[, seq(2, ncol(gen_import), 2)])
  # (1,2)列合并成一个长列，(3,4)列合并成一个长列...
  
  xa_converter <- function(gen_col, numSamples, maf_limit){
    geno_count <- table(gen_col)
    if(min(geno_count)/length(gen_col) <= maf_limit){
      return(rep(3, numSamples))
    }
    minor_allele <- names(geno_count[geno_count == min(geno_count)])
    xa <- (gen_col[1:numSamples] == minor_allele) + (gen_col[(numSamples+1):length(gen_col)] == minor_allele)
    xa <- xa-1
    return(xa)
  }
  
  xa_mat  <- apply(gen_input, 2, xa_converter, nrow(gen_import), 0.05)
  xa_mat <- xa_mat[, xa_mat[1, ] != 3]
  xd_mat <- 1 - 2*abs(xa_mat)
  
  return(list(xa_mat, xd_mat))
}

codes <- genotype_coder(gen_import, 0)
xa_mat <- codes[[1]]
xd_mat <- codes[[2]]


# Define the function to run GWAS under linear regression framework and return p-values
pval_calculator_w_covars <- function(pheno_input, xa_input, xd_input, xz_input){
  n_samples <- length(xa_input) 
  X_mx <- cbind(rep(1, length(xa_input)), xa_input, xd_input, xz_input) 
  
  MLE_beta <- ginv(t(X_mx) %*% X_mx) %*% t(X_mx) %*% pheno_input 
  
  x_h0 =  cbind(rep(1, length(xa_input)), xz_input) 
  MLE_h0 = ginv(t(x_h0) %*% x_h0) %*% t(x_h0) %*% pheno_input 
  y_hat_0 = x_h0 %*% MLE_h0 
  y_hat_1 = X_mx%*% MLE_beta 
  
  SSE_theta_0 = sum((pheno_input - y_hat_0)^2) 
  SSE_theta_1 = sum((pheno_input - y_hat_1)^2) 
  
  df_M <- 2
  df_E <- n_samples - 4 
  
  numerator <- (SSE_theta_0 - SSE_theta_1) / df_M 
  denom <- SSE_theta_1 / df_E
  Fstatistic <- numerator / denom
  
  pval <- pf(Fstatistic, df_M, df_E, lower.tail = FALSE) 
  return(pval)
}


# Run the function and organize the results
results <- lapply(1:ncol(xa_mat), function(column.counter){
  data.table(pval_calculator_w_covars(pheno_input = sim_pheno_mx[, 2],
                                      xa_input = xa_mat[, column.counter],
                                      xd_input = xd_mat[, column.counter],
                                      xz_input = xc_mat[, 1]))
}) %>% 
  rbindlist() %>% 
  mutate(p = V1, index = 1:ncol(xa_mat))


# Make the Manhattan Plot
my.alpha = 0.05/ncol(xa_mat)
man <- ggplot(results, aes(x = index, y = -log10(p))) +
  geom_point() + 
  geom_hline(yintercept = -log10(my.alpha), color = 'red', lty = 2) +
  labs(x = 'Index', y = expression(-log[10]~p), 
       title = 'GWAS Manhattan Plot', subtitle = 'Covariates Included')


# Make the QQ plot
observed_pvals = sort(results$p)
expected_pvals = qunif(seq(0, 1, length.out = length(observed_pvals) + 2), min = 0, max = 1)  
expected_pvals = expected_pvals[expected_pvals != 0 & expected_pvals != 1]  

p_df = data.frame(observed = -log10(observed_pvals),
                  expected = -log10(expected_pvals))

qq <- ggplot(p_df, aes(x = expected, y = observed)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  labs(x = '-log10 Expected p-val',
       y = '-log10 Observed p-val',
       title = 'GWAS QQ plot',
       subtitle = 'Covariates Included')
grid.arrange(man, qq, ncol = 2)
