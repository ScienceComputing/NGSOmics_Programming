# Calculate the error term (gamma inverse)
gamma_inv_calc <- function(X_mx, beta_t){
  ## Initialize gamma
  K <- X_mx %*% beta_t
  gamma_inv <- exp(K)/(1 + exp(K))
  return(gamma_inv)
}

# Calculate the variance of our error
W_calc <- function(gamma_inv){
  W <- diag(as.vector(gamma_inv * (1 - gamma_inv)))
  return(W)
}

# Calculate beta estimates given the error term and variance of error 
beta_update <- function(X_mx, W, Y, gamma_inv, beta){
  beta_up <- beta + ginv(t(X_mx) %*% W %*% X_mx) %*% t(X_mx) %*% (Y - gamma_inv)
  return(beta_up)
}

# Calculate re-weighting model deviance and log-likelihood of final estimates
dev_calc <- function(Y, gamma_inv){
  deviance <- 2*( sum(Y[Y==1]*log(Y[Y==1]/gamma_inv[Y==1])) + 
                  sum((1-Y[Y==0])*log((1-Y[Y==0])/(1-gamma_inv[Y==0]))) )  
  return(deviance)
}

loglik_calc <- function(Y, gamma_inv){
  loglik <- sum(Y*log(gamma_inv)+(1-Y)*log(1-gamma_inv))
  return(loglik)
}

# Wrap up everything into one function that returns the betas and log-likelihoods
logistic.IRLS.covars <- function(Xa, Xd, Xz, Y = Y, 
                                 beta.initial.vec = c(0, 0, 0, 0), 
                                 d.stop.th = 1e-6, it.max = 100) {
  
  ## Create initial values
  X_mx <- cbind(rep(1, length(Xa)), Xa, Xd, Xz)
  beta_t <- beta.initial.vec
  dt <- 0
  gamma_inv <- gamma_inv_calc(X_mx, beta_t)
  
  ## Start the optimization loop
  for(i in 1:it.max) {
    
    ## Store the previous deviance
    dpt1 <- dt 
    
    ## Calculate W (variance of errors)
    W <- W_calc(gamma_inv)
    
    ## Update beta
    beta_t <- beta_update(X_mx, W, Y, gamma_inv, beta_t)
    
    ## Update gamma_inv (error term)
    gamma_inv <- gamma_inv_calc(X_mx, beta_t)
    
    ## Calculate new deviance
    dt <- dev_calc(Y, gamma_inv)
    absD <- abs(dt - dpt1)
    
    ## Check if deviance is smaller than the threshold
    if(absD < d.stop.th) {
      logl <- loglik_calc(Y, gamma_inv) ## Log likelihood goes here
      return(list(beta_t, logl)) ## Return a list that has beta.t and logl saved
    }   
  }
  
  ## If the algorithm does not converge when the iteration step reaches the max
  return(list(beta_t = c(NA, NA, NA, NA), logl = NA)) ## Return NA values 
}

# Create one function that returns the p-values
logistic.IRLS.pval.covars <- function(Xa, Xd, Xz, Y = Y, 
                                      beta.initial.vec = c(0, 0, 0, 0),
                                      d.stop.th = 1e-6, it.max = 100) {
  
  h1 <- logistic.IRLS.covars(Xa, Xd, Xz, Y = Y)
  
  h0 <- logistic.IRLS.covars(Xa = rep(0, nrow(Y)), Xd = rep(0, nrow(Y)), Xz, Y = Y)
  
  LRT <- 2*h1[[2]] - 2*h0[[2]] ## Likelihood ratio test statistic
  pval <- pchisq(LRT, 2, lower.tail = F)
  return(pval)
}

# Run the p-value calculator on all SNPs
xz.mat <- ancestry.data$ancestry |> as.matrix()
allPvals.covars <- mclapply(1:ncol(xa.mat), 
                            function(SNP.i) {
                              logistic.IRLS.pval.covars(Xa = xa.mat[, SNP.i],
                                                        Xd = xd.mat[, SNP.i],
                                                        Xz = xz.mat,
                                                        Y = Y) },
                            mc.cores = 4) %>% 
  do.call(rbind, .)
