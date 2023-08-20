#####Estimate the coverage#####
# C: coverage; G: haploid genome length (bp); L: read length (bp); N: number of reads

G <- 3.05*(10^9) # Human genome size
L <- 5000
N <- 189*(10^6)

C <- (L*N)/G
C
# Interpret: each base in the human genome will be sequenced between 309 and 310 times on average

# Use the Poisson distribution to estimate the probability of a base being sequenced 3 times or less
## Manual calculation
P_y0 <- (C^(0)*exp(-C))/factorial(0) # 2.753635e-135: precentage of the bases in the genome not yet sequenced
P_y1 <- (C^(1)*exp(-C))/factorial(1)
P_y2 <- (C^(2)*exp(-C))/factorial(2)
P_y3 <- (C^(3)*exp(-C))/factorial(3)
P_final <- P_y3+P_y2+P_y1+P_y0
P_final # 1.378361e-128
## Use the default function in R
ppois(q = 3, lambda = C) # 1.378361e-128

# Calculate the number of bases not yet sequenced = total gap length
P_y0 * G # 8.398588e-126

# Calculate the number of gaps
P_y0 * N # 5.204371e-127

# Extended calculation: if there is no possibility of a base being sequenced 3 times or less, what coverage is expected?
# Function to calculate P_final for a given C
calculate_P_final <- function(C) {
  P_y0 <- (C^(0)*exp(-C))/factorial(0)
  P_y1 <- (C^(1)*exp(-C))/factorial(1)
  P_y2 <- (C^(2)*exp(-C))/factorial(2)
  P_y3 <- (C^(3)*exp(-C))/factorial(3)
  P_final <- P_y3 + P_y2 + P_y1 + P_y0
  return(P_final)
}

# Binary search to find the minimum C for P_final = 0
lower_bound <- 0
upper_bound <- 1000 # You can adjust this upper bound as needed

while (upper_bound - lower_bound > 1e-10) {
  mid_point <- (lower_bound + upper_bound) / 2
  P_final <- calculate_P_final(mid_point)
  
  if (P_final == 0) {
    return(mid_point)
  } else if (P_final < 0) {
    upper_bound <- mid_point
  } else {
    lower_bound <- mid_point
  }
}

cat("Minimum C for P_final =", mid_point, "\n")
# Minimum C for P_final = 750 


# If there is 0.5% probability of a base being sequenced 3 times or less, what coverage is expected?
lower_bound <- 0
upper_bound <- 1000 
while (upper_bound - lower_bound > 1e-10) {
  mid_point <- (lower_bound + upper_bound) / 2
  P_final <- calculate_P_final(mid_point)
  
  if (P_final == 0.005) {
    return(mid_point)
  } else if (P_final < 0.005) {
    upper_bound <- mid_point
  } else {
    lower_bound <- mid_point
  }
}

cat("Minimum C for P_final =", mid_point, "\n")
# Minimum C for P_final = 10.97748 


# If there is 15% probability of a base being sequenced 3 times or less, what coverage is expected?
lower_bound <- 0
upper_bound <- 1000 
while (upper_bound - lower_bound > 1e-10) {
  mid_point <- (lower_bound + upper_bound) / 2
  P_final <- calculate_P_final(mid_point)
  
  if (P_final == 0.15) {
    return(mid_point)
  } else if (P_final < 0.15) {
    upper_bound <- mid_point
  } else {
    lower_bound <- mid_point
  }
}

cat("Minimum C for P_final =", mid_point, "\n") 
# Minimum C for P_final = 6.013537 
