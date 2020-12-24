# Spectral Regularization Algorithms for Learning Large Incomplete Matrices
# Algorithm 1: SOFT-IMPUTE

# X, observed data, incomplete matrix, unobserved entries are represented as NA
# Z, solution

softImpute <- function(X, lambda, converge = 1e-5, MAX = 100){
  
  # indices of observed entries
  Omega <- X
  Omega[!is.na(Omega)] <- 1
  Omega[is.na(Omega)] <- 0
  
  ProX <- Projection(X, Omega)
  
  m <- nrow(X)
  n <- ncol(X)
  # initialize
  Z <- matrix(0, nrow = m, ncol = n)
  
  t <- 1 # counter
  
  repeat{
    
    Z_new <- SoftThresholding(ProX + Complementary(Z, Omega), lambda)
    
    if(Frobenius(Z_new - Z) / Frobenius(Z_new) < converge || t == MAX){
      break
    }
    
    Z <- Z_new
    t <- t + 1
    
  }
  
  Z <- Z_new
  NuclearNorm <- sum(svd(Z)$d)
  Rank <- rankMatrix(Z)[1]
  
  return(list(Z_hat = Z, NuclearNorm = NuclearNorm, Rank = Rank, lambda = lambda))
  
}