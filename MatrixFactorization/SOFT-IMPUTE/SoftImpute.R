# Spectral Regularization Algorithms for Learning Large Incomplete Matrices
# Algorithm 1: SOFT-IMPUTE

# X, observed data, incomplete matrix, unobserved entries are represented as NA
# Z, solution

SoftImpute <- function(X, lambda0, K = 100, converge = 0.01, MAX = 20){
  
  # indices of observed entries
  Omega <- X
  Omega[!is.na(Omega)] <- 1
  Omega[is.na(Omega)] <- 0
  
  ProX <- Projection(X, Omega)
  
  m <- nrow(X)
  n <- ncol(X)
  lambda <- seq(from = lambda0, to = 0, length.out = K)
  Z <- array(dim = c(K, m, n))
  
  s <- 1 # counter
  
  while(s <= K){
    
    Z_old <- matrix(0, nrow = m, ncol = n) # initialize
    t <- 1 # counter
    
    repeat{
      
      Z_new <- SoftThresholding(ProX + Complementary(Z_old, Omega), lambda[s])
      
      if((Frobenius(Z_new - Z_old) / Frobenius(Z_new) < converge) | t == MAX){
        break
      }
      
      Z_old <- Z_new
      t <- t + 1
      
    }
    
    Z[s, , ] <- Z_new
    s <- s + 1
    
  }
  
  return(Z)
  
}