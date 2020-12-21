# Spectral Regularization Algorithms for Learning Large Incomplete Matrices
# Algorithm 1: SOFT-IMPUTE

# X, observed data, incomplete matrix, unobserved entries are represented as NA
# Z, solution

SoftImpute <- function(X, lambda, K = 1000, converge = 1e-5, MAX = 100){
  
  # indices of observed entries
  Omega <- X
  Omega[!is.na(Omega)] <- 1
  Omega[is.na(Omega)] <- 0
  
  ProX <- Projection(X, Omega)
  
  m <- nrow(X)
  n <- ncol(X)
  
  lambda <- seq(from = lambda, to = 0, length.out = K)
  # lambda <- exp(seq(from = log(lambda), to = 0, length.out = K)) ## log-space
  Z <- array(dim = c(K, m, n))
  NuclearNorm <- vector()
  Rank <- vector()
  
  s <- 1 # counter
  
  while(s <= K){
    
    # initialize, warm starts
    if(s > 1){
      Z_old <- Z[s-1, , ]
    }else{
      Z_old <- matrix(0, nrow = m, ncol = n)
    }
    
    t <- 1 # counter
    
    repeat{
      
      Z_new <- SoftThresholding(ProX + Complementary(Z_old, Omega), lambda[s])
      
      if(Frobenius(Z_new - Z_old) / Frobenius(Z_new) < converge || t == MAX){
        break
      }
      
      Z_old <- Z_new
      t <- t + 1
      
    }
    
    Z[s, , ] <- Z_new
    NuclearNorm <- c(NuclearNorm, sum(svd(Z_new)$d))
    Rank <- c(Rank, rankMatrix(Z_new)[1])
    s <- s + 1
    
  }
  
  return(list(Z_hat = Z, NuclearNorm = NuclearNorm, Rank = Rank, lambda = lambda))
  
}