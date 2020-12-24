# Matrix Completion and Low-Rank SVD via Fast Alternating Least Squares
# Algorithm 4.1: softImpute-ALS

# X, observed data, incomplete matrix, unobserved entries are represented as NA
# Z, solution

softImpute_ALS <- function(X, lambda, converge = 1e-7, MAX = 100){
  
  # indices of observed entries
  Omega <- X
  Omega[!is.na(Omega)] <- 1
  Omega[is.na(Omega)] <- 0
  
  ProX <- Projection(X, Omega)
  
  r <- sum(svd(ProX)$d > lambda) + 1
  if(r < 2){
    r <- 2
  }
  m <- nrow(X)
  n <- ncol(X)
  
  # initialize
  # A <- matrix(rnorm(m * r), nrow = m, ncol = r)
  # B <- matrix(rnorm(n * r), nrow = n, ncol = r)
  A <- diag(1, nrow = m, ncol = r)
  B <- diag(1, nrow = n, ncol = r)
    
  t <- 1 # counter
  
  repeat{
    
    X_star <- ProX + Complementary(A %*% t(B), Omega)
    A_new <- X_star %*% B %*% solve(t(B) %*% B + lambda * diag(r))
    X_star <- ProX + Complementary(A_new %*% t(B), Omega)
    B_new <- t(X_star) %*% A_new %*% solve(t(A_new) %*% A_new + lambda * diag(r))
    
    if(Frobenius(A_new - A) / Frobenius(A) + Frobenius(B_new - B) / Frobenius(B) < converge || t == MAX){
      break
    }
    
    A <- A_new
    B <- B_new
    t <- t + 1
    
  }
  
  rank1 <- sum(svd(A%*%t(B))$d > lambda)
  rank2 <- rankMatrix(A%*%t(B))[1]
  return(list(A = A, B = B, rank1 = rank1, rank2 = rank2, r = r))
  
}