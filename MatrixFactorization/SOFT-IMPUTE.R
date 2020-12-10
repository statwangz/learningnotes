# Spectral Regularization Algorithms for Learning Large Incomplete Matrices
# Algorithm 1: SOFT-IMPUTE

# X, observed data, incomplete matrix, unobserved entries are represented as NA
# Z, solution

# projection of the matrix Y onto omega (observed entries).
Projection <- function(Y, omega){
  pro <- Y
  pro[omega == 0] <- 0
  return(pro)
}

# complementary projection
Complementary <- function(Y, omega){
  com <- Y
  com[omega == 1] <- 0
  return(com)
}

# SVD <- function(A){
#   light <- eigen(A %*% t(A))
#   right <- eigen(t(A) %*% A)
#   U <- light$vectors
#   V <- right$vectors
#   if(nrow(A) <= ncol(A)){
#     sigma <- sqrt(light$values)
#   }else{
#     sigma <- sqrt(right$values)
#   }
#   return(list(sigma = sigma, U = U, V = V))
# }

# soft-thresholding
SoftThresholding <- function(W, lamda){
  s <- svd(W)
  d_lamda <- s$d
  d_lamda[d_lamda < lamda] <- 0
  S <- s$u %*% diag(d_lamda) %*% t(s$v)
  return(S)
}

# Frobenius norm of matrix
Frobenius <- function(A){
  a <- as.vector(A)
  f <- sum(a^2)
  return(f)
}

SoftImpute <- function(X, lamda, converge){
  
  # indices of observed entries
  omega <- X
  omega[!is.na(omega)] <- 1
  omega[is.na(omega)] <- 0
  
  m <- nrow(X)
  n <- ncol(X)
  
  K <- length(lamda)
  Z <- array(dim = c(K, m, n))
  
  s <- 1 # counter
  
  while(s <= K){
    
    Z_old <- matrix(0, nrow = m, ncol = n) # initialize
    
    repeat{
      
      Z_new <- SoftThresholding(Projection(X, omega) + Complementary(Z_old, omega), lamda[s])
      
      if(Frobenius(Z_new - Z_old) / Frobenius(Z_new) < converge){
        break
      }
      
      Z_old <- Z_new
      
    }
    
    Z[s, , ] <- Z_new
    s <- s + 1
    
  }
  
  return(Z)
  
}