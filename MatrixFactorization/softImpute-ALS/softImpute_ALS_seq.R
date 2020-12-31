# Matrix Completion and Low-Rank SVD via Fast Alternating Least Squares
# Algorithm 4.1: softImpute-ALS

# X, observed data, incomplete matrix, unobserved entries are represented as NA
# Z, solution

softImpute_ALS_seq <- function(X, lambda, K = 1000, converge = 1e-5, MAX = 2000){
  
  # indices of observed entries
  Omega <- X
  Omega[!is.na(Omega)] <- 1
  Omega[is.na(Omega)] <- 0
  
  ProX <- Projection(X, Omega)
  d <- svd(ProX)$d
  
  m <- nrow(X)
  n <- ncol(X)
  
  # lambda <- seq(from = lambda, to = 0, length.out = K)
  lambda <- exp(seq(from = log(lambda), to = 0, length.out = K)) # log-space
  Z_hat <- array(dim = c(K, m, n))
  Z_hat2 <- array(dim = c(K, m, n))
  Z_hat3 <- array(dim = c(K, m, n))
  A <- list()
  B <- list()
  NuclearNorm <- vector()
  Rank <- vector()
  r <- vector()
  
  k <- 1 # counter
  
  while (k <= K) {
    
    r_k <- sum(d > lambda[k]) + 1
    if(r_k < 2){
      r_k <- 2
    }
    
    # initialize, warm starts
    if(k > 1 && r_k == r[k - 1]){
      A_old <- A[[k - 1]]
      B_old <- B[[k - 1]]
    }else{
      A_old <- diag(1, nrow = m, ncol = r_k)
      B_old <- diag(1, nrow = n, ncol = r_k)
    }
    
    t <- 1 # counter
    
    repeat{
      
      X_star <- ProX + Complementary(A_old %*% t(B_old), Omega)
      if(rankMatrix(t(B_old) %*% B_old + lambda[k] * diag(r_k)) < r_k){
        A_new <- A_old
        B_new <- B_old
        break
      }
      A_new <- X_star %*% B_old %*% solve(t(B_old) %*% B_old + lambda[k] * diag(r_k))
      X_star <- ProX + Complementary(A_new %*% t(B_old), Omega)
      if(rankMatrix(t(A_new) %*% A_new + lambda[k] * diag(r_k)) < r_k){
        B_new <- B_old
        break
      }
      B_new <- t(X_star) %*% A_new %*% solve(t(A_new) %*% A_new + lambda[k] * diag(r_k))
      
      if(Frobenius(A_new - A_old) / Frobenius(A_old) + Frobenius(B_new - B_old) / Frobenius(B_old) < converge || t == MAX){
        break
      }
      
      A_old <- A_new
      B_old <- B_new
      t <- t + 1
      
    }
    
    Z <- A_new %*% t(B_new)
    A <- c(A, list(A_new))
    B <- c(B, list(B_new))
    svd_k <- svd(Z)
    d_k <- svd_k$d
    
    q <- sum(d_k > lambda[k])
    if(q > 1){
      Z_hat[k, , ] <- svd_k$u[ , 1 : q] %*% diag(d_k[1 : q]) %*% t(svd_k$v)[1 : q, ]
      Z_hat2[k, , ] <- svd_k$u[ , 1 : q] %*% diag(d_k[1 : q] - lambda[k]) %*% t(svd_k$v)[1 : q, ]
    }else if(q == 1){
      Z_hat[k, , ] <- d_k[1] * svd_k$u[ , 1] %*% t(svd_k$v)[1, ]
      Z_hat2[k, , ] <- (d_k[1] - lambda[k]) * svd_k$u[ , 1] %*% t(svd_k$v)[1, ]
    }else{
      Z_hat[k, , ] <- d_k[1] * svd_k$u[ , 1] %*% t(svd_k$v)[1, ]
      Z_hat2[k, , ] <- Z_hat[k, , ]
    }
    Z_hat3[k, , ] <- Z
    
    NuclearNorm <- c(NuclearNorm, sum(d_k[d_k > lambda[k]]))
    Rank <- c(Rank, sum(d_k > lambda[k]))
    r <- c(r, r_k)
    k <- k + 1
    
  }
  
  return(list(A = A, B = B, Z_hat = Z_hat, Z_hat2 = Z_hat2, Z_hat3 = Z_hat3, NuclearNorm = NuclearNorm, Rank = Rank, lambda = lambda, r = r))
  
}