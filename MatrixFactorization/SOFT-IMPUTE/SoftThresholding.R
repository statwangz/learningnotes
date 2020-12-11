# soft-thresholding

SoftThresholding <- function(W, lambda){
  s <- svd(W) # SVD
  d_lambda <- s$d
  d_lambda[d_lambda < lambda] <- 0
  S <- s$u %*% diag(d_lambda) %*% t(s$v)
  return(S)
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