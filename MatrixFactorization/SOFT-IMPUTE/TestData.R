# Generate test data, a matrix with unobserved entries
# m, nrows
# n, ncols
# r, rank
# p, the proportion of unobserved entries to total
# SNR, the ratio of the standard deviation of the entries to the standard deviation of the noise

TestData <- function(m = 100, n = 100, r = 10, p = 0.5, SNR = 1){
  U <- matrix(rnorm(m*r), m, r)
  V <- matrix(rnorm(n*r), n, r)
  Z <- U %*% t(V) + matrix(rnorm(m*n, sd = 1/SNR), m, n) # complete matrix with noise
  miss <- sample(1:(m*n), size = m*n*p, replace = FALSE)
  X <- Z
  X[miss] <- NA # observed entries
  return(list(U = U, V = V, Z = Z, X = X))
}