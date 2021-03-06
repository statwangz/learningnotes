# compute the test error and the training error
# U, original matrix
# V, original matrix
# Z, complete matirx with noise
# Z_hat, Approximation of the complete matrix
# Omega, indices of observed entries

Error <- function(U, V, Z, Z_hat, X){
  Omega <- X
  Omega[!is.na(Omega)] <- 1
  Omega[is.na(Omega)] <- 0
  test <- Frobenius(Complementary(U %*% t(V) - Z_hat, Omega)) / Frobenius(Complementary(U %*% t(V), Omega))
  training <- Frobenius(Projection(Z - Z_hat, Omega)) / Frobenius(Projection(Z, Omega))
  return(c(test, training))
}