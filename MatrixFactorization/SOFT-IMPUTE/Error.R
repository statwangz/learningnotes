# compute the test error and the training error
# U, original matrix
# V, original matrix
# Z, complete matirx with noise
# Z_hat, Approximation of the complete matrix
# Omega, indices of observed entries

Error <- function(U, V, Z, Z_hat, Omega){
  test <- Frobenius(Complementary(U %*% t(V) - Z_hat, Omega)) / Frobenius(Complementary(U %*% t(V), Omega))
  training <- Frobenius(Projection(Z - Z_hat, Omega)) / Frobenius(Projection(Z, Omega))
  return(list(test.error = test, training.error = training))
}