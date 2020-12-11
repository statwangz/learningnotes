# Frobenius norm of matrix

Frobenius <- function(A){
  f <- sum((as.vector(A))^2)
  return(f)
}