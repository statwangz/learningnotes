# compute the smallest value for lambda such that SoftImpute(X,lambda) returns the zero solution

Lambda0 <- function(X){
  X[is.na(X)] <- 0
  lambda <- svd(X)$d[1]
  if(lambda > 0){
    return(lambda)
  }else{
    return(0)
  }
}