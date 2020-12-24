# projection of the matrix Y onto Omega (observed entries).

Projection <- function(Y, Omega){
  pro <- Y
  pro[Omega == 0] <- 0
  return(pro)
}