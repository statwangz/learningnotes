# complementary projection

Complementary <- function(Y, omega){
  com <- Y
  com[omega == 1] <- 0
  return(com)
}