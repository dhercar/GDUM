SS_calc <- function(x, LCBD = FALSE){
  
  stopifnot(is.matrix(x), isSymmetric(x))
  
  n <- nrow(x)
  H <- diag(n) - matrix(1, n, n)/n 
  A <- -0.5 * H %*% x %*% H
  SS <- diag(A)
  
  if (LCBD) {
   return(SS/sum(SS)) 
  } else {
   return(SS)
  }
}
