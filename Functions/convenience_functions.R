# This script contains convenience functions used in other scripts (including other functions)

# Numerator and denominator of bray curtis distance
bray1 <- function(x,y){
  num = sum(abs(x - y))
  den = sum(x,y)
  return(c(num=num, den=den))
}

# Logit
logit <- function(x){
  log(x/(1-x))
}

# Inverse logit
inv_logit <- function(x){
  exp(x) / (1 + exp(x))
} 

# Reparametrise mu / phi into shape1 (alpha) / shape2 (beta) of beta distribution
mu_phi_to_a_b <- function(mu, phi){
  c(a = mu*phi, 
    b = (1-mu)*phi)
}

# Vectorised version of rpois and rbinom functions
vrpois <- Vectorize(rpois)
vrbinom <- Vectorize(rbinom)

# scale between 0 and 1
scale01 <- function(x){
  (x - min(x)) / (max(x) - min(x))
}





