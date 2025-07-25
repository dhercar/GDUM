print_title <- function(text, width = getOption("width"), symb = '-') {
  # Full-width line of dashes
  separator_line <- strrep(symb, width)
  
  # Centered text
  padding <- floor((width - nchar(text)) / 2)
  if (padding < 0) padding <- 0
  centered_text <- paste0(strrep(" ", padding), text)
  
  # Print all
  cat(separator_line, "\n")
  cat(centered_text, "\n")
  cat(separator_line, "\n\n")
}

print_title2 <- function(text, width = getOption("width"), symb = '-') {
  total_padding <- width - nchar(text)
  if (total_padding < 0) total_padding <- 0
  
  left_padding <- floor(total_padding / 2)
  right_padding <- ceiling(total_padding / 2)
  
  decorated_text <- paste0(strrep(symb, left_padding), text, strrep(symb, right_padding))
  cat(decorated_text, "\n\n")
}

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


