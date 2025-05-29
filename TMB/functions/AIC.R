#' Calculate AICc
#' 
#' This function calculates the Akaike Information Criterion corrected for small sample sizes (AICc).
#' 
#' @param log_likelihood The log-likelihood of the model.
#' @param n The number of samples.
#' @param k The number of parameters estimated in the model.
#' 
#' @return A list containing the AIC and AICc values.
#' @export
#' @examples
#' # Example usage
#' log_likelihood <- -100
#' n <- 50
#' k <- 5
#' result <- AICc(log_likelihood, n, k)
#' print(result)
#' 

AICc <- function(log_likelihood, # model log-likelihood, 
                 n,              # n samples
                 k) {            # k parameters estimated
  
  AIC <- -2 * log_likelihood + 2 * k
  correct <- (2 * k * (k + 1)) / (n - k - 1)
  
  return(list(AIC = AIC, AICc = AIC + correct))
}

