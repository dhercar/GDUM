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


#' Bayesian Information Criterion (BIC) Calculation
#' 
#' This function calculates the Bayesian Information Criterion (BIC) for a given model.
#' 
#' @param log_likelihood The log-likelihood of the model.
#' @param n The number of samples.
#' @param k The number of parameters estimated in the model.
#' 
#' @return The BIC value for the model.
#' #' @details The BIC is calculated using the formula:
#' \deqn{BIC = -2 * log_likelihood + k * log(n)}
#' #' where:
#' \itemize{
#'  \item log_likelihood: The log-likelihood of the model.
#' #  \item n: The number of samples.
#' \item k: The number of parameters estimated in the model.
#' #' }
#' #' The BIC is a criterion for model selection among a finite set of models. It is based on the likelihood function and includes a penalty term for the number of parameters in the model.
#' #' A lower BIC value indicates a better model fit, taking into account the complexity of the model.
#' #'
#' @examples
#' # Example usage of BICcalc function
#' log_likelihood <- -150.5
#' n <- 100
#' k <- 5
#' BIC_value <- BICcalc(log_likelihood, n, k)
#' print(BIC_value)
#' 
#'

BIC2 <- function(log_likelihood, # model log-likelihood, 
                 n,              # n samples
                 k) {             # k parameters estimated
  
  BIC <- -2 * log_likelihood + k * log(n)
  return(BIC)
}

