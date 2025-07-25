#' Obtain LOO-IC, WAIC or RMSE from \code{greta} Output
#'
#' @param observed Vector of the observed response variable
#' @param posterior MCMC posterior from a \code{greta} model
#' @param mean Name of the mean/location parameter, defaults to \code{mu}
#' @param scale Name of the scale/variance parameter, defaults to \code{sd}
#' @param family Distribution of the observational model. Currently, only
#' \code{normal}, \code{lognormal}, and \code{student} are implemented
#' @param method
#' \code{"loo"} to calculate the Leave-One-Out Information Criteria or
#' \code{"waic"} for the Widely Applicable Information Criteria in the \code{loo} package
#' @param moment_match Logical, whether to use moment_match for loo.
#' @param ... additional arguments passed to distribution functions
#'
#' @return LOO-IC or WAIC; see \code{?loo::loo} for more details.
#'
#' @examples
#' obs <- iris$Sepal.Length
#' loo_greta(obs, posterior_lognormal, mean = "mean", scale = "sd",
#'           family = "lognormal", method = "loo")
#' loo_greta(obs, posterior_normal, mean = "mean", scale = "sd",
#'           family = "normal", method = "loo")
#'
#' @import loo
#' @export
#' 

library(loo)
loo_greta <- function(observed,
                      posterior,
                      chains = 4,
                      mean = "mu",
                      scale = NULL,
                      family = c("normal", "student", 'beta.binomial', 'binomial', 'beta'),
                      method = c("loo", "waic", 'RMSE'),
                      trials = NULL,
                      moment_match = FALSE,
                      ...) {
  # check that arguments are consistent with the implemented options
  family <- match.arg(family)
  method <- match.arg(method)
  
  posterior <- do.call(rbind, posterior)
  # extract the relative variables locally
  mu_cols <-  grepl(paste0(mean), colnames(posterior))
  mu <- posterior[, mu_cols]
  
  if(!is.null(scale)){
  scale_cols <-  grepl(paste0(scale), colnames(posterior))
  sd <- posterior[,scale_cols]
  } 
  
  nsim <- nrow(mu)
  
  # generate a matrix of log-likelihood values for each posterior predicted value
  if (family == "normal") {
    LL_mat <- lapply(
      seq_len(nsim),
      function(i){
        dnorm(observed, mu[i,], sd[i], log = TRUE, ...)
      }
    )
  } else if (family == "student") {
    LL_mat <- lapply(
      seq_len(nsim),
      function(i){
        dlst(x = observed, df = 2, mu = mu[i,], sigma = sd[i], log = TRUE, ...)
      }
    )
  } else if (family == 'beta.binomial') {
    LL_mat <- lapply(
      seq_len(nsim),
      function(i){
        shape1 = mu[i,]*sd[i]
        shape2 = (1 - mu[i,])*sd[i]
        VGAM::dbetabinom.ab(x = observed, size = trials, shape1, shape2, log = TRUE)
      }
    )
  } 
  else if (family == 'binomial') {
    LL_mat <- lapply(
    seq_len(nsim),
    function(i){
      stats::dbinom(x = observed, size = trials, mu[i,], log = TRUE)
    })
  }else if (family == 'beta') {
    LL_mat <- lapply(
      seq_len(nsim),
      function(i){
        shape1 = mu[i,]*sd[i]
        shape2 = (1 - mu[i,])*sd[i]
        dbeta(x = observed, shape1 = shape1, shape2 = shape2, log = TRUE)
      }
    )
  }
  
  LL_mat <- do.call(rbind, LL_mat)
  
  # compute the relative inefficiencies of the posterior predicted values
  rel_eff <- loo::relative_eff(exp(LL_mat), chain_id = rep(1:chains, each = nsim/chains))
  
  # estimate the LOO-IC, WAIC, etc
  out <- switch(method,
                loo  = loo::loo(LL_mat, r_eff = rel_eff, moment_match = moment_match),
                waic = loo::waic(LL_mat, r_eff = rel_eff),
                RMSE = loo::loo_predictive_metric(x = mu, y = observed, LL_mat, r_eff = rel_eff, metric = 'rmse')
  )
  
  return(out)
}


