# bayesian r2

r2_greta <- function(observed, 
                     estimated, 
                     latent = FALSE, # Set to TRUE if r2 is computed on u. If set to TRUE, 'observed' should corresponds to the estimated response and 'estimated' to the estimated response without the site-level re 
                     mean = 'mu',
                     trials = NULL,
                     summary = FALSE,
                     summary_quantiles = c(0.025, 0.5, 0.975)) {
  
  estimated <- do.call(rbind, estimated)
  mu_cols <-  grepl(paste0(mean), colnames(estimated))
  estimated <- estimated[, mu_cols]
  
  if (!is.null(trials)) {
   estimated <- sweep(estimated, 2, trials, FUN = '*')
  }
  
  # calculate residual 
  if (latent) {
    residual = estimated - observed
  } else {
   residual <- sweep(estimated, 2, observed) * -1
  }

  # Variances
  var_r <- apply(residual, 1, var) # residual
  var_m <- apply(estimated, 1, var) # estimated
  
  # r2
  r2 <- var_m / (var_m + var_r)
  
  if(summary == TRUE){
    return(quantile(r2, probs = summary_quantiles))
  }else {
    return(r2)
  }
}


