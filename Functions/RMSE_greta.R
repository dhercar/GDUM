
RMSE_greta <- function(observed, 
                       posterior, 
                       mean = 'mu',
                       trials = NULL,
                       summary = FALSE,
                       summary_quantiles = c(0.025, 0.5, 0.975)){
  
  posterior <- do.call(rbind, posterior)
  mu_cols <-  grepl(paste0(mean), colnames(posterior))
  estimated <- posterior[, mu_cols]
  
  if (!is.null(trials)) {
    observed <- observed/trials
  }
  
  residual <- sweep(estimated, 2, observed) * -1
  RMSE <- apply(residual, 1, function(x) sqrt(sum((abs(x)^2))/length(observed)))
  
  if(summary == TRUE){
    return(quantile(RMSE, probs = summary_quantiles))
  }else {
    return(RMSE)
  }
  
}

