summary.gdmm <- function(obj,
                         print_result = TRUE) {
    
    # Get coefficient table with p-values
    coef_table <- summary(obj$sdrep, select = c('fixed'), p.value = TRUE)
    
    # Add significance stars
    p_values <- coef_table[, "Pr(>|z^2|)"]
    sig_stars <- ifelse(p_values < 0.001, "***",
                        ifelse(p_values < 0.01, "**",
                               ifelse(p_values < 0.05, "*",
                                      ifelse(p_values < 0.1, ".", " "))))
    
    # Create data.frame with lm-style column names and add significance stars
    result <- as.data.frame(coef_table)

    
    # Calculate AIC 
    m_AIC <- AICc(log_likelihood = -obj$opt$objective,
                  n = length(obj$Y_diss),
                  k = length(obj$opt$par))
    m_BIC <- BIC2(log_likelihood = -obj$opt$objective,
                  n = length(obj$Y_diss),
                  k = length(obj$opt$par))
    
    out <- list(call = obj$call,
                table = coef_table, 
                p_values = p_values,
                names_beta =  paste0('diss: ', colnames(obj$form_X$predictors)),
                names_lambda =  paste0('uniq: ', colnames(obj$form_W$predictors)),
                sig = sig_stars,
                AIC = m_AIC$AIC,
                AICc = m_AIC$AICc,
                BIC = m_BIC,
                logL = -obj$opt$objective)
    class(out) <- 'summary.gdmm'
    return(out)
}


summary.bbgdmm <- function(obj, 
                         quantiles = c(0.025, 0.5, 0.975),
                         print_result = TRUE,
                         null_value = 0) {
  
  # Basic statistics 
  samples <- obj$boot_samples[, (colnames(obj$boot_samples) != 'logLikelihood') & (substr(colnames(obj$boot_samples), 1,5) != 'u_re_')]
  logL <- mean(obj$boot_samples[,'logLikelihood'])
  
  # Calculate AIC 
  m_AIC <- AICc(log_likelihood = logL,
                n = length(obj$Y_diss),
                k = ncol(samples))
  m_BIC <- BIC2(log_likelihood = logL,
                n = length(obj$Y_diss),
                k = ncol(samples))
  
  means <- apply(samples, 2, mean)
  sds <- apply(samples, 2, sd)
  
  # Quantiles
  CI <- t(apply(samples,2, function(x) quantile(x, probs = quantiles)))
  
  # Pseudo Z
  pseudo_z <- (means - null_value) / sds
  
  # Pseudo P
  pseudo_p <- apply(samples, 2, function(x) {
    n <- length(x)
    obs_mean <- mean(x)
    side <- sign(obs_mean - null_value)
    return(sum( (x - null_value)*side < null_value) / n)
  })
  
  
  sig_stars <- ifelse(pseudo_p < 0.001, "***",
                      ifelse(pseudo_p < 0.01, "**",
                             ifelse(pseudo_p < 0.05, "*",
                                    ifelse(pseudo_p < 0.1, ".", " "))))
  
  out <- list(call = obj$call, 
              estimate = means,
              sds = sds,
              pseudo_zval = pseudo_z,
              pseudo_pval = pseudo_p,
              CI = CI,
              n_boot = obj$n_boot,
              quantiles = quantiles,
              AIC = m_AIC$AIC,
              AICc = m_AIC$AICc,
              BIC = m_BIC,
              names_beta =  paste0('diss: ', colnames(obj$form_X$predictors)),
              names_lambda =  paste0('uniq: ', colnames(obj$form_W$predictors)),
              sig = sig_stars,
              logL = logL)
  
  class(out) <- 'summary.bbgdmm'
  return(out)
}

print.summary.gdmm <- function(x, ...) {
  print_title('GDMM SUMMARY', symb = '—')  # call
  cat("call:\n\n")
  call_str <- deparse(x$call)
  cat(" ", call_str, "\n\n", sep = '')
  
  # table
  cat("coeff. summary table:\n\n")
  rownames(x$table)[rownames(x$table) == 'beta'] <- x$names_beta
  rownames(x$table)[rownames(x$table) == 'lambda'] <- x$names_lambda
  rownames(x$table)[rownames(x$table) == 'intercept'] <- '(Intercept)'
  x$table <- data.frame(x$table, x$sig)
  names(x$table) <-  c("Estimate", "Std.Err.", "z value", "Pr(>|z^2|)", " ")
  print(x$table, digits = 3)
  
  cat("---\n")
  cat("signif. codes: '***' <0.001 '**' <0.01 '*' <0.05 '.' <0.1 \n")
  cat("---\n\n")
  # logL, AIC, AICc and BIC
  cat("Marginal log-likelihood: ", x$logL, "\nAIC: ", x$AIC, ", AICc: ", x$AICc, ", BIC: ", x$BIC, "\n")
}


print.summary.bbgdmm <- function(x, ...) {
  print_title('BBGDMM SUMMARY', symb = '—')  # call
  cat("call:\n\n")
  call_str <- deparse(x$call)
  cat(" ", call_str, "\n\n", sep = '')
  
  # table
  cat("coeff. summary table:\n\n")
  rownames(x$CI)[rownames(x$CI) == 'beta'] <- x$names_beta
  rownames(x$CI)[rownames(x$CI)  == 'lambda'] <- x$names_lambda
  rownames(x$CI)[rownames(x$CI)  == 'intercept'] <- '(Intercept)'
  
  colnames(x$CI) <-  paste0(x$quantiles*100, '%')
  
  table <- data.frame(x$estimate,
                      x$sds,
                      x$CI,
                      x$pseudo_zval,
                      ifelse(x$pseudo_pval == 0, paste('<', 1/x$n_boot), x$pseudo_pval),
                      x$sig)
  names(table) <- c('Estimate', 'Std.Err.', paste0(x$quantiles*100,'%'), 'pseudo-Z', 'pseudo-Pval', ' ')
  print(table, digits = 4)
  
  cat("---\n")
  cat("signif. codes: '***' <0.001 '**' <0.01 '*' <0.05 '.' <0.1 \n")
  cat("---\n\n")
  # logL, AIC, AICc and BIC
  cat("Marginal log-likelihood (average): ", x$logL, "\nAIC: ", x$AIC, ", AICc: ", x$AICc, ", BIC: ", x$BIC, "\n")

}





