
#' Title
#'
#' @param m 
#' @param var 
#' @param n 
#' @param CI 
#' @param n_sim 
#' @param CI_quant 
#'
#' @return
#' @export
#'
#' @examples
diss_gradient <- function(m, 
                          var = 'all',
                          n = 100,
                          CI = TRUE,
                          n_sim = NULL,
                          CI_quant = c(0.95)){
  
  all_vars <- all.vars(m$diss_formula)
  
  if (var == 'all') {
    var <- all.vars(m$diss_formula)
  }
  
  stopifnot(all(var %in% all.vars(m$diss_formula)))
  
  # Mean beta
  if (class(m) == 'gdmm') {
    mean_beta <- m$sdrep$value[names(m$sdrep$value) == 'e_beta']
  } else if (class(m) == 'bbgdmm') {
    mean_beta <- apply(m$boot_samples[,colnames(m$boot_samples) %in% c('e_beta'), drop = F], 2, mean)
  }
  
  # SIM
  if (CI) {
    # Param combinations 
    if (class(m) == 'bbgdmm') {
      if (is.null(n_sim)) n_sim = m$n_boot
      sims <- m$boot_samples[,colnames(m$boot_samples) %in% c('e_beta'), drop = F]
      
    } else if (class(m) == 'gdmm')  { 
      if (is.null(n_sim)) n_sim = 1000
      sims <- MASS::mvrnorm(n_sim, m$sdrep$value, m$sdrep$cov)
    }
    sims <- sims[,colnames(sims) == 'e_beta', drop = F]
  }
  
  # prep x data 
  x_mean <- apply(m$X[,all_vars, drop = F], 2, min)
  x_min <- matrix(rep(x_mean, each = n), nrow = n, ncol =  length(x_mean), byrow = F)
  colnames(x_min) <- all_vars
  
  # mean f(x)
  out_list <- list()
  
  for (v in var) {
    x_new <- x_min
    x_new[,v] <- seq(min(m$X[,v]), max(m$X[,v]), length.out = n)
    x_new <- apply(x_new, 2, as.numeric)
    suppressWarnings({
      f_x <- as.matrix(forge(x_new, m$form_X$blueprint)$predictors) %*% mean_beta
    })
    
    out_i <- data.frame(var = v,
                        f_x = f_x,
                        x = x_new[,v])
    
    if (CI) {
      quantiles <- sort(c((1 - CI_quant)*0.5, CI_quant + (1 - CI_quant)*0.5))
      samples <- do.call(cbind, lapply(1:n_sim, function(i){
        beta_i <- sims[i,,drop = F]
        suppressWarnings({
          as.matrix(forge(x_new, m$form_X$blueprint)$predictors) %*% c(beta_i)
        })
      }))
      
      CI_i <- t(apply(samples, 1, function(x){quantile(x, probs = quantiles)}))
      colnames(CI_i) <- c(paste0('CI ', colnames(CI_i)))
      
      out_i <- cbind(out_i, CI_i)
    }
    out_list[[v]] <- out_i
  }
  
  return(out_list)
}






