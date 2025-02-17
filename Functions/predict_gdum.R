predict_gdum <- function(fit, 
                         X_new = NULL, 
                         W_new = NULL, 
                         D_new = NULL, 
                         re = TRUE, 
                         samples = 1000,
                         response = 'dissimilarity',
                         quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  # Use original data only if BOTH X and W are missing
  if (is.null(X_new) && is.null(W_new)) {
    X_new <- fit$data$X
    W_new <- fit$data$W
    D_new <- fit$data$D
    if(re == FALSE){
      fit$greta_arrays$e_s <- zeros(length(unique(c(D_new[,1], D_new[,2]))))
    }
  } 
  
  if (is.null(D_new)) {
    if (!is.null(X_new)) {
      D_new <- matrix(1:(dim(X_new)[1] * 2), ncol = 2)
    } else if (!is.null(W_new)) {
      D_new <- t(combn(1:dim(W_new)[1], 2))
    }
  }
  
  if (is.na(re))  {
    fit$greta_arrays$e_s <-
      normal(0, fit$greta_arrays$SD_s, length(unique(c(D_new[, 1], D_new[, 2]))))
  } else if (re == FALSE) {
    fit$greta_arrays$e_s <-
      zeros(length(unique(c(D_new[, 1], D_new[, 2]))))
    
  }
  

  # Ensure column order matches the original model
  if (!is.null(X_new)) {
    X_mat <- model.matrix(eval(fit$diss_formula), data = X_new)
    if ("(Intercept)" %in% colnames(X_mat)) {
      X_mat <- X_mat[, -1, drop = FALSE]
    }
  } else {
    X_mat = NULL
  }
  
  if (!is.null(W_new)) {
    W_mat <- model.matrix(eval(fit$site_formula), data = W_new)
    if("(Intercept)" %in% colnames(W_mat)){
      W_mat <- W_mat[,-1, drop = FALSE]
    }
  } else{
    W_mat = NULL
  }
  
  # Compute pairwise effect
  new_h_eta <- if (!is.null(X_mat) && !is.null(fit$greta_arrays$beta)) X_mat %*% fit$greta_arrays$beta else 0
  
  # Compute site-level effect
  new_v_eta <- if (!is.null(W_mat) && !is.null(fit$greta_arrays$lambda)) {
    new_v <- W_mat %*% fit$greta_arrays$lambda + fit$greta_arrays$e_s
  } else fit$greta_arrays$e_s
  
  # Compute predictions
  if(response == 'dissimilarity'){
  new_mu <- fit$ilink(new_v_eta[D_new[,1]] + new_v_eta[D_new[,2]] + new_h_eta + fit$greta_arrays$alpha)
  draws_mu <- t(apply(data.frame(calculate(new_mu, values = fit$draws, nsim = samples)), 2, function(x) quantile(x, probs = quantiles)))
  return(draws_mu)
  }
  
  if(response == 'uniqueness'){
  D_s_m <- zeros(length(unique(c(D_new[,1], D_new[,2]))), length(unique(c(D_new[,1], D_new[,2]))))
  D_s_m[cbind(D_new[,1], D_new[,2])]  <-  new_h_eta
  D_s_m[cbind(D_new[,2], D_new[,1])]  <-  new_h_eta
  D_s_s <- greta::rowSums(D_s_m) / nrow(D_s_m)
  u_eta <- fit$greta_arrays$alpha/2 + D_s_s + new_v_eta
  draws_u <- t(apply(data.frame(calculate(u_eta, values = fit$draws, nsim = samples)), 2, function(x) quantile(x, probs = quantiles)))
  return(draws_u)
  }
  
  return(list(mu = draws_mu, u = draws_u))
}

