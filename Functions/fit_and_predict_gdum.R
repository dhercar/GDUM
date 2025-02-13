fit_gdum <- function(Y, # Response (dissimilarity)
                     X = NULL, # Pairwise level predictors
                     W = NULL, # Site level predictors
                     D, # Design (combinations of rows and columns)
                     trials = NULL, # Denominator in sorensen or jaccard dissimilarity index
                     diss_formula =~.,
                     site_formula =~.,
                     priors = list(),
                     pos_beta = TRUE,
                     family = 'betabinomial',
                     link = 'logit',
                     ...) {
                   
  # basic checks
  stopifnot(is.null(X) | is.data.frame(data.frame(X)), 
            is.null(W) | is.data.frame(W),
            is.data.frame(data.frame(D)),
            family %in% c('gaussian', 'binomial', 'beta', 'betabinomial'),
            link %in% c('identity', 'logit', 'probit'))
  
  # inverse link
  ilink <- switch(link,
                  'identity' = function(x) x,
                  'logit' = function(x) ilogit(x),
                  'probit' = function(x) iprobit(x))
  
  # X and W model matrices
  if(!is.null(X)){
    X_mat <- model.matrix(diss_formula, data = data.frame(X))
    if("(Intercept)" %in% colnames(X_mat)){
      X_mat <- as.matrix(X_mat[,-1])
    }
  }

  if(!is.null(W)){
    W_mat <- model.matrix(site_formula, data = data.frame(W))
    if("(Intercept)" %in% colnames(W_mat)){
      W_mat <- as.matrix(W_mat[,-1])
    }
  }

  # Define default priors
  # Define priors as separate objects (Greta arrays cannot be inside lists)
  alpha <- normal(0, 1)
  SD_s <- gamma(1, 1)
  
  phi <- if (family %in% c('beta', 'betabinomial')) exponential(0.01) else NULL
  sigma <- if (family == 'gaussian') normal(0, 1, truncation = c(0, Inf)) else NULL
  
  beta <- if (!is.null(X)) {
    if (pos_beta) normal(0, 1, dim = ncol(X_mat), truncation = c(0, Inf))
    else normal(0, 1, dim = ncol(X_mat))
  } else NULL
  
  lambda <- if (!is.null(W)) normal(0, 1, dim = ncol(W_mat)) else NULL
  
  # User-defined priors (overriding defaults where necessary)
  if (!is.null(priors$alpha)) alpha <- priors$alpha
  if (!is.null(priors$SD_s)) SD_s <- priors$SD_s
  if (!is.null(priors$phi)) phi <- priors$phi
  if (!is.null(priors$sigma)) sigma <- priors$sigma
  if (!is.null(priors$beta)) beta <- priors$beta
  if (!is.null(priors$lambda)) lambda <- priors$lambda
  
  # Site level re
  s1 <- D[,1]
  s2 <- D[,2]
  e_s <- normal(0,1, dim = length(unique(c(s1, s2))))*SD_s
  
  # site_level effect
  if(!is.null(W)){
   v <- W_mat%*%lambda + e_s
   v_eta <- v[s1] + v[s2]
  } else { 
   v <- e_s
   v_eta <- v[s1] + v[s2]
  }
  
  # pairwise effect
  if(!is.null(X)){
    h_eta <- X_mat%*%beta 
  } else { 
    h_eta <- zeros(1)
  }
  
  # Linear predictor
  eta <- v_eta + h_eta + alpha
  
  # mu
  mu <- ilink(eta)
  
  if(family %in% c('beta', 'betabinomial')){
    shape1 <- mu*phi
    shape2 <- (1 - mu)*phi
  }
  
  distribution(Y) <- switch(family,
                            'gaussian' = normal(mu,sigma),
                            'beta' = greta::beta(shape1, shape2),
                            'betabinomial' = beta_binomial(trials, shape1, shape2),
                            'binomial' = binomial(trials, mu))
  
  # Model call: 
  arg_strings <- c() 
  if (!is.null(alpha))  arg_strings <- c(arg_strings, "alpha = alpha")
  if (!is.null(SD_s))   arg_strings <- c(arg_strings, "SD_s = SD_s")
  if (!is.null(phi))    arg_strings <- c(arg_strings, "phi = phi")
  if (!is.null(sigma))  arg_strings <- c(arg_strings, "sigma = sigma")
  if (!is.null(beta))   arg_strings <- c(arg_strings, "beta = beta")
  if (!is.null(lambda)) arg_strings <- c(arg_strings, "lambda = lambda")
  if (!is.null(e_s))    arg_strings <- c(arg_strings, "e_s = e_s")
  
  # Combine the arguments into a call 
  call_str <- paste0("model(", paste(arg_strings, collapse = ", "), ")")
  
  # Evaluate the string \
  model_obj <- eval(parse(text = call_str))
  draws <- greta::mcmc(model_obj, ...)
  
  result <- list(
    draws = draws,
    greta_arrays = list(
      alpha = alpha,
      SD_s  = SD_s,
      phi   = phi,
      sigma = sigma,
      beta  = beta,
      lambda = lambda,
      e_s   = e_s
    ),
    data = list(
      X      = X,
      W      = W,
      D      = D,
      trials = trials
    ),
    link = link,
    ilink = ilink,
    diss_formula = substitute(diss_formula),
    site_formula = substitute(site_formula)
  )
  
  return(result)
}

predict_gdum <- function(fit, 
                         X_new = NULL, 
                         W_new = NULL, 
                         D_new = NULL, 
                         re = FALSE, 
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
  
  if(is.null(D_new)){
    if (!is.null(X_new)) {
      D_new <- matrix(1:(dim(X_new)[1]*2), ncol = 2)
      if (re == FALSE) {
        fit$greta_arrays$e_s <- zeros(prod(dim(D_new)))
      } else {
        fit$greta_arrays$e_s <- normal(0, fit$greta_arrays$SD_s, prod(dim(D_new)))
      }
    } else if (!is.null(W_new)) {
      D_new <- t(combn(1:dim(W_new)[1],2))
      if (re == FALSE) {
      fit$greta_arrays$e_s <- normal(0, fit$greta_arrays$SD_s, dim(W_new)[1])
      } else {
      fit$greta_arrays$e_s <- zeros(dim(W_new)[1])
    }
    }
    }
  

  # Ensure column order matches the original model
  if (!is.null(X_new)) {
    X_mat <- model.matrix(eval(fit$diss_formula), data = X_new)
    X_mat <- matrix(if ("(Intercept)" %in% colnames(X_mat)) X_mat[, -1] else X_mat)
  }else{
    X_mat = NULL
  }
  
  if (!is.null(W_new)) {
    W_mat <- model.matrix(fit$site_formula, data = W_new)
    W_mat <- matrix(if ("(Intercept)" %in% colnames(W_mat)) W_mat[, -1] else W_mat)
  }else{
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
  draws_mu <- cbind(draws_mu, X_new)
  return(draws_mu)
  }
  
  if(response == 'uniqueness'){
  D_s_m <- zeros(length(unique(c(D_new[,1], D_new[,2]))), length(unique(c(D_new[,1], D_new[,2]))))
  D_s_m[cbind(D_new[,1], D_new[,2])]  <-  new_h_eta
  D_s_m[cbind(D_new[,2], D_new[,1])]  <-  new_h_eta
  D_s_s <- greta::rowSums(D_s_m) / nrow(D_s_m)
  u_eta <- fit$greta_arrays$alpha/2 + D_s_s + new_v_eta
  draws_u <- t(apply(data.frame(calculate(u_eta, values = fit$draws, nsim = samples)), 2, function(x) quantile(x, probs = quantiles)))
  draws_u <- cbind(draws_u, W_new)
  return(draws_u)
  }
  
  return(list(mu = draws_mu, u = draws_u))
}

