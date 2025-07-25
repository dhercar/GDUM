fit_gdum <- function(Y, # Response (dissimilarity)
                     X = NULL, # Pairwise level predictors
                     W = NULL, # Site level predictors
                     D, # Design (combinations of rows and columns)
                     Y_den = NULL, # Denominator in sorensen or jaccard dissimilarity index
                     diss_formula =~.,
                     site_formula =~.,
                     priors = list(),
                     pos_beta = TRUE,
                     family = 'gaussian',
                     link = 'identity',
                     chains = 4,
                     n_cores = NULL,
                     n_samples = 4000,
                     warmup = 4000,
                     Lmin = 10,
                     Lmax = 15,
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
      X_mat <- X_mat[,-1, drop = FALSE]
    }
  }
  
  if(!is.null(W)){
    W_mat <- model.matrix(site_formula, data = data.frame(W))
    
    if("(Intercept)" %in% colnames(W_mat)){
      W_mat <- W_mat[,-1, drop = FALSE]
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
                            'betabinomial' = beta_binomial(Y_den, shape1, shape2),
                            'binomial' = binomial(Y_den, mu))
  
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
  draws <- greta::mcmc(model_obj, 
                       hmc(Lmin = Lmin, Lmax = Lmax),
                       warmup = warmup,
                       n_samples = n_samples,
                       chains = chains,
                       n_cores = n_cores)
  
  # Change names to match predictor names
  for(i in 1:chains){
    if(!is.null(X_mat)){
    colnames(draws[[i]])[grepl("beta", colnames(draws[[i]]), fixed=TRUE)] <- paste0('beta_', colnames(X_mat))
    }
    
    if(!is.null(W_mat)){
      colnames(draws[[i]])[grepl("lambda", colnames(draws[[i]]), fixed=TRUE)] <- paste0('lambda_', colnames(W_mat))
    }
  }
  
  result <- list(
    model = model_obj,
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
      Y_den  = Y_den
    ),
    link = link,
    ilink = ilink,
    diss_formula = substitute(diss_formula),
    site_formula = substitute(site_formula)
  )
  
  return(result)
}



s <- function(x, degree = 2, df = 3, ... ){
  splines2::isp(x, degree = degree, ...)
}



