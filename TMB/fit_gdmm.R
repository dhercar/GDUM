
gdmm <- function(Y,
                 X, 
                 diss_formula= ~.,
                 uniq_formula= NULL,
                 mono = FALSE,
                 family = 'normal',
                 link = 'identity',
                 binary = FALSE,
                 method = 'bray',
                 control = NULL,
                 trace = FALSE,
                 bboot = FALSE,
                 n_boot = 1000,
                 n_cores = NULL){
  require(doSNOW)
  # Checks
  # ToDo: 
  # - No factors / re in diss_formula
  # - valid link
  # - valid family
  
  # Dissimilarity model matrix
  if (!is.null(diss_formula) & !is.null(X)) {
  form_X <- hardhat::mold(diss_formula, 
                          X, 
                          blueprint = default_formula_blueprint(intercept = FALSE))
  } else {
  form_X <- list(predictors =  matrix(ncol = 0, nrow = 0))
  }
  
  # Uniqueness model matrix
  if (!is.null(uniq_formula) & !is.null(X)) {
    fix_terms_W <- lme4::nobars(uniq_formula)
    re_vars <- all.vars(uniq_formula)[!(all.vars(uniq_formula) %in% all.vars(fix_terms_W))] 
    
    if ( length(attr(terms(fix_terms_W), 'variables')) > 1 ) {
      form_W <- hardhat::mold(fix_terms_W, 
                              X, 
                              blueprint = default_formula_blueprint(intercept = FALSE, indicators = 'traditional'))
    } else {
      form_W <- list(predictors =  matrix(ncol = 0, nrow = 0))
    }
    
  } else {
    re_vars <- character(0)
    form_W <- list(predictors =  matrix(ncol = 0, nrow = 0))
  }
  
  # Re model matrix
  if (length(re_vars) > 0) {
    Z_design <- sparse.model.matrix(eval(parse(text = paste0('~ 0 + ', paste0(re_vars, collapse = ' + ')))),
                             data = X)
    has_re = 1
    map_re <- as.numeric(as.factor(unlist(lapply(re_vars, function(x) rep(x, length(unique(X[[x]])))))))
    
  }else{
    Z_design <- Matrix::sparseMatrix(
      i = integer(0),  
      j = integer(0), 
      x = numeric(0), 
      dims = c(ncol(Y), 1))
    has_re = 0
    map_re = 0
  }

  if (family %in% c('binomial')) {
    Y_pair <- make_y_df(com = Y, method = method, binary = binary, num_den = TRUE)
    Y_diss <- Y_pair[,3]
    Y_den <- Y_pair[,4]
    D <- Y_pair[,1:2]
    map = list(log_scale = factor(NA))
    
  } else {
    Y_pair <- make_y_df(com = Y, method = method, binary = binary, num_den = FALSE)
    Y_diss <- Y_pair[,3]
    Y_den <- 0
    D <- Y_pair[,1:2]
    map = list()
  }
  
  family_num <- switch(family,
                       'normal' = 0,
                       'binomial' = 1,
                       'beta' = 2)
  link_num <- switch(link,
                     'identity' = 0,
                     'logit' = 1,
                     'neg_exp' = 2,
                     'neg_gaus' = 3)
  
  data <- list(
    Y = Y_diss,
    Y_den = Y_den,
    D = as.matrix(D - 1), # remove 1 to match c++
    X = as.matrix(form_X$predictors),
    W = as.matrix(form_W$predictors),
    Z = Z_design,
    has_random = has_re,
    map_re = map_re - 1,  # remove 1 to match c++ indexing
    mono = 0,
    link = link_num,
    family = family_num)
  
  parameters <- list(
    intercept = 0,
    beta = rep(0, ncol(form_X$predictors)),
    lambda = rep(0, ncol(form_W$predictors)),
    log_sigma_re = rep(0, length(re_vars)),
    u = rep(0, ncol(Z_design)*has_re),
    log_scale = 0)
  
  control_def <- list(rel.tol = 1e-10,
                      eval.max = 1000,
                      iter.max = 1000,
                      trace = 0)
  
  if (!is.null(control)) {
    control_def[[names(control)]] <- control
  }
  
  
  if (bboot == FALSE) {
    data[['weights']] <- rep(1, length(Y_diss))
    
    obj <- MakeADFun(
      data = data,
      parameters = parameters,
      DLL = 'gdum',
      map = map,
      silent = !trace,
      random = c('u'),
    )
    
    n_par <- length(obj$par)
    lower_bounds <- rep(-Inf, n_par)
    upper_bounds <- rep( Inf, n_par)
    
    if (mono) {
      lower_bounds[names(obj$par) == 'beta'] <- 0  
    }

    opt <- nlminb(
      start = obj$par, 
      objective = obj$fn, 
      gradient = obj$gr,
      lower = lower_bounds,
      upper = upper_bounds,
      control = control_def)
      
      print(opt$message)
      
      out <- list(Y = Y,
                  Y_diss = Y_diss,
                  Y_den = Y_den,
                  X = X, 
                  D = D,
                  form_X = form_X,
                  form_W = form_W,
                  Z_design = Z_design,
                  map_re = map_re,
                  re_vars = re_vars,
                  diss_formula = diss_formula,
                  uniq_formula = uniq_formula,
                  link = link,
                  family = family,
                  sdrep = sdreport(obj),
                  obj = obj,
                  opt = opt,
                  boot = FALSE,
                  mono = mono,
                  call = match.call())
      class(out) <- 'gdmm'
      return(out)
    
  } else if (bboot == TRUE) {
    
    if (is.null(n_cores))  n_cores = parallel::detectCores(logical = TRUE) - 2 
    
    cat('Running bayesian bootstrapping on', n_cores, 'cores ...')
    
    # Make cluster
    cl <- parallel::makeCluster(n_cores, type = "SOCK")
    doSNOW::registerDoSNOW(cl)
    on.exit(parallel::stopCluster(cl))  # Ensure cleanup
    
    # Add necessary things to cluster
    parallel::clusterEvalQ(cl, {
    dyn.load(TMB::dynlib('TMB/gdum'))
    })
    
    # Initialise progress bar
    pb <- txtProgressBar(max = n_boot, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    results <- foreach(i = 1:n_boot, .options.snow = opts, .packages = c('TMB', 'gtools'), .combine = rbind) %dopar% {
      # calculate random weights
      weights_sample <-  c(gtools::rdirichlet(1,rep(1,nrow(Y))))
      weights_pair <- (weights_sample[D[,1]] + weights_sample[D[,2]]) 
      data[['weights']] <- weights_pair/sum(weights_pair)*length(Y_diss) # standardise weights 
      # objective
      obj <- MakeADFun(
        data = data,
        parameters = parameters,
        DLL = 'gdum',
        map = map,
        silent = !trace,
        random = c('u'),
      )
      
      n_par <- length(obj$par)
      lower_bounds <- rep(-Inf, n_par)
      upper_bounds <- rep( Inf, n_par)
      
      if (mono) {
        lower_bounds[names(obj$par) == 'beta'] <- 0  
      }
      
      opt <- nlminb(
        start = obj$par, 
        objective = obj$fn, 
        gradient = obj$gr,
        lower = lower_bounds,
        upper = upper_bounds,
        control = control_def)
      
      c(opt$par, u_re_ = obj$report()$u_0, logLikelihood = -opt$objective)
    }
    
    close(pb)
    
    out <- list(Y = Y,
                Y_diss = Y_diss,
                Y_den = Y_den,
                X = X, 
                D = D,
                form_X = form_X,
                form_W = form_W,
                Z_design = Z_design,
                map_re = map_re,
                re_vars = re_vars,
                diss_formula = diss_formula,
                uniq_formula = uniq_formula,
                link = link,
                family = family,
                boot_samples = results,
                n_boot = n_boot,
                boot = TRUE,
                mono = mono,
                call = match.call()) 
    class(out) <- 'bbgdmm'
    return(out)
  }
}



print.gdmm <- function(m, ...){
  print_title('Generalized dissimilarity mixed model (GDMM)', symb = '—')  # call
  
  cat('convergence:', ifelse(m$opt$convergence == 0, 'succesful convergence\n', 'optimization has not reached succesful convergence\n'))
  cat('nlminb message: ')
  cat(m$opt$message, '\n\n')
  
  cat('function call:\n')
  print(m$call)
  cat('\n')
  
  print_title2(' Details ', symb = '-')
  if (m$mono) cat('- monotonic dissimilarity effects (mono = TRUE) \n')
  cat('- distribution (family):', m$family, '\n')
  cat('- link:', m$link, '\n\n')
  
  print_title2(' Predictors ', symb = '-')
  
  cat('- dissimilarity gradient(s): ')
  if (!is.null(m$diss_formula)) {
    print(m$diss_formula)
  } else {
    cat('no dissimilarity gradients included \n')
  }
  
  cat('- effect(s) on uniqueness: ')
  if (!is.null(m$uniq_formula)) {
    print(lme4::nobars(m$uniq_formula))
  } else {
    cat('no predictors on uniqueness included \n')
  }
  
  cat('- random effects: ')
  cat(ifelse(length(m$re_vars) > 0, paste(m$re_vars,sep = ' ,'), 'no random effects included \n'))
  
}

print.bbgdmm <- function(m, ...){
  print_title('Generalized dissimilarity mixed model (GDMM)', symb = '—')  # call
  
  cat('Bayesian Bootstrapping with', m$n_boot, 'random samples\n\n')
  
  cat('function call:\n')
  print(m$call)
  cat('\n')
  
  print_title2(' Details ', symb = '-')
  if (m$mono) cat('- monotonic dissimilarity effects (mono = TRUE) \n')
  cat('- distribution (family):', m$family, '\n')
  cat('- link:', m$link, '\n\n')
  
  print_title2(' Predictors ', symb = '-')
  
  cat('- dissimilarity gradient(s): ')
  if (!is.null(m$diss_formula)) {
    print(m$diss_formula)
  } else {
    cat('none included \n')
  }
  
  cat('- effect(s) on uniqueness: ')
  if (!is.null(m$uniq_formula)) {
    print(lme4::nobars(m$uniq_formula))
  } else {
    cat('none included \n')
  }
  
  cat('- random effects: ')
  cat(ifelse(length(m$re_vars) > 0, paste(m$re_vars,sep = ' ,'), 'none included \n'))
  
  cat('\n')
  print_title2('', symb = '—')
}



