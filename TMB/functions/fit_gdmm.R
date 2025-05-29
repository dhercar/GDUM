
gdmm <- function(Y,
                 X, 
                 diss_formula= ~.,
                 uniq_formula= NULL,
                 mono = FALSE,
                 family = 'binomial',
                 link = 'logit',
                 binary = FALSE,
                 method = 'bray',
                 control = NULL,
                 trace = FALSE){
  
  # Checks
  # ToDo: 
  # - No factors / re in diss_formula
  # - valid link
  # - valid family
  
  # Dissimilarity model matrix
  form_X <- hardhat::mold(diss_formula, 
                          X, 
                          blueprint = default_formula_blueprint(intercept = FALSE))
  
  # Uniqueness model matrix
  fix_terms_W <- lme4::nobars(uniq_formula)
  re_vars <- all.vars(uniq_formula)[!(all.vars(uniq_formula) %in% all.vars(fix_terms_W))] 
  form_W <- hardhat::mold(fix_terms_W, 
                          X, 
                          blueprint = default_formula_blueprint(intercept = FALSE, indicators = 'traditional'))
  # Re model matrix
  if (length(re_vars)>0){
    Z_design <- sparse.model.matrix(eval(parse(text = paste0('~ 0 + ', paste0(re_vars, collapse = ' + ')))),
                             data = env)
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

  
  if(family %in% c('binomial')){
    Y_pair <- make_y_df(com = Y, method = method, num_den = TRUE)
    Y_diss <- Y_pair[,3]
    Y_den <- Y_pair[,4]
    D <- Y_pair[,1:2]
    map = list(log_scale = factor(NA))
  } else {
    Y_pair <- make_y_df(com = Y, method = method)
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
    mono = 1,
    link = link_num,
    family = family_num)
  
  parameters <- list(
    intercept = 0,
    beta = rep(0, ncol(form_X$predictors)),
    lambda = rep(0, ncol(form_W$predictors)),
    log_sigma_re = rep(0, length(re_vars)),
    u = rep(0, ncol(Z_design)*has_re),
    log_scale = 0
  )

  obj <- MakeADFun(
    data = data,
    parameters = parameters,
    DLL = 'gdum',
    map = map,
    silent = !trace,
    random = c('u'),
  )
  
  control_def <- list(rel.tol = 1e-10,
                      eval.max = 1000,
                      iter.max = 1000,
                      trace = 0)
  
  if(!is.null(control)){
    control_def[[names(control)]] <- control
  }
  
  opt <- nlminb(
    start = obj$par, 
    objective = obj$fn, 
    gradient = obj$gr,
    control = control_def
  )
  
  print(opt$message)
  
  out <- list(Y = Y,
              Y_den = Y_den,
              X = X, 
              D = D,
              Y_diss = Y_diss,
              form_X = form_X,
              form_W = form_W,
              Z_design = Z_design,
              diss_formula = diss_formula,
              uniq_formula = uniq_formula,
              link = link,
              family = family,
              sdrep = sdreport(obj),
              obj = obj,
              opt = opt)
  
  return(out)
}




