predict_gdmm <- function(obj,
                         new_X = NULL,
                         new_W = NULL,
                         conditional = FALSE,
                         re = NULL,
                         component = 'dissimilarity',
                         type = 'link'){
  if (is.null(new_X)) {
    new_X <- obj$X
  }
  
  if (is.null(new_W)) {
    new_W <- obj$X
  }
  
  if (is.null(sites)) {
    sites = 1:nrow(obj$Y)
  }

  n <- nrow(new_X)
  
  form_X_new <- as.matrix(forge(new_X,obj$form_X$blueprint)$predictors)
  form_W_new <- as.matrix(forge(new_W,obj$form_W$blueprint)$predictors)
  
  beta <- obj$obj$report()$e_beta
  lambda <- obj$obj$report()$lambda
  intercept <- obj$obj$report()$intercept
  
  D <- t(combn(n, 2))

  uniq_comp = form_W_new %*% lambda
  diss_comp = abs(form_X_new[D[,1], , drop = F] - form_X_new[D[,2],,drop = F]) %*% beta
  
  if (conditional) {
    uniq_comp = uniq_comp + (obj$Z_design %*%  m$obj$report()$u) [sites]
  }
  
  if (component == 'dissimilarity') {
    pred <- intercept + diss_comp + uniq_comp[D[,1],] +  uniq_comp[D[,2],]
    # TO DO: if(type == 'response') pred <- obj$inv_link(pred_eta)
  } else if (component == 'uniqueness') {
    
    diss_matrix <- matrix(ncol = nrow(new_X),nrow = nrow(new_X))
    diss_matrix[cbind(D[,1], D[,2])] <- diss_comp
    diss_matrix[cbind(D[,2], D[,1])] <- diss_comp
    diag(diss_matrix) <- 0
    pred <- (intercept)/2 + rowSums(diss_matrix)/(nrow(diss_matrix) - 1) + uniq_comp  - (sum(diss_matrix)/((nrow(diss_matrix) - 1)^2))/2
    
  }
  
  return(pred)
}

