# Function to simulate community data 

sim_com_data <- function(n_sp = 50,              # number of species in species pool
                         isolation = NULL,       # vector of isolation values (w) for each site
                         isolation_effect = 0.5, # Effect of isolation on species pool (max % reduction at maximum isolation value)
                         env = NULL,             # Vector of values for environmental gradient (x) 
                         sigma_range = NULL,     # vector with two values indicating the range of species tolerances (sigma of Gaussian function)
                         mu_range = NULL,        # vector with two values indicating the range of species optima (mu of Gaussian function)
                         ab_mu_range = NULL,     # range of species abundance at optima (i.e., species rarity)
                         seed = NULL             # set seed for reproducibility
                         ) {
  
  stopifnot(!is.null(env),
            !is.null(sigma_range),
            !is.null(mu_range),
            !is.null(ab_mu_range))
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  env <- data.frame(env)
  
  # maximum abundance for each species
  max_ab <- runif(n_sp, ab_mu_range[1],ab_mu_range[2])
  names(max_ab) <- paste0('sp_', 1:n_sp)
  
  # simulate species niche
  env_opt <- replicate(ncol(env), runif(n_sp, mu_range[1], mu_range[2])) # Env optima
  colnames(env_opt) <- names(env)
  
  env_tol <- replicate(ncol(env), runif(n_sp, sigma_range[1], sigma_range[2])) # Env tolerance
  colnames(env_tol) <- names(env)
  
  # Average expected abundance on each site
  # Calculate probability of occurrence based on product of functions in each environmental axis
  sp_ab <- lapply(1:n_sp, function(i){
    p_i = 1
    for (j in ncol(env)) {
      p_env_j <- gaussfunc(env[j], 
                           mu = env_opt[i,j], 
                           sigma = env_tol[i,j])
      p_i = p_i*p_env_j
    }
    p_i * max_ab[i] # scale based on species rarity
  })
  sp_ab <-  data.frame(do.call(cbind, sp_ab)) 
  
  # Effect of isolation (random reduction in species pool)
  sp_p0 <- t(vrbinom(n = n_sp, size = 1, prob = isolation_effect*scale01(isolation))) 
  sp_ab <- sp_ab*(1 - sp_p0)
  sp_sample <- data.frame(apply(sp_ab,2, function(x){
    vrbinom(1,1,x)
  }))
  
  # Generate output
  names(sp_sample) <- paste0('sp_', 1:n_sp)
  list(com_data = sp_sample,
       sp_mu = env_opt,
       sp_sigma = env_tol,
       sp_ab_mu = max_ab,
       n_sp_p0 = rowSums(sp_p0))
}



