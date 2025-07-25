# Function to simulate environmental data (x and w)
sim_spatial_data <- function(n_sites = 50,         # number of sites
                             field_size = 1000,    # Size of the grid obtained from random Gaussian field
                             cov_pars = c(1, 0.1), # Parameters for generating random Gaussian field (passed to cov.pars in geoR::grf)
                             n_env = 1,            # number of environmental gradients
                             n_neighbours = NULL,  # number of nearest neighbors used to compute site isolation index (w)
                             seed = NULL           # Set seed for reproducibility
                             ) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (is.null(n_neighbours)) {
    n_neighbours <- n_sites - 1
  }
  
  # Generate Gaussian random field
  field <- grf(n = field_size, 
               grid = 'reg', 
               cov.model = "matern", 
               cov.pars = cov_pars)
  
  # bound grid values between 0-1
  df <- data.frame(field$coords, 
                   z = inv_logit(field$data*5)) 
  
  # draw n_sites samples with probability proportional to z (i.e.,intensity of the field)
  sites <- sample(1:nrow(df),size = n_sites, prob = df$z)
  
  # Generate n_env uncorrelated environmental gradients 
  env <- data.frame(replicate(n_env,runif(n_sites)))
  names(env) <- paste0('env_', 1:n_env)
  
  # Compute isolation as average distance to nearest neighbors
  isolation <- apply(df[sites,1:2], 1, function(x){
    dist <- sqrt((x[1] - df[sites,1])^2 + (x[2] - df[sites,2])^2 )
    dist <- sort(dist[dist > 0])[1:n_neighbours]
    mean(dist)
  })
  
  # Output
  list(data = data.frame(df[sites,1:2],isolation, env),
       field = df)
}
