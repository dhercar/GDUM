library(greta)

# Load functions
dir <- './functions/'
files.sources = list.files(dir)
sapply(paste0(dir,files.sources), source)
set.seed(1)

m0   <- readRDS('models/m0.rds')
m1   <- readRDS('models/m1.rds')
m2_1 <- readRDS('models/m2_1.rds')
m2_2 <- readRDS('models/m2_2.rds')
m3   <- readRDS('models/m3.rds')
m4   <- readRDS('models/m4.rds')
m5   <- readRDS('models/m5.rds')
m6   <- readRDS('models/m6.rds')
m7   <- readRDS('models/m7.rds')

#### CHECK ####
# Check m0
with(m0, {
  print(bayesplot::mcmc_trace(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('e_s'))) #  ok
  coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE) #ok
})

# Check m1
with(m1, {
  print(bayesplot::mcmc_trace(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('e_s'))) #  ok
  coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE) #ok
})

# Check m2.1
with(m2_1, {
  print(bayesplot::mcmc_trace(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('e_s'))) #  ok
  coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE) #ok
})

# Check m2.2
with(m2_2, {
  print(bayesplot::mcmc_trace(draws, regex_pars = c('beta', 'alpha', 'sd'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('beta', 'alpha', 'sd'))) #ok
  coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE) #ok
})

# Check m3
with(m3, {
  print(bayesplot::mcmc_trace(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('e_s'))) #  ok
  coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE) #ok
})

# Check m4
with(m4, {
  print(bayesplot::mcmc_trace(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('e_s'))) #  ok
  coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE) #ok
})

# Check m5
with(m5, {
  print(bayesplot::mcmc_trace(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('e_s'))) #  ok
  coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE) #ok
})

# Check m6
with(m6, {
  print(bayesplot::mcmc_trace(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('e_s'))) #  ok
  coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE) #ok
})

# Check m6
with(m7, {
  print(bayesplot::mcmc_trace(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('beta', 'alpha', 'SD'))) #ok
  print(bayesplot::mcmc_rank_overlay(draws, regex_pars = c('e_s'))) #  ok
  coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE) #ok
})

#### PARAMS ####
with(m0, {
  summary(draws)
})


