# Load packages
library(geoR)
library(greta)
library(rvmethod)
library(gridExtra)
library(coda)
library(vegan)
library(adespatial)
library(ggplot2)

# Load functions
dir <- './functions/'
files.sources = list.files(dir)
sapply(paste0(dir,files.sources), source)
set.seed(1)

# Create output directories
dir.create('./models')
dir.create('./SI')

# ---- SIMULATED DATA ----

# Set parameters 
n_sites = 50 # Number on unique sites / samples
n_sp = 100 # Unique species in species pool
field_size = 1000 # For spatial representation
w_effect = 0.5 # probability of removing a species from species pool at max w value
sigma_range = c(0.1,0.5) # range of environmental tolerances
mu_range = c(-0.2, 1.2) # range of environmental optima
ab_mu_range = c(1, 1) # range of peak abundances (mean abundance at env = mu)
n_neighbours = 10 # Neighbours used to compute isolation index

# SIMULATE SPATIAL DATA
env_data <- sim_spatial_data( n_sites = n_sites,
                              field_size = field_size,
                              cov_pars = c(1, 0.1),
                              n_neighbours = n_neighbours,
                              seed = 1)

(map <- ggplot(env_data$data, aes(x=x, y=y)) + 
  geom_tile(data = env_data$field, aes(x = x, y=y, fill =z)) +
  scale_fill_gradientn(colours = c('white', 'grey90')) +
  geom_point(aes(size = -isolation, col = env_1)) +
  geom_point(shape = 21, aes(size = -isolation)) +
  scale_size_continuous(range = c(0.5,3))+
  scale_colour_gradientn(colours = c('wheat','gold','darkorange', 'brown')) +
  theme(aspect.ratio = 1, legend.position = '') +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)))

ggsave('SI/SI_map.png', map, width = 10, height = 10, units = 'cm', dpi = 600)

#  SIMULATE COMMUNITY
com_data <- sim_com_data(n_sp = n_sp,
                         isolation = env_data$data$isolation,
                         env = env_data$data[,substr(names(env_data$data),1,3) == 'env'],
                         sigma_range = sigma_range, 
                         mu_range = mu_range,
                         ab_mu_range = ab_mu_range,
                         isolation_effect = w_effect,
                         seed = 1)


head(com_data$com_data[,1:10])

# PREPARE DATA
n_pred = 100 # number of points calculated to inspect partial effects

# Predictors 
X_e_dis <- scale(make_x_df(data.frame(env_data$data$env_1))[, 3]) # Delta x (pairwise level predictor)
X_iso <- scale(env_data$data$isolation)  # w (site-level predictor)
X_env <- splines::ns(scale(env_data$data$env_1), df = 2) # x (site-level predictor)\
X_e_uni <- scale(adespatial::LCBD.comp(dist(env_data$data$env_1))$LCBD) # u(x) (site-level predictor)

# ---- M0: FULL MODEL - GAUSSIAN ----
# Response
Y <- make_y_df( com_data$com_data > 0, num_den = FALSE, method = 'sorensen')[, 3] # Dissimilarity as b + c / 2*a + b + c
D <- make_y_df(com_data$com_data > 0)[,1:2] # Design matrix (combinations of rows and columns)

# Index variables
row <- D[, 1]
col <- D[, 2]

# PRIORS
# Regression coefficients
beta_e_dis <- normal(0, 1, truncation = c(0, Inf)) 
beta_iso <- normal(0, 1)
alpha <- normal(0, 0.5)
SD_s <- normal(0, 0.5, truncation = c(0, Inf)) 

# Site-level random effect
e_s  <- normal(0, 1, dim = n_sites)*SD_s 
sd <- normal(0, 1, truncation = c(0, Inf))

# Linear predictor
S_s <- X_iso %*% beta_iso + e_s # Uniqueness component
D_s <- X_e_dis %*% beta_e_dis
eta <- alpha + D_s + S_s[row] + S_s[col] # Linear predictor

# (inverse) link
mu <- eta # link 

# # Distribution
distribution(Y) <- normal(mu, sd)

# MCMC
mod0 <- model(beta_e_dis, beta_iso, alpha, e_s, S_s, SD_s)
draws0 <- greta::mcmc(mod0,
                      warmup = 2000,
                      n_samples = 3000,
                      hmc(Lmin = 15, Lmax = 20))

saveRDS(list(mod = mod0, 
             draws = draws0, # model and draws 
             beta_e_dis = beta_e_dis, 
             beta_iso = beta_iso,
             alpha = alpha, 
             e_s = e_s, 
             S_s = S_s, 
             SD_s = SD_s,
             sd = sd,
             row = row, 
             col = col, 
             Y = Y, 
             X_e_dis = X_e_dis,
             X_iso = X_iso,
             n_pred = n_pred),'models/m0.rds')


# ---- M1: NO UNIQUENESS PREDICTORS ----
# Response
Y <- make_y_df( com_data$com_data > 0, num_den = TRUE, method = 'sorensen')[, c(3,4)] # Dissimilarity as two columns (b + c and 2*a + b + c)
D <- make_y_df(com_data$com_data > 0)[,1:2] # Design matrix (combinations of rows and columns)

# dimensions
trials <- Y[,2]
n_realisations <- nrow(Y)

# Priors
beta_e_dis <- normal(0, 1, truncation = c(0, Inf)) # Regression coefficients
alpha <- normal(0, 0.5)
SD_s <- normal(0, 0.5, truncation = c(0, Inf)) # Site re SD
e_s  <- normal(0, 1, dim = n_sites)*SD_s # Site re

# Linear predictor
S_s <- e_s # Uniqueness component
D_s <- X_e_dis %*% beta_e_dis # Dissimilarity component 
eta <- alpha + D_s + S_s[row] + S_s[col]  # Linear predictor

# (inverse) link
mu <- ilogit(eta)

# Distribution
distribution(Y[,1]) <- binomial(trials, mu)

#MCMC
mod1 <- model(beta_e_dis, alpha, e_s, S_s, SD_s)
draws1 <- greta::mcmc(mod1,
              warmup = 2000,
              n_samples = 3000,
              hmc(Lmin = 15, Lmax = 20))


saveRDS(list(mod = mod1, 
             draws = draws1, 
             beta_e_dis = beta_e_dis, 
             alpha = alpha, 
             e_s = e_s, 
             S_s = S_s, 
             SD_s = SD_s,
             row = row, 
             col = col, 
             Y = Y, 
             X_e_dis = X_e_dis,
             n_pred = n_pred),'models/m1.rds')



# ---- M3 FULL MODEL (LCBD) ----
Y <- make_y_df( com_data$com_data > 0, num_den = TRUE)[, 3:4]

# dimensions
trials <- Y[,2]
n_realisations <- nrow(Y)

# Priors
beta_e_dis <- normal(0, 1, truncation = c(0, Inf)) 
beta_iso <- normal(0, 1)
alpha <- normal(0, 0.5)
SD_s <- normal(0, 0.5, truncation = c(0, Inf)) # Site re SD
e_s  <- normal(0, 1, dim = n_sites)*SD_s #Site re

# Linear predictor
S_s <- X_iso %*% beta_iso + e_s # Uniqueness component
D_s <- X_e_dis %*% beta_e_dis # Pairwise component component
eta <- alpha + D_s + S_s[row] + S_s[col] # Linear predictor

# (inverse) link
mu <- ilogit(eta) # link 

# # Distribution
distribution(Y[,1]) <- binomial(trials, mu)

# MCMC
mod3 <- model(beta_e_dis, beta_iso, alpha, e_s, S_s, SD_s)
draws3 <- greta::mcmc(mod3,
                      warmup = 2000,
                      n_samples = 3000,
                      hmc(Lmin = 15, Lmax = 20))

saveRDS(list(mod = mod3, 
             draws = draws3, # model and draws 
             beta_e_dis = beta_e_dis, 
             beta_iso = beta_iso,
             alpha = alpha, 
             e_s = e_s, 
             S_s = S_s, 
             SD_s = SD_s,
             row = row, 
             col = col, 
             Y = Y, 
             X_e_dis = X_e_dis,
             X_iso = X_iso,
             n_pred = n_pred),'models/m3.rds')

# ---- M2.1: GDUM METHOD ----
Y <- make_y_df( com_data$com_data > 0, num_den = TRUE)[, 3:4]

# dimensions

trials <- Y[,2]
n_realisations <- nrow(Y)
# Priors
beta_iso <- normal(0, 1)
alpha <- normal(0, 0.5)
SD_s <- normal(0, 0.5, truncation = c(0, Inf)) # Site re SD
e_s  <- normal(0, 1, dim = n_sites)*SD_s #Site re

# Linear predictor
S_s <- X_iso %*% beta_iso + e_s # Uniqueness component
eta <- alpha +S_s[row] + S_s[col] # Linear predictor

# (inverse) link
mu <- ilogit(eta) # link 

# # Distribution
distribution(Y[,1]) <- binomial(trials, mu)

# MCMC
mod2.1 <- model(beta_iso, alpha, e_s, S_s, SD_s)
draws2.1 <- greta::mcmc(mod2.1,
                      warmup = 2000,
                      n_samples = 3000,
                      hmc(Lmin = 15, Lmax = 20))

saveRDS(list(mod = mod2.1, 
             draws = draws2.1, # model and draws 
             beta_iso = beta_iso,
             alpha = alpha, 
             e_s = e_s, 
             S_s = S_s, 
             SD_s = SD_s,
             row = row, 
             col = col, 
             Y = Y, 
             X_iso = X_iso,
             n_pred = n_pred),'models/m2_1.rds')

# ---- M2.2: CLASSIC METHOD ----

# Ontain SSi as LCBD * SS_tot (LCBD = SS_i / SS_tot)
SS_tot <- adespatial::LCBD.comp(vegan::vegdist(com_data$com_data > 0, method = 'bray'))$beta[1]
Y <- adespatial::LCBD.comp(vegan::vegdist(com_data$com_data > 0, method = 'bray'))$LCBD*SS_tot

beta_iso <- normal(0,1)
alpha <- normal(0,1)
sd <- normal(0, 0.5, truncation = c(0, Inf))

eta <- alpha + X_iso %*% beta_iso
distribution(Y) <- normal(eta, sd)

# MCMC
mod2.2 <- model(beta_iso, alpha, sd)
draws2.2 <- greta::mcmc(mod2.2,
                      warmup = 2000,
                      n_samples = 2000,
                      hmc(Lmin = 10, Lmax = 15))

saveRDS(list(mod = mod2.2, 
             draws = draws2.2,
             beta_iso = beta_iso,
             alpha = alpha, 
             X_iso = X_iso,
             sd = sd,
             SS_tot = SS_tot,
             Y = Y,
             n_pred = n_pred),'models/m2_2.rds')

# ---- ENVIRONMETNAL UNIQUENESS ----

# ---- M4 ----
Y <- make_y_df( com_data$com_data > 0, num_den = TRUE)[, 3:4]

# priors
beta_e_uni <- normal(0, 1)
beta_iso <- normal(0, 1)
alpha <- normal(0, 0.5)
SD_s <- normal(0, 0.5, truncation = c(0, Inf))
e_s  <- normal(0, 1, dim = n_sites) * SD_s

# Linear predictor
S_s <- X_iso %*% beta_iso + X_e_uni %*% beta_e_uni + e_s  # Uniqueness component
D_s <- 0
eta <- alpha + D_s + S_s[row] + S_s[col] # Linear predictor

# link
mu <- ilogit(eta)  


trials <- Y[,2]
n_realisations <- nrow(Y)

distribution(Y[,1]) <- binomial(trials, mu)

# MCMC
mod4 <- model(beta_e_uni, beta_iso, alpha, e_s, S_s, SD_s)

draws4 <- greta::mcmc(mod4,
                      warmup = 2000,
                      n_samples = 3000,
                      hmc(Lmin = 15, Lmax = 20))


saveRDS(list(mod = mod4, 
             draws = draws4, # model and draws 
             beta_e_uni = beta_e_uni, 
             beta_iso = beta_iso,
             alpha = alpha, 
             e_s = e_s, 
             S_s = S_s, 
             SD_s = SD_s,
             row = row, 
             col = col, 
             Y = Y, 
             X_iso = X_iso,
             X_e_uni = X_e_uni
),'models/m4.rds')


# ---- M5 ----

# priors
beta_e_dis <- normal(0, 1)
beta_e_uni <- normal(0, 1)
beta_iso <- normal(0, 1)
alpha <- normal(0, 0.5)
SD_s <- normal(0, 0.5, truncation = c(0, Inf))
e_s  <- normal(0, 1, dim = n_sites) * SD_s

# Linear predictor
S_s <- X_iso %*% beta_iso + X_e_uni %*% beta_e_uni + e_s # Uniqueness component
D_s <- X_e_dis %*% beta_e_dis
eta <- alpha + D_s + S_s[row] + S_s[col] # Linear predictor
mu <- ilogit(eta) # link 

trials <- Y[,2]
n_realisations <- nrow(Y)

distribution(Y[,1]) <- binomial(trials, mu)

# MCMC
mod5 <- model(beta_e_uni, beta_iso, beta_e_dis, alpha, e_s, S_s, SD_s)
draws5 <- greta::mcmc(mod5,
                      warmup = 2000,
                      n_samples = 3000,
                      hmc(Lmin = 15, Lmax = 20))

saveRDS(list(mod = mod5, 
             draws = draws5, # model and draws 
             beta_e_uni = beta_e_uni, 
             beta_e_dis = beta_e_dis,
             beta_iso = beta_iso,
             alpha = alpha, 
             e_s = e_s, 
             S_s = S_s, 
             SD_s = SD_s,
             row = row, 
             col = col, 
             Y = Y, 
             X_iso = X_iso,
             X_e_uni = X_e_uni,
             X_e_dis = X_e_dis),
        'models/m5.rds')

# ---- RAW ENVIRONMETNAL VALUES ----

# ---- M6 ----

# priors
beta_env  <- normal(0,1, dim = ncol(X_env))
beta_iso <- normal(0, 1)
alpha <- normal(0, 0.5)
SD_s <- normal(0, 0.5, truncation = c(0, Inf))
e_s  <- normal(0, 1, dim = n_sites) * SD_s

# Linear predictor
S_s <- X_iso %*% beta_iso + X_env %*% beta_env + e_s  # Uniqueness component
D_s <- 0
eta <- alpha + D_s + S_s[row] + S_s[col] # Linear predictor
mu <- ilogit(eta) # link 

trials <- Y[,2]
n_realisations <- nrow(Y)

distribution(Y[,1]) <- binomial(trials, mu)

# MCMC
mod6 <- model(beta_env, beta_iso, alpha, e_s, S_s, SD_s)
draws6 <- greta::mcmc(mod6,
                      warmup = 2000,
                      n_samples = 3000,
                      hmc(Lmin = 15, Lmax = 20))

saveRDS(list(mod = mod6, 
             draws = draws6, # model and draws 
             beta_env = beta_env, 
             beta_iso = beta_iso,
             alpha = alpha, 
             e_s = e_s, 
             S_s = S_s, 
             SD_s = SD_s,
             row = row, 
             col = col, 
             Y = Y, 
             env_raw = env_data$data$env_1,
             X_iso = X_iso,
             X_env = X_env),
        'models/m6.rds')

# ---- M7 ----
# priors
beta_e_dis <- normal(0, 1)
beta_env  <- normal(0,1, dim = ncol(X_env))
beta_iso <- normal(0, 1)
alpha <- normal(0, 0.5)
SD_s <- normal(0, 0.5, truncation = c(0, Inf))
e_s  <- normal(0, 1, dim = n_sites) * SD_s

# Linear predictor
S_s <- X_iso %*% beta_iso + X_env %*% beta_env + e_s  # Uniqueness component
D_s <- X_e_dis %*% beta_e_dis
eta <- alpha + D_s + S_s[row] + S_s[col] # Linear predictor
mu <- ilogit(eta) # link 

trials <- Y[,2]
n_realisations <- nrow(Y)

# Distribution
distribution(Y[,1]) <- binomial(trials, mu)

# MCMC
mod7 <- model(beta_env, beta_iso, beta_e_dis, alpha, e_s, S_s, SD_s)
draws7 <- greta::mcmc(mod7,
                      warmup = 2000,
                      n_samples = 3000,
                      hmc(Lmin = 15, Lmax = 20))

saveRDS(list(mod = mod7, 
             draws = draws7, # model and draws 
             beta_env = beta_env, 
             beta_e_dis = beta_e_dis,
             beta_iso = beta_iso,
             alpha = alpha, 
             e_s = e_s, 
             S_s = S_s, 
             SD_s = SD_s,
             row = row, 
             col = col, 
             Y = Y, 
             X_iso = X_iso,
             X_e_dis = X_e_dis,
             env_raw = env_data$data$env_1,
             X_env = X_env),
        'models/m7.rds')

