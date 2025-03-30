# Load packages
library(geoR)
library(greta)
library(rvmethod)
library(cowplot)
library(coda)
library(vegan)
library(adespatial)
library(ggplot2)
library(showtext)
library(latex2exp)

# Load functions
dir <- './functions/'
files.sources = list.files(dir)
sapply(paste0(dir,files.sources), source)
set.seed(1)

# plotting preferences
theme_set(theme_bw())
theme_update(panel.grid = element_blank(), text = element_text(size = 7),
             panel.background = element_rect(fill = 'grey90', colour = 'transparent'),
             panel.border = element_rect(colour = 'transparent'))
set.seed(1)

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
#  SIMULATE COMMUNITY
com_data <- sim_com_data(n_sp = n_sp,
                         isolation = env_data$data$isolation,
                         env = env_data$data[,substr(names(env_data$data),1,3) == 'env'],
                         sigma_range = sigma_range, 
                         mu_range = mu_range,
                         ab_mu_range = ab_mu_range,
                         isolation_effect = w_effect,
                         seed = 1)

# PREPARE DATA
n_pred = 100 # number of points calculated to inspect partial effects

# Predictors 
X <- make_x_df(data.frame(env1 = env_data$data$env_1))
X$s_dist_env1 <- scale(X$dist_env1)[,1]

W <- env_data$data
W$s_isolation <- scale(W$isolation)[,1]

# ---- GAUSSIAN GDUM ----
Y <- make_y_df( com_data$com_data > 0, num_den = FALSE, method = 'sorensen')[, 3] # Dissimilarity as b + c / 2*a + b + c
D <- make_y_df(com_data$com_data > 0)[,1:2] # Design matrix (combinations of rows and columns)

m0 <- fit_gdum(Y = Y, 
               X = X, 
               D = D, 
               W = W, 
               diss_formula = ~ s_dist_env1, 
               site_formula = ~s_isolation,
               Lmin = 15,
               Lmax = 20,
               n_samples = 8000,
               warmup = 2000)

# Check chains
bayesplot::mcmc_trace(m0$draws, regex_pars = c('beta', 'lambda', 'SD', 'sigma', 'alpha'))
gelman.diag(m0$draws)

u <- data.frame(predict_gdum(m0, response = 'uniqueness'))
u$X50_no_re <- t(predict_gdum(m0, response = 'uniqueness', re = FALSE, quantiles = c(0.5)))
u$X50_no_x <- t(predict_gdum(m0, W_new = W, response = 'uniqueness', re = FALSE, quantiles = c(0.5)))
 
u$LCBD <- adespatial::LCBD.comp(vegdist(com_data$com_data))$LCBD
u$SS <- u$LCBD*adespatial::LCBD.comp(vegdist(com_data$com_data))$beta['SStotal']

plot_grid(
  ggplot(u, aes(x = LCBD, y = X50_no_x/sum( X50_no_x) )) + 
    geom_point(shape = 21) + 
    theme(aspect.ratio = 1) +
    ylab('E(LCBD)') + 
    ggtitle(TeX('$E(LCBD_i) | \\Delta x_{ij} = \\bar{\\Delta x}, w_i, \\alpha_i = 0$')) +
    ylim(range(c(u$X50./sum( u$X50.) , u$X50_no_re/sum( u$X50_no_re), u$X50_no_x/sum( u$X50_no_x)))) +
    geom_abline(lty = 5),
ggplot(u, aes(x = LCBD, y = X50_no_re/sum( X50_no_re) )) + 
  geom_point(shape = 21) + 
  theme(aspect.ratio = 1) +
  ylab('E(LCBD)') + 
  ylim(range(c(u$X50./sum( u$X50.) , u$X50_no_re/sum( u$X50_no_re), u$X50_no_x/sum( u$X50_no_x)))) +
  ggtitle(TeX('$E(LCBD_i) | \\Delta x_{ij}, w_i, \\alpha_i = 0$')) +
  geom_abline(lty = 5),
ggplot(u, aes(x = LCBD, y = X50./sum( X50.) )) + 
  geom_point(shape = 21) + 
  ylab('E(LCBD)') + 
  theme(aspect.ratio = 1) +
  ylim(range(c(u$X50./sum( u$X50.) , u$X50_no_re/sum( u$X50_no_re), u$X50_no_x/sum( u$X50_no_x)))) +
  ggtitle(TeX('$E(LCBD_i) | \\Delta x_{ij}, w_i, \\alpha_i$')) +
  geom_abline(lty = 5), labels = c('A', 'B', 'C'),
ncol = 3)

ggsave('./SI/one_to_one_plot.png', width = 19, height = 8, dpi = 600, units = 'cm')
