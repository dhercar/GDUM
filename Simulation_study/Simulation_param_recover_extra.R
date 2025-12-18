# Author: Daniel Hernandez Carrasco
# Email: dani.hc97@gmail.com
# Created: 17/12/2025
# License: MIT (see LICENSE file for details)

# Description: This script performs additional simulations comparing GDUM to conventional models of community uniquenes
# including: (1) Correlated vs. uncorrelated variables
#            (2) With vs. without gradient  

# packages
library(gdmmTMB)
library(doParallel)
library(ggplot2)
library(tidyverse)
library(adespatial)

n_cores = 2  # specify cores

#plotting preferences
theme_set(theme_bw())
theme_update(strip.background = element_blank(),
             text = element_text(size = 8),
             panel.grid = element_blank())

# simulating data and fitting models given specified parameters
sim_gdmm_data <- function(n = NULL, # sample size
                          lambda = NULL, # effect size (uniqueness)
                          beta = NULL, # effect size (dissimilarity),
                          include_beta = FALSE,
                          m = NULL,
                          beta_0 = NULL,
                          sd_pair = NULL,
                          sd_site= NULL, 
                          skew_X = FALSE,
                          rho = 0,
                          n_cores = 10,
                          mod = c(1,2,3)) {
  
  beta = rep(beta, m)
  lambda = rep(lambda, m)
  
  D <- t(combn(n,2))
  n_pair = nrow(D)
  
  re_site = scale(rnorm(n))*sd_site
  err_pair = scale(rnorm(n_pair))*sd_pair
  
  Sigma <- matrix(rho, m, m)
  diag(Sigma) <- 1
  
  X = MASS::mvrnorm(n = n, mu = rep(0,m), Sigma = Sigma, empirical = TRUE)
  if(skew_X) X = apply(X, 2, function(x) scale(log1p(x - min(x))))
  
  diss_comp <- abs(X[D[,1],,drop = F] - X[D[,2],,drop = F]) %*% beta
  uniq_comp <- X %*% lambda + re_site
  
  obs <- beta_0 + diss_comp + uniq_comp[D[,1]] + uniq_comp[D[,2]] + err_pair
  
  obs_matrix <- matrix(0,n, n)
  obs_matrix[D] <- obs
  obs_matrix[cbind(D[,2], D[,1])] <- obs
  
  obs_SS <- SS_calc(obs_matrix)
  
  X_data <- data.frame(site = as.character(1:n), X)
  formula = paste0(colnames(X)[-1])
  
  diss_formula = paste0('~',paste0(colnames(X_data)[-1], collapse = ' + '))
  uniq_formula = paste0('~',paste0(colnames(X_data)[-1], collapse = ' + '))
  uniq_formula_h = paste0(uniq_formula, '+ (1|site)')
  uniq_formula_lcbd = paste0('SS ',uniq_formula)
  
  ##### FULL RE MODEL ####
  if (1%in%mod) {
    mod_h <- gdmm(
      Y_diss = obs,
      X = X_data,
      diss_formula = eval(parse(text = diss_formula)),
      uniq_formula = eval(parse(text = uniq_formula_h)),
      D = D)
    sdrep = TMB::sdreport(mod_h$obj)
    h_sum <- data.frame(
      n = n,
      sd_pair = sd_pair,
      sd_site = sd_site,
      m = m,
      par = names(sdrep$value),
      par_true = c(beta, lambda, beta_0),
      par_recov = sdrep$value,
      par_sd = sdrep$sd
    ) 
  } else {
    h_sum = NULL
  }
  
  if(2%in%mod) {
    ##### BAYESIAN BOOTSTRAPPING ####
    mod_bb <- gdmm(
      Y_diss = obs,
      X = X_data,
      diss_formula = eval(parse(text = diss_formula)),
      uniq_formula = eval(parse(text = uniq_formula)),
      D = D,
      bboot = TRUE,
      n_boot = 1000,
      n_cores = n_cores)
    
    bb_sum <- data.frame(
      n = n,
      sd_pair = sd_pair,
      sd_site = sd_site,
      m = m,
      par = colnames(mod_bb$boot_samples[,-ncol(mod_bb$boot_samples)]),
      par_true = c(beta, lambda, beta_0),
      par_recov = colMeans(mod_bb$boot_samples[,-ncol(mod_bb$boot_samples)]),
      par_low = apply(mod_bb$boot_samples[,-ncol(mod_bb$boot_samples)], 2, function(x) quantile(x, 0.025)),
      par_high = apply(mod_bb$boot_samples[,-ncol(mod_bb$boot_samples)], 2, function(x) quantile(x, 0.975))
    )
  } else {
    bb_sum = NULL
  }

  ##### CLASSIC LCBD MODEL ####
  
  if (3%in%mod) {
    mod_ss <- glmmTMB::glmmTMB(eval(parse(text = uniq_formula_lcbd)), data = cbind(SS = obs_SS, X_data))
    
    ss_sum <- data.frame(
      n = n,
      sd_pair = sd_pair,
      sd_site = sd_site,
      m = m,
      par = c('intercept', rep('lambda', m)),
      par_true = c(beta_0, lambda),
      par_recov = summary(mod_ss)$coefficients$cond[, "Estimate"],
      par_sd = summary(mod_ss)$coefficients$cond[, "Std. Error"]
    )
  } else {
    ss_sum = NULL
  }
  
  return(list(re_model = h_sum, bb_model = bb_sum, ss_model = ss_sum))
}

#### SCENARIOS 1: NORMAL ####
scenarios <-  expand.grid(n = 50, 
                          lambda = 1, # combination of effect size (uniqueness)
                          beta = 1, # combination of effect size (dissimilarity),
                          m = 2, # number of predictors
                          sim = 1:1000,
                          beta_0 = 0,
                          sd_pair = 0.5,
                          sd_site = 0.5,
                          skew_X = c(TRUE, FALSE),
                          rho = c(0, 0.5, 0.8))

# id for different scenarios
scenarios$scenario <- 1:nrow(scenarios)

if(!file.exists('Simulation_study/results_cor.RDS')){
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  results_corr <- foreach(i = 1:nrow(scenarios), .combine = 'c') %dopar% {
    scenario <- scenarios[i,]
    library(gdmmTMB)
    
    sim_gdmm_data(
      n = scenario$n,
      lambda = scenario$lambda,
      beta = scenario$beta,
      m = scenario$m,
      beta_0 = scenario$beta_0,
      sd_pair = scenario$sd_pair,
      sd_site = scenario$sd_site,
      rho = scenario$rho,
      skew_X = scenario$skew_X
    )
      
  }
  
  stopCluster(cl) 
  saveRDS(results_corr, 'Simulation_study/results_cor.RDS')
}

results_corr <- readRDS('Simulation_study/results_cor.RDS')

results_re <- do.call(rbind,results_corr[names(results_corr) == 're_model'])
results_re <- results_re %>% mutate(scenario = rep(1:nrow(scenarios), each = 5)) %>% 
  inner_join(scenarios) %>% mutate(skewed = ifelse(skew_X, 'predictor\nskewed', 'predictor\nnormal'))

results_ss <- do.call(rbind,results_corr[names(results_corr) == 'ss_model'])
results_ss <- results_ss %>% mutate(scenario = rep(1:nrow(scenarios), each = 3)) %>% 
  inner_join(scenarios)%>% mutate(skewed = ifelse(skew_X, 'predictor\nskewed', 'predictor\nnormal'))

results_bb <- do.call(rbind,results_corr[names(results_corr) == 'bb_model'])
results_bb <- results_bb %>% mutate(scenario = rep(1:nrow(scenarios), each = 5)) %>%
  inner_join(scenarios)%>% mutate(skewed = ifelse(skew_X, 'predictor\nskewed', 'predictor\nnormal'))

results_ss$par_recov_scale <- ifelse(results_ss$par == 'lambda', (results_ss$n / (results_ss$n - 2)) * results_ss$par_recov, results_ss$par_recov)
results_ss$par_true_scale <- ifelse(results_ss$par == 'lambda', ((results_ss$n -2) / (results_ss$n)) * results_ss$par_true, results_ss$par_true)
results_ss$par_min <- (results_ss$par_recov - (1.96 * results_ss$par_sd))
results_ss$par_max <- results_ss$par_recov + (1.96 * results_ss$par_sd)
results_ss$coverage <- (results_ss$par_true_scale <= results_ss$par_max) & (results_ss$par_true_scale >= results_ss$par_min)
results_re$par_min = results_re$par_recov - (1.96*results_re$par_sd)
results_re$par_max = results_re$par_recov + (1.96*results_re$par_sd)
results_re$coverage <- (results_re$par_true <= results_re$par_max) & (results_re$par_true >= results_re$par_min)
results_bb$coverage <- (results_bb$par_true <= results_bb$par_high) & (results_bb$par_true >= results_bb$par_low)

#### BIAS ####
ss_BIAS <- results_ss %>%
  mutate(model = 'conventional uniqueness model') %>% 
  mutate(beta = ifelse(beta > 0, 'dissimilarity\ngradient', ' no dissimilarity\ngradient')) %>%
  mutate(BIAS = par_recov_scale - par_true) %>%
  filter(par == 'lambda')

re_BIAS <- results_re %>% 
  mutate(model = 'GDUM (site random effects)') %>% 
  mutate(BIAS = par_recov - par_true) %>%
  mutate(beta = ifelse(beta > 0, 'dissimilarity\ngradient', ' no dissimilarity\ngradient')) %>%
  filter(par == 'lambda')

bb_BIAS <- results_bb %>% 
  mutate(model = 'GDUM (Bayesian bootstrapping)') %>% 
  mutate(BIAS = par_recov - par_true) %>%
  mutate(beta = ifelse(beta > 0, 'dissimilarity\ngradient', ' no dissimilarity\ngradient')) %>%
  filter(par == 'lambda')

BIAS_plot_normal <- ss_BIAS %>% bind_rows(re_BIAS) %>% bind_rows(bb_BIAS) %>% 
  ggplot(aes(x = as.factor(rho), y = BIAS)) +
  geom_boxplot(outlier.size = 0.2,col = 'darkred') +
  geom_hline(yintercept = 0, lty = 5, col = 'grey')+
  facet_grid(skewed~model) +
  theme(legend.position = '',
        legend.key.spacing = unit(0.1, 'cm'),
        legend.background = element_blank()) +
  xlab('correlation') +
  ylab('bias') 

### ERROR ####
ss_RMSE <- results_ss %>%
  mutate(model = 'conventional uniqueness model') %>% 
  group_by(n, par, par_true, beta, model, rho, skewed) %>%
  mutate(beta = ifelse(beta > 0, 'dissimilarity\ngradient', ' no dissimilarity\ngradient')) %>%
  summarise(RMSE = sqrt(sum((par_recov_scale - par_true)^2))) %>%
  filter(par == 'lambda')

re_RMSE <- results_re %>% 
  mutate(model = 'GDUM (site random effects)') %>% 
  group_by(n,par, par_true, beta, model, rho, skewed) %>%
  mutate(beta = ifelse(beta > 0, 'dissimilarity\ngradient', ' no dissimilarity\ngradient')) %>%
  summarise(RMSE = sqrt(sum((par_recov - par_true)^2))) %>%
  filter(par == 'lambda')

bb_RMSE <- results_bb %>% 
  mutate(model = 'GDUM (Bayesian bootstrapping)') %>% 
  group_by(n,par, par_true, beta, model, rho, skewed) %>%
  mutate(beta = ifelse(beta > 0, 'dissimilarity\ngradient', ' no dissimilarity\ngradient')) %>%
  summarise(RMSE = sqrt(sum((par_recov - par_true)^2))) %>%
  filter(par == 'lambda')

RMSE_plot_cor <- ss_RMSE %>% bind_rows(re_RMSE) %>% bind_rows(bb_RMSE) %>% 
  ggplot(aes(x = rho, y = RMSE, colour = as.factor(par_true))) +
  geom_path(col = 'darkred') +
  geom_jitter(width = 0, height = 0, col = 'darkred') +
  facet_grid(skewed~model) +
  xlab('correlation') + 
  theme(legend.position = 'inside',
        legend.position.inside = c(0.95,0.8),
        legend.title.position = 'left',
        legend.key.spacing = unit(0.1, 'cm'),
        legend.background = element_blank()) +
  scale_shape_manual(values = c(21, 24)) +
  scale_linetype_manual(values = c(1,1)) 


## COVERAGE #### 
re_COV <- results_re %>%
  mutate(model = 'GDUM (site random effects)') %>% 
  group_by(rho, skewed, par, par_true, beta, model) %>% 
  mutate(beta = ifelse(beta > 0, 'dissimilarity\ngradient', ' no dissimilarity\ngradient')) %>%
  summarise(coverage = sum(coverage)/n()) %>% 
  filter(par == 'lambda') 

ss_COV <- results_ss %>% 
  mutate(model = 'conventional uniqueness model') %>% 
  group_by(rho, skewed, par, par_true, beta, model) %>% 
  mutate(beta = ifelse(beta > 0, 'dissimilarity\ngradient', ' no dissimilarity\ngradient')) %>%
  summarise(coverage = sum(coverage)/n()) %>% 
  filter(par == 'lambda')

bb_COV <- results_bb %>% 
  mutate(model = 'GDUM (Bayesian bootstrapping)') %>% 
  group_by(rho, skewed, par, par_true, beta, model) %>% 
  mutate(beta = ifelse(beta > 0, 'dissimilarity\ngradient', ' no dissimilarity\ngradient')) %>%
  summarise(coverage = sum(coverage)/n()) %>% 
  filter(par == 'lambda')

COV_plot_cor <-ss_COV %>% bind_rows(re_COV) %>% bind_rows(bb_COV) %>% 
  ggplot(aes(x = rho, y = coverage, col = as.factor(par_true))) +
  geom_jitter(width =0, height = 0, colour = 'darkred') +
  ylim(0,1) +
  xlab('correlation') + 
  facet_grid(skewed~model)+ 
  theme(legend.position = '') +
  ylab('coverage (%)') +
  #geom_hline(yintercept = qbinom(c(0.025, 0.975), size = 1000, prob = 0.95) / 1000, lty = 2, colour = 'grey40') 
  geom_hline(yintercept = 0.95, lty = 5, colour = 'grey')



cowplot::plot_grid(RMSE_plot_cor, BIAS_plot_normal, COV_plot_cor, ncol = 1, align = 'h', labels = c('A', 'B', 'C'))

ggsave('figs/simulation_study_cor.png', width = 15, height = 20, units = 'cm', dpi = 600)


#### COMPUTATION TIME ####
scenarios <-  expand.grid(n = seq(10,100,5), 
                          lambda = 1, # combination of effect size (uniqueness)
                          beta = 1, # combination of effect size (dissimilarity),
                          m = 1, # number of predictors
                          sim = 1,
                          beta_0 = 0,
                          sd_pair = 0.5,
                          sd_site = 0.5,
                          skew_X = FALSE,
                          mod = c(1,2),
                          rho = 0)


scenarios$time = NULL 

for(i in 1:nrow(scenarios)){
  start <- Sys.time()
  
  scenario = scenarios[i,]
  print(scenario$n)
  sim_gdmm_data(
    n = scenario$n,
    lambda = scenario$lambda,
    beta = scenario$beta,
    m = scenario$m,
    beta_0 = scenario$beta_0,
    sd_pair = scenario$sd_pair,
    sd_site = scenario$sd_site,
    rho = scenario$rho,
    skew_X = scenario$skew_X,
    mod = scenario$mod
  )
  scenarios$time[i] <- Sys.time() - start
}


mod = c('GDUM: Site random effects', 'GDUM: Bayesian bootstrapping')
scenarios$model = mod[scenarios$mod]

ggplot(scenarios, aes(x = n, y = time, col = as.factor(model))) + 
  geom_point()+
  geom_line() + 
  ylab('seconds') + 
  xlab('sample size') +
  scale_colour_manual(values = c('darkred', 'pink2'))+
  theme(legend.title = element_blank())

ggsave('figs/comp_time.png', width = 15, height = 10, dpi = 600, units = 'cm')




