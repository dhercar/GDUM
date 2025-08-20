# Author: Daniel Hernaandez Carrasco
# Email: dani.hc97@gmail.com
# Created: 8/13/2025
# License: MIT (see LICENSE file for details)

# Description: This script performs a simulation study comparing different approaches used to model xommunity uniqueness

library(gdmmTMB)
library(doParallel)
library(ggplot2)
theme_set(theme_bw())
theme_update(strip.background = element_blank(),
             text = element_text(size =7),
             panel.grid = element_blank())

sim_gdmm_data <- function(n = NULL, # sample size
                          lambda = NULL, # effect size (uniqueness)
                          beta = NULL, # effect size (dissimilarity),
                          include_beta = FALSE,
                          m = NULL,
                          beta_0 = NULL,
                          sd_pair = NULL,
                          sd_site= NULL, 
                          skew_X = FALSE) {
  
  beta = rep(beta, m)
  lambda = rep(lambda, m)

  D <- t(combn(n,2))
  n_pair = nrow(D)
  
  re_site = scale(rnorm(n))*sd_site
  err_pair = scale(rnorm(n_pair))*sd_pair

  X = scale(MASS::mvrnorm(n = n, mu = rep(0,m), Sigma = as.matrix(Matrix::Diagonal(m,1))))
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
  
  ##### BAYESIAN BOOTSTRAPPING ####
  mod_bb <- gdmm(
    Y_diss = obs,
    X = X_data,
    diss_formula = eval(parse(text = diss_formula)),
    uniq_formula = eval(parse(text = uniq_formula)),
    D = D,
    bboot = TRUE,
    n_boot = 1000,
    n_cores = 4)

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

  
  ##### CLASSIC LCBD MODEL ####
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
  return(list(re_model = h_sum, bb_model = bb_sum, ss_model = ss_sum))
}


scenarios <-  expand.grid(n = c(20, 50, 100), 
                          lambda = c(0, 0.5, 1), # effect size (uniqueness)
                          beta = c(0, 1), # effect size (dissimilarity),
                          m = 1,
                          sim = 1:1000,
                          beta_0 = 0,
                          sd_pair = 0.5,
                          sd_site = 0.5) 

scenarios$scenario <- 1:nrow(scenarios)

library(doParallel)

# Set up cluster
cl <- makeCluster(5)
registerDoParallel(cl)

results_normal <- foreach(i = 1:nrow(scenarios), .combine = 'c') %dopar% {
  scenario <- scenarios[i,]
  library(gdmmTMB)
  sim_gdmm_data(
    n = scenario$n,
    lambda = scenario$lambda,
    beta = scenario$beta,
    m = scenario$m,
    beta_0 = scenario$beta_0,
    sd_pair = scenario$sd_pair,
    sd_site = scenario$sd_site
  )
}

stopCluster(cl)
saveRDS(results_normal, 'Simulation_study/results.RDS')

cl <- makeCluster(5)
registerDoParallel(cl)

results_skewed <- foreach(i = 1:nrow(scenarios), .combine = 'c') %dopar% {
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
    skew_X = TRUE
  )
}

saveRDS(results_skewed, 'Simulation_study/results_skewed.RDS')

stopCluster(cl)

results_re <- do.call(rbind,results_normal[names(results_normal) == 're_model'])
results_re <- results_re %>% mutate(scenario = rep(1:nrow(scenarios), each = 3)) %>% inner_join(scenarios)

results_ss <- do.call(rbind,results_normal[names(results_normal) == 'ss_model'])
results_ss <- results_ss %>% mutate(scenario = rep(1:nrow(scenarios), each = 2)) %>% inner_join(scenarios)


results_bb <- do.call(rbind,results_normal[names(results_normal) == 'bb_model'])
results_bb <- results_bb %>% mutate(scenario = rep(1:nrow(scenarios), each = 3)) %>% inner_join(scenarios)


results_ss$par_recov_scale <- ifelse(results_ss$par == 'lambda', (results_ss$n / (results_ss$n - 2)) * results_ss$par_recov, results_ss$par_recov)
results_ss$par_true_scale <- ifelse(results_ss$par == 'lambda', ((results_ss$n -2) / (results_ss$n)) * results_ss$par_true, results_ss$par_true)
results_ss$par_min <- (results_ss$par_recov - (1.96 * results_ss$par_sd))
results_ss$par_max <- results_ss$par_recov + (1.96 * results_ss$par_sd)
results_ss$coverage <- (results_ss$par_true_scale <= results_ss$par_max) & (results_ss$par_true_scale >= results_ss$par_min)

results_re$par_min = results_re$par_recov - (1.96*results_re$par_sd)
results_re$par_max = results_re$par_recov + (1.96*results_re$par_sd)
results_re$coverage <- (results_re$par_true <= results_re$par_max) & (results_re$par_true >= results_re$par_min)

results_bb$coverage <- (results_bb$par_true <= results_bb$par_high) & (results_bb$par_true >= results_bb$par_low)


##### FULL COMBINATIONS ####

#### BIAS ####
ss_BIAS <- results_ss %>%
  mutate(model = ' Site uniqueness') %>% 
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  mutate(BIAS = par_recov_scale - par_true) %>%
  filter(par == 'lambda')

re_BIAS <- results_re %>% 
  mutate(model = 'GDUM\n(site random effects)') %>% 
  mutate(BIAS = par_recov - par_true) %>%
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  filter(par == 'lambda')

bb_BIAS <- results_bb %>% 
  mutate(model = 'GDUM\n(Bayesian bootstrapping)') %>% 
  mutate(BIAS = par_recov - par_true) %>%
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  filter(par == 'lambda')

BIAS_plot_normal <- ss_BIAS %>% bind_rows(re_BIAS) %>% bind_rows(bb_BIAS) %>% 
  ggplot(aes(x = as.factor(n), y = BIAS, colour = as.factor(par_true))) +
  geom_boxplot(outlier.size = 0.2) +
  geom_hline(yintercept = 0, lty = 5, col = 'red')+
  facet_grid(beta~model) +
  theme(legend.position = '',
        legend.key.spacing = unit(0.1, 'cm'),
        legend.background = element_blank()) +
  xlab('n') +
  ylim(-1, 1) +
  scale_colour_manual(expression(lambda), values = c('steelblue', 'khaki3', 'orange'), labels = c(0, 0.5, 1))   

(BIAS_plot_normal_sub <- ss_BIAS %>% 
  bind_rows(re_BIAS) %>% 
  bind_rows(bb_BIAS) %>% 
  filter(n %in% c(20, 100), lambda == 1) %>%
  ggplot(aes(x = as.factor(n), 
             y = BIAS, 
             fill = as.factor(model))) +
  geom_hline(yintercept = 0, lty = 5, col = 'grey') +
  stat_summary(aes(lty = beta, shape = beta), 
               geom = 'pointrange',
               fun.data = mean_sdl,
               position = position_dodge(width = 0.7)) +
  theme(legend.position = '',
        legend.key.spacing = unit(0.1, 'cm'),
        legend.background = element_blank()) +
  xlab('n') +
  ylim(-1,1) + 
  scale_shape_manual(values = c(21, 24)) +
  scale_linetype_manual(values = c(1,5)) +
  scale_fill_manual('model', 
                      values = c('steelblue', 'khaki3', 'coral'), 
                      labels = c(0, 0.5, 1)))

### ERROR ####
ss_RMSE <- results_ss %>%
  mutate(model = ' Site uniqueness') %>% 
  group_by(n, par, par_true, beta, model) %>%
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  summarise(RMSE = sqrt(sum((par_recov_scale - par_true)^2))) %>%
  filter(par == 'lambda')

re_RMSE <- results_re %>% 
  mutate(model = 'GDUM\n(site random effects)') %>% 
  group_by(n,par, par_true, beta, model) %>%
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  summarise(RMSE = sqrt(sum((par_recov - par_true)^2))) %>%
  filter(par == 'lambda')

bb_RMSE <- results_bb %>% 
  mutate(model = 'GDUM\n(Bayesian bootstrapping)') %>% 
  group_by(n,par, par_true, beta, model) %>%
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  summarise(RMSE = sqrt(sum((par_recov - par_true)^2))) %>%
  filter(par == 'lambda')

RMSE_plot_normal <- ss_RMSE %>% bind_rows(re_RMSE) %>% bind_rows(bb_RMSE) %>% 
  ggplot(aes(x = n, y = RMSE, colour = as.factor(par_true))) +
  geom_path() +
  geom_jitter(width = 0.1, height = 0) +
  facet_grid(beta~model) +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.95,0.8),
        legend.title.position = 'left',
        legend.key.spacing = unit(0.1, 'cm'),
        legend.background = element_blank()) +
  ylim(1,11.5) + 
  scale_shape_manual(values = c(21, 24)) +
  scale_linetype_manual(values = c(1,5)) +
  scale_colour_manual(expression(lambda), 
                    values = c('steelblue', 'khaki3', 'coral'), 
                    labels = c(0, 0.5, 1))


(RMSE_plot_normal_sub <- ss_RMSE %>% 
    bind_rows(re_RMSE) %>% 
    bind_rows(bb_RMSE) %>% 
    filter(n %in% c(20, 100), par_true == 1) %>%
    ggplot(aes(x = as.factor(n), y = RMSE, colour = model)) +
    geom_hline(yintercept = 0, colour = 'grey50') +
    geom_linerange(aes(lty = beta, ymin = 0, ymax = RMSE), position = position_dodge(width = 0.6)) +
    geom_point(aes(shape = beta), position = position_dodge(width = 0.6)) +
    theme(legend.position = 'inside',
          legend.position.inside = c(0.95,0.8),
          legend.title.position = 'left',
          legend.key.spacing = unit(0.1, 'cm'),
          legend.background = element_blank()) +
    ylim(0,11.5) + 
    scale_shape_manual(values = c(16, 17)) +
    scale_linetype_manual(values = c(1,5)) +
    scale_colour_manual('model', 
                      values = c('steelblue', 'khaki3', 'coral'), 
                      labels = c(0, 0.5, 1)))


## COVERAGE #### COVERAGE ##RMSE_plot_normal
re_COV <- results_re %>%
  mutate(model = 'GDUM\n(site random effects)') %>% 
  group_by(n, par, par_true, beta, model) %>% 
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  summarise(coverage = sum(coverage)/n()) %>% 
  filter(par == 'lambda') 

ss_COV <- results_ss %>% 
  mutate(model = ' Site uniqueness') %>% 
  group_by(n, par, par_true, beta, model) %>% 
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  summarise(coverage = sum(coverage)/n()) %>% 
  filter(par == 'lambda')

bb_COV <- results_bb %>% 
  mutate(model = 'GDUM\n(Bayesian bootstrapping)') %>% 
  group_by(n, par, par_true, beta, model) %>% 
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  summarise(coverage = sum(coverage)/n()) %>% 
  filter(par == 'lambda')

COV_plot_normal <-ss_COV %>% bind_rows(re_COV) %>% bind_rows(bb_COV) %>% 
  ggplot(aes(x = n, y = coverage, col = as.factor(par_true))) +
  geom_jitter(width =1, height = 0) +
  ylim(0,1) +
  facet_grid(beta~model)+ 
  theme(legend.position = '') +
  ylab('coverage (%)') +
  scale_colour_manual(expression(lambda), values = c('steelblue', 'khaki3', 'orange'), labels = c(0, 0.5, 1)) +
  #geom_hline(yintercept = qbinom(c(0.025, 0.975), size = 1000, prob = 0.95) / 1000, lty = 2, colour = 'grey40') 
  geom_hline(yintercept = 0.95, lty = 2, colour = 'grey')


(COV_plot_normal_sub <- ss_COV %>% 
    bind_rows(re_COV) %>% 
    bind_rows(bb_COV) %>% 
    filter(n %in% c(20, 100), par_true == 1) %>%
    ggplot(aes(x = as.factor(n), y = coverage, colour = model)) +
    geom_hline(yintercept = 0.95, colour = 'grey50') +
    geom_point(aes(shape = beta), position = position_dodge(width = 0.6)) +
    theme(legend.position = 'inside',
          legend.position.inside = c(0.95,0.8),
          legend.title.position = 'left',
          legend.key.spacing = unit(0.1, 'cm'),
          legend.background = element_blank()) +
    ylim(0,1) + 
    scale_shape_manual(values = c(16, 17)) +
    scale_colour_manual('model', 
                        values = c('steelblue', 'khaki3', 'coral'), 
                        labels = c(0, 0.5, 1))
)


cowplot::plot_grid(RMSE_plot_normal, BIAS_plot_normal, COV_plot_normal, ncol = 1, align = 'h', labels = c('A', 'B', 'C'))
ggsave('figs/simulation_study.png', width = 15, height = 20, units = 'cm', dpi = 600)


#### SKEWED DATA ###

results_re <- do.call(rbind,results_skewed[names(results_skewed) == 're_model'])
results_re <- results_re %>% mutate(scenario = rep(1:nrow(scenarios), each = 3)) %>% inner_join(scenarios)

results_ss <- do.call(rbind,results_skewed[names(results_skewed) == 'ss_model'])
results_ss <- results_ss %>% mutate(scenario = rep(1:nrow(scenarios), each = 2)) %>% inner_join(scenarios)

results_bb <- do.call(rbind,results_skewed[names(results_skewed) == 'bb_model'])
results_bb <- results_bb %>% mutate(scenario = rep(1:nrow(scenarios), each = 3)) %>% inner_join(scenarios)


results_ss$par_recov_scale <- ifelse(results_ss$par == 'lambda', (results_ss$n / (results_ss$n - 2)) * results_ss$par_recov, results_ss$par_recov)
results_ss$par_true_scale <- ifelse(results_ss$par == 'lambda', ((results_ss$n -2) / (results_ss$n)) * results_ss$par_true, results_ss$par_true)
results_ss$par_min <- (results_ss$par_recov - (1.96 * results_ss$par_sd))
results_ss$par_max <- results_ss$par_recov + (1.96 * results_ss$par_sd)
results_ss$coverage <- (results_ss$par_true_scale <= results_ss$par_max) & (results_ss$par_true_scale >= results_ss$par_min)

results_re$par_min = results_re$par_recov - (1.96*results_re$par_sd)
results_re$par_max = results_re$par_recov + (1.96*results_re$par_sd)
results_re$coverage <- (results_re$par_true <= results_re$par_max) & (results_re$par_true >= results_re$par_min)

results_bb$coverage <- (results_bb$par_true <= results_bb$par_high) & (results_bb$par_true >= results_bb$par_low)



### BIAS ####
ss_BIAS <- results_ss %>%
  mutate(model = ' Site uniqueness') %>% 
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  mutate(BIAS = par_recov_scale - par_true) %>%
  filter(par == 'lambda')

re_BIAS <- results_re %>% 
  mutate(model = 'GDUM\n(site random effects)') %>% 
  mutate(BIAS = par_recov - par_true) %>%
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  filter(par == 'lambda')

bb_BIAS <- results_bb %>% 
  mutate(model = 'GDUM\n(Bayesian bootstrapping)') %>% 
  mutate(BIAS = par_recov - par_true) %>%
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  filter(par == 'lambda')

BIAS_plot_skewed <- ss_BIAS %>% bind_rows(re_BIAS) %>% bind_rows(bb_BIAS) %>% 
  ggplot(aes(x = as.factor(n), y = BIAS, colour = as.factor(par_true))) +
  geom_boxplot(outlier.size = 0.2) +
  geom_hline(yintercept = 0, lty = 5, col = 'red')+
  facet_grid(beta~model) +
  theme(legend.position = '',
        legend.key.spacing = unit(0.1, 'cm'),
        legend.background = element_blank()) +
  xlab('n') +
  ylim(-1, 1) +
  scale_colour_manual(expression(lambda), values = c('steelblue', 'khaki3', 'orange'), labels = c(0, 0.5, 1))   

(BIAS_plot_skewed_sub <- ss_BIAS %>% 
    bind_rows(re_BIAS) %>% 
    bind_rows(bb_BIAS) %>% 
    filter(n %in% c(20, 100), lambda == 1) %>%
    ggplot(aes(x = as.factor(n), 
               y = BIAS, 
               fill = as.factor(model))) +
    geom_hline(yintercept = 0, lty = 5, col = 'grey') +
    stat_summary(aes(lty = beta, shape = beta), 
                 geom = 'pointrange',
                 fun.data = mean_sdl,
                 position = position_dodge(width = 0.7)) +
    theme(legend.position = '',
          legend.key.spacing = unit(0.1, 'cm'),
          legend.background = element_blank()) +
    xlab('n') +
    ylim(-1,1) + 
    scale_shape_manual(values = c(21, 24)) +
    scale_linetype_manual(values = c(1,5)) +
    scale_fill_manual('model', 
                      values = c('steelblue', 'khaki3', 'coral'), 
                      labels = c(0, 0.5, 1)))

### ERROR ####
ss_RMSE <- results_ss %>%
  mutate(model = ' Site uniqueness') %>% 
  group_by(n, par, par_true, beta, model) %>%
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  summarise(RMSE = sqrt(sum((par_recov_scale - par_true)^2))) %>%
  filter(par == 'lambda')

re_RMSE <- results_re %>% 
  mutate(model = 'GDUM\n(site random effects)') %>% 
  group_by(n,par, par_true, beta, model) %>%
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  summarise(RMSE = sqrt(sum((par_recov - par_true)^2))) %>%
  filter(par == 'lambda')

bb_RMSE <- results_bb %>% 
  mutate(model = 'GDUM\n(Bayesian bootstrapping)') %>% 
  group_by(n,par, par_true, beta, model) %>%
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  summarise(RMSE = sqrt(sum((par_recov - par_true)^2))) %>%
  filter(par == 'lambda')

RMSE_plot_skewed <- ss_RMSE %>% bind_rows(re_RMSE) %>% bind_rows(bb_RMSE) %>% 
  ggplot(aes(x = n, y = RMSE, colour = as.factor(par_true))) +
  geom_path() +
  geom_jitter(width = 0.1, height = 0) +
  facet_grid(beta~model) +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.95,0.8),
        legend.title.position = 'left',
        legend.key.spacing = unit(0.1, 'cm'),
        legend.background = element_blank()) +
  ylim(1,11.5) + 
  scale_shape_manual(values = c(21, 24)) +
  scale_linetype_manual(values = c(1,5)) +
  scale_colour_manual(expression(lambda), 
                      values = c('steelblue', 'khaki3', 'coral'), 
                      labels = c(0, 0.5, 1))


(RMSE_plot_skewed_sub <- ss_RMSE %>% 
    bind_rows(re_RMSE) %>% 
    bind_rows(bb_RMSE) %>% 
    filter(n %in% c(20, 100), par_true == 1) %>%
    ggplot(aes(x = as.factor(n), y = RMSE, colour = model)) +
    geom_hline(yintercept = 0, colour = 'grey50') +
    geom_linerange(aes(lty = beta, ymin = 0, ymax = RMSE), position = position_dodge(width = 0.6)) +
    geom_point(aes(shape = beta), position = position_dodge(width = 0.6)) +
    theme(legend.position = 'inside',
          legend.position.inside = c(0.95,0.8),
          legend.title.position = 'left',
          legend.key.spacing = unit(0.1, 'cm'),
          legend.background = element_blank()) +
    ylim(0,11.5) + 
    scale_shape_manual(values = c(16, 17)) +
    scale_linetype_manual(values = c(1,5)) +
    scale_colour_manual('model', 
                        values = c('steelblue', 'khaki3', 'coral'), 
                        labels = c(0, 0.5, 1)))


## COVERAGE ##
re_COV <- results_re %>%
  mutate(model = 'GDUM\n(site random effects)') %>% 
  group_by(n, par, par_true, beta, model) %>% 
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  summarise(coverage = sum(coverage)/n()) %>% 
  filter(par == 'lambda') 

ss_COV <- results_ss %>% 
  mutate(model = ' Site uniqueness') %>% 
  group_by(n, par, par_true, beta, model) %>% 
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  summarise(coverage = sum(coverage)/n()) %>% 
  filter(par == 'lambda')

bb_COV <- results_bb %>% 
  mutate(model = 'GDUM\n(Bayesian bootstrapping)') %>% 
  group_by(n, par, par_true, beta, model) %>% 
  mutate(beta = ifelse(beta > 0, 'dissimilarity gradient', ' no dissimilarity gradient')) %>%
  summarise(coverage = sum(coverage)/n()) %>% 
  filter(par == 'lambda')

COV_plot_skewed <-ss_COV %>% bind_rows(re_COV) %>% bind_rows(bb_COV) %>% 
  ggplot(aes(x = n, y = coverage, col = as.factor(par_true))) +
  geom_jitter(width =1, height = 0) +
  ylim(0,1) +
  facet_grid(beta~model)+ 
  theme(legend.position = '') +
  ylab('coverage (%)') +
  scale_colour_manual(expression(lambda), values = c('steelblue', 'khaki3', 'orange'), labels = c(0, 0.5, 1)) +
  #geom_hline(yintercept = qbinom(c(0.025, 0.975), size = 1000, prob = 0.95) / 1000, lty = 2, colour = 'grey40') 
  geom_hline(yintercept = 0.95, lty = 2, colour = 'grey')


(COV_plot_skewed_sub <- ss_COV %>% 
    bind_rows(re_COV) %>% 
    bind_rows(bb_COV) %>% 
    filter(n %in% c(20, 100), par_true == 1) %>%
    ggplot(aes(x = as.factor(n), y = coverage, colour = model)) +
    geom_hline(yintercept = 0.95, colour = 'grey50') +
    geom_point(aes(shape = beta), position = position_dodge(width = 0.6)) +
    theme(legend.position = 'inside',
          legend.position.inside = c(0.95,0.8),
          legend.title.position = 'left',
          legend.key.spacing = unit(0.1, 'cm'),
          legend.background = element_blank()) +
    ylim(0,1) + 
    scale_shape_manual(values = c(16, 17)) +
    scale_colour_manual('model', 
                        values = c('steelblue', 'khaki3', 'coral'), 
                        labels = c(0, 0.5, 1))
)



cowplot::plot_grid(RMSE_plot_skewed, BIAS_plot_skewed, COV_plot_skewed, ncol = 1, align = 'h', labels = c('A', 'B', 'C'))
ggsave('figs/simulation_study_skew.png', width = 15, height = 20, units = 'cm', dpi = 600)


set.seed(123)
(dist_normal <- data.frame(x = scale(rnorm(10000))) %>% 
    ggplot(aes(x = x)) + geom_histogram(fill = 'grey80', col = 'grey50') + 
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()) + 
    xlab('w') +
    xlim(-6,4) +
    scale_y_continuous('frequency', expand = c(0,0)))
(dist_skewed <- data.frame(x = scale(rnorm(10000))) %>%
    ggplot(aes(x = scale(log1p(x - min(x))))) + 
    geom_histogram(fill = 'grey80', col = 'grey50') + 
    theme_classic() +
    xlab('w') +
    xlim(-6,4) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()) + 
    scale_y_continuous('frequency', expand = c(0,0)))


cowplot::plot_grid(dist_normal, dist_skewed,
                   RMSE_plot_normal_sub, RMSE_plot_skewed_sub,
                   BIAS_plot_normal_sub, BIAS_plot_skewed_sub,
                   ncol = 2, 
                   align = 'h', 
                   rel_heights = c(2,3,3),
                   labels = c('A', 'B', 'C', 'D', 'E', 'F'))
ggsave('figs/simulation_study_skew.png', width = 15, height = 20, units = 'cm', dpi = 600)

### ERROR BETA ####

re_RMSE <- results_re %>% 
  mutate(model = 'GDUM (site random effects)') %>% 
  group_by(n,par, par_true, beta, model, lambda) %>%
  summarise(RMSE = sqrt(sum((par_recov - par_true)^2))) %>%
  filter(par == 'e_beta')

bb_RMSE <- results_bb %>% 
  mutate(model = 'GDUM (Bayesian bootstrapping)') %>% 
  group_by(n,par, par_true, beta, model, lambda) %>%
  summarise(RMSE = sqrt(sum((par_recov - par_true)^2))) %>%
  filter(par == 'e_beta')

RMSE_plot <- re_RMSE %>% bind_rows(bb_RMSE) %>% 
  ggplot(aes(x = n, y = RMSE, colour = as.factor(par_true))) +
  geom_path() +
  geom_jitter(width =1, height = 0) +
  facet_grid(lambda~model) +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.1,0.95),
        legend.background = element_blank()) +
  scale_colour_manual('', values = c('steelblue4', 'orange'), labels = paste0('\03bb = ', c(0,1))) 




