
sim_gdmm_data <- function(n = 100, 
                          m = 3, 
                          X_mean = 0,
                          X_sd = NULL,
                          beta_range = c(0, 1),
                          lambda_range = c(-0.5,0.5),
                          sigma_re_range = c(0.001, 0.1),
                          sd_range = c(0.001, 0.1),
                          intercept_range = c(-1,1)) {
  
  if (is.null(X_sd)) X_sd = rep(1, m)
  
  D <- t(combn(n,2))
  X = MASS::mvrnorm(n, rep(0,m), Sigma = Diagonal(m))
  
  # parameters
  real_int = runif(1, intercept_range[1], intercept_range[2])
  real_beta = runif(m, beta_range[1], beta_range[2])
  real_lambda = runif(m, lambda_range[1], lambda_range[2])
  real_sigma_re = runif(1,sigma_re_range[1], sigma_re_range[2])
  real_re <- scale(rnorm(n, 0, 1))*real_sigma_re # ensure mean = 0
  real_sd <- runif(1, sd_range[1], sd_range[2])
  
  # Generate predictions
  diss_comp <- abs(X[D[,1],,drop = F] - X[D[,2],,drop = F]) %*% real_beta
  uniq_comp <- X %*% real_lambda + real_re
  eta <- diss_comp + real_int + uniq_comp[D[,1]] + uniq_comp[D[,2]]
  obs <- eta + scale(rnorm(length(eta))) * real_sd
  
  # Matrix predictions 
  diss_matrix <- matrix(ncol = n, nrow = n)
  diss_matrix[cbind(D[,1], D[,2])] <- obs
  diss_matrix[cbind(D[,2], D[,1])] <- obs
  diag(diss_matrix) <- 0
  Y_SS <- SS_calc(diss_matrix, LCBD = FALSE)
  
  
  # Generate predictions
  X <- data.frame(X, site = as.factor(1:n))
  colnames(X) <- c(paste0('x_', 1:m), 'site')
  
  out <- list(real_param = list(intercept = real_int, 
                                beta = real_beta,
                                lambda = real_lambda),
              Y_diss = obs,
              Y_lcbd = Y_SS,
              X = X,
              D = D)
  return(out)  
}

#############################################################################
############################### with beta lambda ############################
#############################################################################

# Simulation without effects on lambda
coef_list <- list()
for (n in c(30, 100)) {
  for (incl_beta in c(0,1)) {
    for (i in 1:1000) {
      print(i)
      sim_i <- sim_gdmm_data(n = n, 
                             m = 1,
                             beta_range = c(0.1,1)*incl_beta, 
                             sigma_re_range = c(0.1, 0.1),
                             sd_range = c(0.1, 0.1))  
      
      # GDMM
      m_gdmm <- gdmm(Y_diss = sim_i$Y_diss,
                     X = sim_i$X,
                     D = sim_i$D,
                     diss_formula = ~ x_1,
                     uniq_formula = ~ x_1 + (1|site)) 
      
      # BB GDMM
      # m_bbgdmm <- gdmm(Y_diss = sim_i$Y_diss,
      #                  X = sim_i$X,
      #                  D = sim_i$D,
      #                  diss_formula = ~ x_1,
      #                  uniq_formula = ~ x_1,
      #                  bboot = T,
      #                  n_boot = 1000)
      
      # GLMM
      m_lcbd <- glmmTMB::glmmTMB(Y_lcbd ~ x_1, data = cbind(Y_lcbd = sim_i$Y_lcbd, sim_i$X)) 
      
      coef_table_gdmm <- summary(m_gdmm$sdrep)
      para_gdmm <- rownames(coef_table_gdmm)
      coef_table_gdmm <- coef_table_gdmm[para_gdmm %in% c('lambda', 'intercept'),]
      colnames(coef_table_gdmm) <- c('estimate_gdmm', 'se_gdmm')
      
      # sum_bbgdm <- summary(m_bbgdmm, quantiles = c(0.975, 0.025))
      # coef_table_bbgdm <- cbind(sum_bbgdm$estimate, sum_bbgdm$CI)
      # coef_table_bbgdm <- coef_table_bbgdm[rownames(coef_table_bbgdm) %in% c('lambda', 'intercept'),]
      # colnames(coef_table_bbgdm) <- c('estimate_bbgdmm', 'low_bbgdmm', 'high_bbgdm')
      
      coef_table_glmm <- summary(m_lcbd$sdr)
      coef_table_glmm <- coef_table_glmm[rownames(coef_table_glmm) %in% c('beta'),]
      rownames(coef_table_glmm)[1] <- 'intercept'
      colnames(coef_table_glmm) <- c('estimate_lcbd', 'se_lcbd')
      
      coef_list[[length(coef_list) + 1]] <- data.frame(
        real = do.call(c, sim_i$real_param[c('intercept', 'lambda')]),
        para = rownames(coef_table_gdmm),
        coef_table_gdmm,
        coef_table_glmm,
        n = n,
        beta = incl_beta,
       # coef_table_bbgdm,
        sim = i
      )
    }
  }
}




coef_table <- data.frame(do.call(rbind, coef_list))
coef_table$estimate_lcbd <- ifelse(coef_table$para == 'intercept',coef_table$estimate_lcbd*2, coef_table$estimate_lcbd)
coef_table$cover_gdmm <- (coef_table$real >= (coef_table$estimate_gdmm - 1.96*coef_table$se_gdmm)) & (coef_table$real <= (coef_table$estimate_gdmm + 1.96*coef_table$se_gdmm))
coef_table$cover_lcbd <- (coef_table$real >= (coef_table$estimate_lcbd - 1.96*coef_table$se_lcbd)) & (coef_table$real <= (coef_table$estimate_lcbd + 1.96*coef_table$se_lcbd))
coef_table$cover_bbgdm <- (coef_table$real >= coef_table$low_bbgdm) & (coef_table$real <= coef_table$high_bbgdm)

coef_table$bias_gdmm <- coef_table$estimate_gdmm - coef_table$real 
coef_table$bias_lcbd <- coef_table$estimate_lcbd - coef_table$real  
coef_table$bias_bbgdm <- coef_table$estimate_bbgdmm - coef_table$real  
coef_table$n <- paste0('n = ', coef_table$n)
coef_table$beta <- ifelse(coef_table$beta == 1, 'dissimilarity gradient', ' no dissimilarity gradient')


coef_table <- subset(coef_table, para == 'lambda')

rmse_fun <- function(x){
  return(data.frame(y = quantile(x, 0.999), label = paste0("RMSE: ", round(sqrt(mean(x^2)),4))))
}

library(ggplot2)
theme_set(theme_bw())
theme_update(strip.background = element_blank(), 
             legend.position = '',
             panel.grid = element_blank())
plot_lcbd <- ggplot(coef_table, aes(x = as.factor(n), y = abs(bias_lcbd), colour = 'LCBD model')) +
  geom_jitter(width = 0.2,alpha = 0.1, size = 0.5) + 
 # geom_boxplot(aes(colour = as.factor(n))) + 
  stat_summary(colour = 'black', fun = mean, fun.max = function(x) quantile(x, 0.975), fun.min = function(x) quantile(x, 0.25),  size = 0.2)  + 
  facet_wrap(~as.factor(beta), ncol = 1) + 
  ylab('|bias|') +
  xlab('') + 
  ylim(0, max(coef_table$bias_gdmm, coef_table$bias_lcbd)) + 
  stat_summary(col = 'black', fun.data = rmse_fun, geom = "text", size = 2) + 
  scale_colour_manual(values = c('orange'))

plot_gdmm <- ggplot(coef_table, aes(x = as.factor(n), y = abs(bias_gdmm), colour = 'LCBD model')) +
  geom_jitter(width = 0.2,alpha = 0.1, size = 0.5) + 
  # geom_boxplot(aes(colour = as.factor(n))) + 
  stat_summary(colour = 'black', fun = mean, fun.max = function(x) quantile(x, 0.975), fun.min = function(x) quantile(x, 0.25), size = 0.2)  + 
  facet_wrap(~as.factor(beta), ncol = 1) + 
  ylab('|bias|') +
  xlab('') + 
  guides(alpha = 'none') + 
  ylim(0, max(coef_table$bias_gdmm, coef_table$bias_lcbd)) + 
  stat_summary(col = 'black', fun.data = rmse_fun, geom = "text", size = 2) + 
  scale_colour_manual(values = c('steelblue'))

library(cowplot)

coef_match <- ggplot(coef_table, aes(x = real, y = estimate_gdmm)) + 
  geom_linerange(aes(colour = cover_gdmm, ymin = estimate_gdmm - 1.96*se_gdmm, ymax = estimate_gdmm + 1.96*se_gdmm ), alpha = 0.2) +
  geom_point(aes(colour = cover_gdmm), size = 0.5) +  
  ylab('lambda (recovered)') +
  xlab('lambda (simulated)') +
  geom_text(data = NULL, aes(x = -0.3, y = 0.5, label = 'coverage: 94.35%')) + 
  geom_abline(slope = 1, intercept = 0, colour = 'grey') + 
  scale_colour_manual(values = c('red', 'black'))

plot_grid(plot_grid(plot_lcbd, plot_gdmm, labels = c('A', 'B')), coef_match, labels = c('', 'C'))

ggsave('Proposal/sim.png', dpi = 600, units = 'cm', width = 19, height = 10)



