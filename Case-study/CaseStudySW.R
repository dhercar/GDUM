# Author: Daniel Hernaandez Carrasco
# Email: dani.hc97@gmail.com
# Created: 8/13/2025
# License: MIT (see LICENSE file for details)

library(tidyverse)
library(gdmmTMB)
library(ggplot2)
library(gdm)
library(adespatial)

# plotting preferences
theme_set(theme_bw())
theme_update(strip.background = element_blank(), 
             text = element_text(size = 8),
             panel.grid = element_blank())

# data
sp <- gdm::southwest %>% dplyr::select(site, species) %>% 
  group_by(site, species) %>%
  summarise(n = 1) %>%
  pivot_wider(id_cols = site, names_from = species,values_from = n, values_fill = 0 ) %>%
  arrange(site)

env <- gdm::southwest %>% dplyr::select(-c(species)) %>%
  distinct() %>% arrange(site)


D <- t(combn(nrow(env), 2))
dist = data.frame(dist = scale(dist(env[,c('bio19')]))[D])

# model with dissimilarity gradients
m2.d <- gdmm(Y = sp[,-1],
          X = env,
          X_pair = dist,
          D = D,
          diss_formula = ~ isp(bio5) + isp(bio6) + isp(bio15) + isp(log(bio19)),
          uniq_formula = ~ scale(bio5) + scale(bio6) + scale(bio15) + scale(log(bio19)),
          pair_formula = ~ dist,
          mono = T,
          mono_pair = T,
          family = 'binomial',
          n_boot = 1000,
          bboot = T)

# model without dissimilarity gradients
m2.nd <- gdmm(Y = sp[,-1],
             X = env,
             X_pair = dist,
             D = D,
            # diss_formula = ~ isp(bio5) + isp(bio6) + isp(bio15) + isp(log(bio19)),
             uniq_formula = ~ scale(bio5) + scale(bio6) + scale(bio15) + scale(log(bio19)),
            # pair_formula = ~ dist,
             mono = T,
             mono_pair = T,
             family = 'binomial',
             n_boot = 1000,
             bboot = T)

#### DISSIMILARITY GRADIENTS INCLUDED ####
pred_diss_d <- data.frame(predict(m2.d))
pred_diss_d$obs <- m2.nd$Y_diss/m2.d$Y_den

pred_lcbd_d <- data.frame(predict(m2.d, component = 'uniqueness', scale_uniq = T))
pred_lcbd_d$obs <- adespatial::LCBD.comp(vegan::vegdist(sp[,-1]))$LCBD

R2_diss_d <- cor(pred_diss_d$mean, pred_diss_d$obs)^2
R2_lcbd_d <- cor(pred_lcbd_d$mean, pred_lcbd_d$obs)^2

pred_diss_plot_d <- ggplot(pred_diss_d, aes(x = mean, y = obs)) + 
  geom_point(col = 'darkred', size = 1) +
  geom_abline(slope = 1) +
  coord_equal() +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, face =2)) +
  ylim(0.1,1) +
  xlim(0.1,1) +
  xlab('predicted') +
  ylab('observed')+
  ggtitle('dissimilarity') + 
  annotate("text", x = 0.1, y = 0.95, label = paste0("R² = ", round(R2_diss_d, 2)), hjust = 0)

pred_lcbd_plot_d <- ggplot(pred_lcbd_d, aes(x = mean, y = obs)) + 
  geom_point(col = 'darkred', size = 1) +
  coord_equal() +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, face =2)) +
  geom_abline(slope = 1) +
  ylim(0.0083,0.0166) +
  xlim(0.0083,0.015)+
  xlab('predicted') +
  ylab('observed') +
  ggtitle('LCBD')+ 
  annotate("text", x = 0.0083, y = 0.016, label = paste0("R² = ", round(R2_lcbd_d, 2)), hjust = 0)


f_x <- diss_gradient(m2.d, CI_quant = c(0.5, 0.95))

env_long <- env %>% dplyr::select(bio5, bio6, bio15, bio19) %>% pivot_longer(c(bio5, bio6, bio15, bio19), names_to = 'var') 

fx_plot <- do.call(rbind, f_x) %>% ggplot() +
  geom_ribbon(aes(x = x, ymax = `CI 2.5%`, ymin = `CI 97.5%`), fill = 'grey90') +
  geom_ribbon(aes(x = x, ymax = `CI 25%`, ymin = `CI 75%`), fill = 'grey80') +
  facet_wrap(~var, scales = 'free_x', ncol = 1, strip.position = 'bottom') +
  theme(strip.placement = 'outside') +
  ylab(expression('f(x)')) +
  xlab('')+
  geom_line(aes(x = x, y = f_x)) +
  geom_point(data = env_long, aes(x = value, y = -0.2), shape = '|', col = 'darkred')
  
lambda_samples_d <- m2.d$boot_samples[,colnames(m2.d$boot_samples) == 'lambda']
colnames(lambda_samples_d) <- names(m2.d$form_W$predictors)

lambda_sum_d <- data.frame(t(apply(lambda_samples_d, 2, function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))))
lambda_sum_d$var <- c('bio5', 'bio6', 'bio15','bio19')
lambda_plot_d <- ggplot(lambda_sum_d, aes(x = X50., y = var)) +
  geom_vline(xintercept = 0) +
  ylab('') +
  xlab(expression(lambda)) +
  xlim(-0.4, 0.4) +
  scale_y_discrete(limits = rev) +
  geom_linerange(aes(xmin = X2.5., xmax = X97.5.), col = 'grey', linewidth = 1) +
  geom_linerange(aes(xmin = X25., xmax = X75.), linewidth = 1) + 
  geom_point(colour = 'black', fill = 'white', shape = 21, size = 3, stroke = 1) 

cowplot::plot_grid(fx_plot, lambda_plot_d, cowplot::plot_grid(pred_diss_plot_d, pred_lcbd_plot_d, ncol = 1, align = 'vh', labels = c('C', 'D')), ncol = 3, labels = c('A', 'B'), rel_widths = c(1,1,1))

ggsave('figs/cs_2_plot.png', width = 15, height = 15, units = 'cm', dpi = 600)
ggsave('figs/cs_2_plot.eps', device = cairo_ps, width = 15, height = 15, units = 'cm', dpi = 600)


#### DISSIMILARITY GRADIENTS NOT INCLUDED ####
pred_diss_nd <- data.frame(predict(m2.nd))
pred_diss_nd$obs <- m2.nd$Y_diss/m2.nd$Y_den

pred_lcbd_nd <- data.frame(predict(m2.nd, component = 'uniqueness', scale_uniq = T))
pred_lcbd_nd$obs <- adespatial::LCBD.comp(vegdist(sp[,-1]))$LCBD

R2_diss_nd <- cor(pred_diss_nd$mean, pred_diss_nd$obs)^2
R2_lcbd_nd <- cor(pred_lcbd_nd$mean, pred_lcbd_nd$obs)^2

pred_diss_plot_nd <- ggplot(pred_diss_nd, aes(x = mean, y = obs)) + 
  geom_point(col = 'darkred', alpha = 0.3) +
  geom_abline(slope = 1) +
  coord_equal() +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, face =2)) +
  ylim(0.1,1) +
  xlim(0.1,1) +
  xlab('predicted') +
  ylab('observed')+
  ggtitle('dissimilarity') + 
  annotate("text", x = 0.1, y = 0.95, label = paste0("R² = ", round(R2_diss_nd, 2)), hjust = 0)

pred_lcbd_plot_nd <- ggplot(pred_lcbd_nd, aes(x = mean, y = obs)) + 
  geom_point(col = 'darkred', alpha = 0.5) +
  coord_equal() +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, face =2)) +
  geom_abline(slope = 1) +
  ylim(0.0083,0.0166) +
  xlim(0.0083,0.015)+
  xlab('predicted') +
  ylab('observed') +
  ggtitle('LCBD')+ 
  annotate("text", x = 0.0083, y = 0.016, label = paste0("R² = ", round(R2_lcbd_nd, 2)), hjust = 0)


lambda_samples_nd <- m2.nd$boot_samples[,colnames(m2.nd$boot_samples) == 'lambda']
colnames(lambda_samples_nd) <- names(m2.nd$form_W$predictors)

lambda_sum_nd <- data.frame(t(apply(lambda_samples_nd, 2, function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))))
lambda_sum_nd$var <- c('bio5', 'bio6', 'bio15','bio19')

lambda_plot_nd <- ggplot(lambda_sum_nd, aes(x = X50., y = var)) +
  geom_vline(xintercept = 0) +
  ylab('') +
  xlab(expression(lambda)) +
  xlim(-0.7, 0.6) +
  scale_y_discrete(limits = rev) +
  geom_linerange(aes(xmin = X2.5., xmax = X97.5.), col = 'grey', linewidth = 1) +
  geom_linerange(aes(xmin = X25., xmax = X75.), linewidth = 1) + 
  geom_point(colour = 'black', fill = 'white', shape = 21, size = 3, stroke = 1) 

cowplot::plot_grid(cowplot::plot_grid(lambda_plot_nd+ xlim(-0.7,0.5),
                                      lambda_plot_d + xlim(-0.7,0.5),
                                      ncol = 2, labels = c('A', 'B'), align = 'vh'),
                   cowplot::plot_grid(pred_diss_plot_nd, 
                                      pred_diss_plot_d, 
                                      ncol = 2, align = 'vh', labels = c('C', 'D')),
                   cowplot::plot_grid(pred_lcbd_plot_nd, 
                                      pred_lcbd_plot_d, 
                                      ncol = 2, align = 'vh', labels = c('E', 'F')), 
                   ncol = 1, rel_widths = c(1,1,1))

ggsave('figs/cs_2_plot_nd.png', width = 11, height = 17, units = 'cm', dpi = 600)
