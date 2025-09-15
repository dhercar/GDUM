# Author: Daniel Hern'andez Carrasco
# Email: dani.hc97@gmail.com
# Created: 8/13/2025
# License: MIT (see LICENSE file for details)

# DESCRIPTION: this script showcases the use of gdmmTMB for recovering directional (centroid) and non-directional (dispersion) environmetnal effects on beta diversity. 
set.seed(2)

# Load data
library(vegan)
library(ggplot2)
library(gdmmTMB)
library(MASS)

# plotting preferences
theme_set(theme_bw())
theme_update(strip.background = element_blank(), 
             text = element_text(size = 8),
             panel.grid = element_blank())

# Microbial data
data('microbialdata', package = 'gllvm') # requires gllvm package installed
env <- microbialdata$Xenv
env$site <- as.factor(1:nrow(env))
sp <- microbialdata$Y

### 1 NAIVE PATTERN ####
sp <- sp[order(env$pH),] # community data
env <- env[order(env$pH),] # environmental data
sp_dist <- as.matrix(vegdist(sp, binary = T)) # Sørensen dissimilarity matrix

# ORDINATION 1
kk_ord <- metaMDS(sp_dist, trymax = 200)
kk_ord$stress
points <- data.frame(kk_ord$points)
points$pH <- env$pH
points$LCBD <- adespatial::LCBD.comp(vegdist(sp, binary = T))$LCBD

(plot_NMDS <- ggplot(points, aes(x = MDS1, y = MDS2, fill = pH )) + 
    theme_void() +
    ylab('') +
    xlab('') +
    theme(legend.position = '',
          plot.margin = unit(c(1,1,1,1)*0.2, 'cm'),
          text = element_text(size = 8),
          axis.title.x = element_text(vjust = -2),
          axis.title.y = element_text(angle = 90, vjust = 0),
          panel.border = element_rect(fill = 'transparent')) +
    geom_hline(yintercept = 0, colour = 'grey') + 
    geom_vline(xintercept = 0, colour = 'grey') + 
    geom_point(shape = 21) + 
    coord_equal() +
    scale_fill_gradientn(colours = c('steelblue4', 'grey90', 'darkred')))


# LCBD vs. pH 

(plot_ushape <- ggplot(points, aes(y = LCBD, x = pH)) + 
  theme(panel.grid = element_blank(), aspect.ratio = 1, legend.position = '') +
  scale_fill_gradientn(colours = c('steelblue4', 'grey90', 'darkred')) +
  geom_smooth(colour = 'black', se = T, fill = 'grey90', alpha = 1, method = 'lm', formula = y~poly(x, 2)) + 
  geom_point(aes(fill = pH), shape = 21, size = 1.5))


### 2 DISSIMILARITY GRADIENT REMOVED ####

## 2 Residual pattern ##
m0 <- gdmm(Y = sp,
           X = env, 
           diss_formula = ~isp(pH, degree = 2), 
           uniq_formula = ~(1|site),#~ splines::bs(pH, df = 3),
           binary = T,
           family = 'normal',
           method = 'bray',
           link = 'identity',
           bboot = F,
           n_boot = 1000,
           mono = T)

pred <- predict(m0, new_X = env, new_W = env, CI = F)

dist_matrix <- matrix(0, nrow = nrow(env), ncol = nrow(env))
dist_matrix[cbind(m0$D[,1], m0$D[,2])] <- c(vegdist(sp, binary = T)) - pred
dist_matrix[cbind(m0$D[,2], m0$D[,1])] <- c(vegdist(sp, binary = T)) - pred
dist_matrix <- dist_matrix + abs(min(dist_matrix)) + 0.001
diag(dist_matrix) <- 0

##### RESIDUAL DISS. ####
kk_ord2 <- metaMDS(as.dist(dist_matrix), trymax = 200)
kk_ord2$stress
points2 <- data.frame(kk_ord2$points)
rot <- procrustes(points2[,c(1,2)], points[,c(1,2)], scale = F)
points2_rot <- data.frame(as.matrix(points2) %*% rot$rotation)
colnames(points2_rot) <- names(points2)
points2_rot$pH <- env$pH
points2_rot$LCBD <-  adespatial::LCBD.comp(as.dist(dist_matrix))$LCBD

(plot_NMDS2 <- ggplot(points2_rot, aes(x = MDS1, y = MDS2, fill = pH )) + 
    theme_void() +
    ylab('') +
    xlab('') +
    theme(legend.position = '',
          plot.margin = unit(c(1,1,1,1)*0.2, 'cm'),
          text = element_text(size = 8),
          axis.title.x = element_text(vjust = -2),
          axis.title.y = element_text(angle = 90, vjust = 0),
          panel.border = element_rect(fill = 'transparent')) +
    geom_hline(yintercept = 0, colour = 'grey') + 
    geom_vline(xintercept = 0, colour = 'grey') + 
    coord_equal() +
    geom_point(shape = 21) + 
    scale_fill_gradientn(colours = c('steelblue4', 'grey90','darkred')))


(plot_ushape2 <- ggplot(points2_rot, aes(y = LCBD, x = pH)) + 
  theme(panel.grid = element_blank(), aspect.ratio = 1, legend.position = '') +
  scale_fill_gradientn(colours = c('steelblue4','grey90','darkred')) +
    geom_smooth(colour = 'black', se = T, fill = 'grey90', alpha = 1, method = 'lm', formula = y~poly(x, 2)) + 
    geom_point(aes(fill = pH), shape = 21, size = 1.5))


cowplot::plot_grid(cowplot::plot_grid(plot_NMDS, plot_ushape, ncol = 1, rel_heights = c(2,3), align = 'v', labels = c('A', 'C')), 
                   cowplot:: plot_grid(plot_NMDS2, plot_ushape2, ncol = 1, rel_heights = c(2,3), align = 'v', labels = c('B', 'D')))


ggsave('figs/microbial_naive.eps', device = cairo_ps, width = 15, height = 8, units = 'cm', dpi = 1200)
ggsave('figs/microbial_naive.png', width = 15, height = 8, units = 'cm', dpi = 1200)


#### FULL MODEL ####
m <- gdmm(Y = sp,
          X = env, 
          diss_formula = ~ isp(pH, degree = 2),
          uniq_formula = ~ isp(pH, degree = 2) + (1|site),
          binary = T,
          family = 'binomial',
          method = 'bray',
          link = 'logit',
          mono = T)

summary(m)

pH <- env$pH
mean_pH <- rep(mean(env$pH), nrow(env))

# TOTAL effect pH (total)
pred_lcbd <- list()

newdata_W =  data.frame(pH = pH)
newdata_X =  data.frame(pH = pH)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     n_sim = 5000,
                     new_X = as.matrix(newdata_X),
                     new_W = as.matrix(newdata_W),
                     component = 'uniqueness',
                     scale_uniq = T,
                     type = 'response'),
             var = 'pH',
             value = pH,
             which = ' full')

# diss effect of pH 
newdata_W =  data.frame(pH = mean_pH)
newdata_X =  data.frame(pH = pH)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     n_sim = 5000,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = T,
                     re_sd = c('site'), # required to obtain site-level uncertainty
                     type = 'response'),
             var = 'pH',
             value = pH,
             which = 'dissimilarity gradient only')


# uniq effect pH
newdata_W =  data.frame(pH = pH)
newdata_X =  data.frame(pH = mean_pH )

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     n_sim = 5000,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = T,
                     re_sd = c('site'), # required to obtain site-level uncertainty
                     type = 'response'),
             var = 'pH',
             value = pH,
             which = 'direct effect only')


pred_lcbd <- do.call(rbind, pred_lcbd)

# FULL
(exp_LCBD_plot <- ggplot(subset(pred_lcbd, which == ' full'), aes(x = value, y = mean)) + 
    ylab('LCBD') +
    xlab('pH') +
    theme(strip.placement = 'outside',
          legend.position = '',panel.grid = element_blank(),
          aspect.ratio = 1) +
    geom_ribbon(aes(ymin = CI.2.5., ymax = CI.97.5.), fill = 'grey90') +
    geom_ribbon(aes(ymin = CI.25., ymax = CI.75.), fill = 'grey70') + 
    geom_line(col = 'black', linewidth = 1) +
    scale_fill_gradientn(colors =c('steelblue4', 'grey90', 'darkred')))

# Partial ###
diss_gradients <- do.call(rbind, diss_gradient(m0, CI_quant = c(0.95, 0.5)))

fx_plot <- ggplot(diss_gradients, aes(x = x, y = f_x)) + 
    ylab('f(pH)') +
    xlab('pH') +
    theme(strip.placement = 'outside',
          panel.grid = element_blank(),
          aspect.ratio = 0.75,
          plot.background = element_blank()) +
    geom_ribbon(aes(ymin = `CI 2.5%`, ymax = `CI 97.5%`),  fill = 'grey90') +
    geom_ribbon(aes(ymin = `CI 25%`, ymax = `CI 75%`),  fill = 'grey70') +
    geom_line(col = 'black', linewidth = 0.4)


(partial_LCBD_plot <- ggplot(subset(pred_lcbd, which != ' full')) + 
    ylab('LCBD') +
    xlab('pH') +
    theme(strip.placement = 'outside',
          strip.background = element_blank(),
          aspect.ratio = 1,
          panel.grid = element_blank(),
          legend.position = '') +
    geom_ribbon(aes(ymin = CI.2.5., ymax = CI.97.5., x = value),  fill = 'grey90') +
    geom_ribbon(aes(ymin = CI.25., ymax = CI.75.,x = value),  fill = 'grey70') +
    geom_line(aes(x = value, y = mean)) +
    facet_wrap(~which, ncol = 2))   # location on y-axis)
  
cowplot::plot_grid(cowplot::plot_grid(
  exp_LCBD_plot +
    geom_point(data = points, size = 2, aes(
      x = pH, y = LCBD, fill = pH), shape = 21),
  cowplot::plot_grid(
    fx_plot,
    partial_LCBD_plot,
    ncol = 1,
    rel_heights = c(1, 1.2),
    labels = c('B', 'C')),
  labels = 'A'
))

ggsave('figs/microbial_model_full.eps', device = cairo_ps, width = 18, height = 9, units = 'cm', dpi = 1200)
ggsave('figs/microbial_model_full.png', width = 18, height = 9, units = 'cm', dpi = 1200)



#### SIMULATE DIFFERENT GRADIENT ####
n <- nrow(env)
pH <- min(env$pH) + (runif(n)^0.5) * (max(env$pH)-min(env$pH))
mean_pH <- rep(mean(env$pH), nrow(env))

# TOTAL effect temmp2 (total)
pred_lcbd <- list()

newdata_W =  data.frame(pH = pH)
newdata_X =  data.frame(pH = pH)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     n_sim = 5000,
                     new_X = as.matrix(newdata_X),
                     new_W = as.matrix(newdata_W),
                     component = 'uniqueness',
                     scale_uniq = T,
                     type = 'response'),
             var = 'pH',
             value = pH,
             which = ' full')

pred_lcbd <- do.call(rbind, pred_lcbd)

# FULL
(exp_LCBD_plot_skewed <- ggplot(subset(pred_lcbd, which == ' full'), aes(x = value, y = mean)) + 
    geom_line(col = 'black', linewidth = 1) +
    ylab('LCBD') +
    xlab('pH') +
    theme(strip.placement = 'outside',
          legend.position = '',panel.grid = element_blank(),
          aspect.ratio = 1) +
    geom_ribbon(aes(ymin = CI.2.5., ymax = CI.97.5.), alpha = 0.1) +
    geom_ribbon(aes(ymin = CI.25., ymax = CI.75.), alpha = 0.3) + 
    geom_point(aes(y = 0.014, x = pH), shape = '|', col = 'darkred'))

cowplot::plot_grid(exp_LCBD_plot + ylim(0.0135, 0.033) +  geom_point(aes(y = 0.014, x = env$pH), shape = '|', col = 'darkred')
                   , exp_LCBD_plot_skewed + ylim(0.0135, 0.033))

ggsave('figs/scenarios.png', width = 15, height = 9, units = 'cm', dpi = 600)
ggsave('figs/scenarios.eps',device = cairo_ps, width = 15, height = 9, units = 'cm', dpi = 600)


#### dissimilaity gradient explained ####
m0 <- gdmm(Y = sp,
           X = env, 
           diss_formula = ~isp(pH, degree = 2), 
           binary = T,
           family = 'normal',
           method = 'bray',
           link = 'identity',
           bboot = T,
           n_boot = 1000,
           mono = T)

diss_gradients <- diss_gradient(m0, CI_quant = c(0.95, 0.5))$pH


fx_plot <- ggplot(diss_gradients, aes(x = x, y = f_x)) + 
  ylab('f(pH)') +
  xlab('pH') +
  theme(
    strip.placement = 'outside',
    panel.grid = element_blank(),
  ) +
  
  # Confidence ribbons and main line
  geom_ribbon(aes(ymin = `CI 2.5%`, ymax = `CI 97.5%`), fill = 'steelblue', alpha = 0.1) +
  geom_ribbon(aes(ymin = `CI 25%`, ymax = `CI 75%`), fill = 'steelblue', alpha = 0.2) +
  geom_line(col = 'steelblue4') + 
  
  # Vertical and horizontal segments for i
  geom_segment(aes(x = 6, y = 0, xend = 6, yend = 0.113)) + 
  geom_segment(aes(x = 8.95, y = 0.113, xend = 6, yend = 0.113)) + 
  
  # Vertical and horizontal segments for j
  geom_segment(aes(x = 7.5, y = 0, xend = 7.5, yend = 0.264)) + 
  geom_segment(aes(x = 7.5, y = 0.264, xend = 8.95, yend = 0.264)) + 
  
  # y-difference arrows (double-headed)
  geom_segment(aes(x = 8.85, y = 0.113, xend = 8.85, yend = 0.264),
               arrow = arrow(type = 'closed', length = unit(0.15, "cm"))) + 
  geom_segment(aes(x = 8.85, y = 0.264, xend = 8.85, yend = 0.113),
               arrow = arrow(type = 'closed', length = unit(0.15, "cm"))) +
  
  # x-difference arrows (double-headed), now closer to y = 0
  geom_segment(aes(x = 6, y = 0.01, xend = 7.5, yend = 0.01),
               arrow = arrow(type = 'closed', length = unit(0.15, "cm"))) +
  geom_segment(aes(x = 7.5, y = 0.01, xend = 6, yend = 0.01),
               arrow = arrow(type = 'closed', length = unit(0.15, "cm"))) +
  
  # y-difference annotation (vertical, rotated)
  annotate("text", x = 9.05, y = (0.113 + 0.264)/2,
           label = expression(h[ij] == group("|", f[pH](pH[i]) - f[pH](pH[j]), "|")),
           angle = 270, vjust = 0.5, size = 2.7) +
  
  # x-difference annotation, moved closer to x-axis
  annotate("text", x = (6 + 7.5)/2, y = 0.015,
           label = expression(group("|", pH[i] - pH[j], "|")),
           vjust = 0, size = 2.7)

ggsave( 'figs/example_fx_gradient.png', fx_plot, width = 12, height = 10, units = 'cm', dpi = 600)





