# Load data
library(gllvm)
library(vegan)

data('microbialdata')
env <- microbialdata$Xenv
env$site <- as.factor(1:nrow(env))
sp <- microbialdata$Y


### 1 NAIVE PATTERN ####
sp <- sp[order(env$pH),]
env <- env[order(env$pH),]
sp_dist <- as.matrix(vegdist(sp, binary = T))


# ORDINATION 1
library(MASS)
kk_ord <- monoMDS(vegdist(sp, binary = T))
points <- data.frame(kk_ord$points)
points$pH <- env$pH
points$LCBD <- adespatial::LCBD.comp(vegdist(sp, binary = T))$LCBD

(plot_NMDS <- ggplot(points, aes(x = MDS1, y = MDS2, fill = pH )) + 
    theme_void() +
    ylab('MDS2') +
    xlab('MDS1') +
    theme(legend.position = '',
          plot.margin = unit(c(1,1,1,1)*0.2, 'cm'),
          text = element_text(size = 8),
          axis.title.x = element_text(vjust = -2),
          axis.title.y = element_text(angle = 90, vjust = 0),
          panel.border = element_rect(fill = 'transparent')) +
    geom_hline(yintercept = 0, colour = 'grey') + 
    geom_vline(xintercept = 0, colour = 'grey') + 
    geom_point(shape = 21, aes(size  = LCBD)) + 
    coord_equal() +
    scale_size_continuous(range = c(0.5,2)) +
    scale_fill_gradientn(colours = c('steelblue4', 'grey90', 'coral2')))


# LCBD vs. pH 

plot_ushape <- ggplot(points, aes(y = LCBD, x = pH)) + 
  theme(panel.grid = element_blank(), aspect.ratio = 1, legend.position = '') +
  geom_smooth(colour = 'black', se = T, method = 'lm', formula = y~poly(x, 2)) + 
  scale_fill_gradientn(colours = c('steelblue4', 'grey90', 'coral2')) +
  geom_point(aes(fill = pH), shape = 21, size = 2)


### 2 DISSIMILARITY GRADIENT REMOVED ####

## 2 Residual pattern ##
m0 <- gdmm(Y = sp,
           X = env, 
           diss_formula = ~isp(pH, degree = 2), # In the future: s( bdo ) + s(nit)
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

kk_ord2 <- monoMDS(as.dist(dist_matrix))
points2 <- data.frame(kk_ord2$points)
rot <- procrustes(points2[,c(1,2)], points[,c(1,2)], scale = F)
points2_rot <- data.frame(as.matrix(points2) %*% rot$rotation)
colnames(points2_rot) <- names(points2)
points2_rot$pH <- env$pH

points2_rot$LCBD <-  adespatial::LCBD.comp(as.dist(dist_matrix))$LCBD

(plot_NMDS2 <- ggplot(points2_rot, aes(x = MDS1, y = MDS2, fill = pH )) + 
    theme_void() +
    ylab('MDS2') +
    xlab('MDS1') +
    theme(legend.position = '',
          plot.margin = unit(c(1,1,1,1)*0.2, 'cm'),
          text = element_text(size = 8),
          axis.title.x = element_text(vjust = -2),
          axis.title.y = element_text(angle = 90, vjust = 0),
          panel.border = element_rect(fill = 'transparent')) +
    geom_hline(yintercept = 0, colour = 'grey') + 
    geom_vline(xintercept = 0, colour = 'grey') + 
    coord_equal() +
    geom_point(shape = 21, aes(size  = LCBD)) + 
    scale_size_continuous(range = c(0.5,2)) +
    scale_fill_gradientn(colours = c('steelblue4', 'grey90','coral2')))


plot_ushape2 <- ggplot(points2_rot, aes(y = LCBD, x = pH)) + 
  theme(panel.grid = element_blank(), aspect.ratio = 1, legend.position = '') +
  geom_smooth(colour = 'black', se = T, method = 'lm', formula = y~poly(x, 2)) + 
  scale_fill_gradientn(colours = c('steelblue4','grey90','coral2')) +
  geom_point(aes(fill = pH), shape = 21, size = 2)


plot_grid(plot_grid(plot_NMDS, plot_ushape, ncol = 1, rel_heights = c(2,3), align = 'v', labels = c('A', 'B')), 
          plot_grid(plot_NMDS2, plot_ushape2, ncol = 1, rel_heights = c(2,3), align = 'v', labels = c('B', 'C')))

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

# TOTAL effect temmp2 (total)
pred_lcbd <- list()

newdata_W =  data.frame(pH = pH)
newdata_X =  data.frame(pH = pH)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     re_sd = c('site'),
                     scale_uniq = T,
                     type = 'response'),
             var = 'pH',
             value = pH,
             which = ' full')

# diss effect of temp (total)
newdata_W =  data.frame(pH = mean_pH)
newdata_X =  data.frame(pH = pH)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = T,
                     re_sd = c('site'),
                     type = 'response'),
             var = 'pH',
             value = pH,
             which = 'diss. gradient only')


# uniq effect temp (total)
newdata_W =  data.frame(pH = pH)
newdata_X =  data.frame(pH = mean_pH )

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = T,
                     re_sd = c('site'),
                     type = 'response'),
             var = 'pH',
             value = pH,
             which = 'direct only')


pred_lcbd <- do.call(rbind, pred_lcbd)

# FULL
(exp_LCBD_plot <- ggplot(subset(pred_lcbd, which == ' full'), aes(x = value, y = mean)) + 
    geom_point(data = points, aes(x = pH, y = LCBD), shape = 21, fill = 'white') +
    geom_line(col = 'black') +
    ylab('E(LCBD)') +
    xlab('pH') +
    theme(strip.placement = 'outside',
          legend.position = '',panel.grid = element_blank(),
          aspect.ratio = 1) +
    geom_ribbon(aes(ymin = CI.2.5., ymax = CI.97.5.), fill = 'transparent', lty = 5, col = 'black'))

# Partial ###

(partial_LCBD_plot <- ggplot(subset(pred_lcbd, which != ' full'), aes(x = value, y = mean)) + 
    geom_point(aes(col = which), shape = 21) +
    geom_line(aes(col = which)) +
    ylab('E(LCBD)') +
    xlab('pH') +
    theme(strip.placement = 'outside',
          strip.background = element_rect(colour = 'grey20', fill = c('transparent')),
          aspect.ratio = 1,
          panel.grid = element_blank(),
          legend.position = '') +
    geom_ribbon(aes(ymin = CI.2.5., ymax = CI.97.5., col = which), fill = 'transparent', lty = 5) +
    scale_colour_manual(values = c('orange2', 'steelblue')) + 
    facet_wrap(~which, scales = 'free_y', ncol = 1))

diss_gradients <- do.call(rbind, diss_gradient(m0, CI_quant = c(0.95, 0.5)))

(fx_plot <- ggplot(diss_gradients, aes(x = x, y = f_x)) + 
    ylab('f(pH)') +
    xlab('pH') +
    theme(strip.placement = 'outside',
          panel.grid = element_blank(),
          aspect.ratio = 1) +
    geom_ribbon(aes(ymin = `CI 2.5%`, ymax = `CI 97.5%`), fill = 'steelblue', alpha = 0.1) +
    geom_ribbon(aes(ymin = `CI 25%`, ymax = `CI 75%`), fill = 'steelblue', alpha = 0.2) +
    geom_line(col = 'steelblue4') )


plot_grid(plot_grid(plot_NMDS, plot_ushape, ncol = 1, rel_heights = c(2,4), align = 'v', labels = c('A', 'B')), 
          plot_grid(plot_NMDS2, plot_ushape2, ncol = 1, rel_heights = c(2,4), align = 'v', labels = c('C', 'D')))

ggsave('figs/microbial_naivee.pdf', width = 15, height = 8, units = 'cm', dpi = 600)
ggsave('figs/microbial_naivee.png', width = 15, height = 8, units = 'cm', dpi = 600)

LCBD_pred <- data.frame(predict(m, component = 'uniqueness'))
LCBD_pred$real_LCBD <- points$LCBD

(plot_1to1 <- ggplot(LCBD_pred, aes(x = real_LCBD, y = mean)) + 
    geom_abline(slope = 1, intercept = 0) + 
    ylab('E(LCBD)') + 
    xlab('LCBD') +
    theme(aspect.ratio = 1)+
    geom_linerange(aes(ymax = `CI.97.5.`, ymin = `CI.2.5.`), colour = 'grey') +
    geom_point(size = 1, shape = 21, fill = 'white'))

plot_grid(plot_grid(exp_LCBD_plot, fx_plot, rel_heights = c(1,1), ncol = 1, align = 'v', labels = c('A', 'B')), partial_LCBD_plot, plot_1to1, ncol  = 3, labels = c('','C', 'D'))
ggsave('figs/microbial_model.pdf', width = 18, height = 10, units = 'cm', dpi = 600)
ggsave('figs/microbial_model.svg', width = 18, height = 10, units = 'cm', dpi = 600)

