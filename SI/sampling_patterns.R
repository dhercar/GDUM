# LOAD FUNCTTIONS
source('./functions/convenience_functions.R')
source('./functions/sim_com_data.R')

library(ggplot2)
library(ggforce)
library(tidyverse)
library(adespatial)
library(glmmTMB)
library(cowplot)
library(rvmethod)
library(vegan)

set.seed(5)

# parameters
n_sp = 100
n_samples = 100
env_range = c(-0.5,1.5)
tol_range = c(0.1, 0.4)



# ---- 1: Distribution of speicies optima / tolerance----
opt_sp <- cbind(opt_x_1 = runif(n_sp, min(env_range), max(env_range)),
                opt_x_2 = runif(n_sp, min(env_range), max(env_range)))

tol_sp <- cbind(tol_x_1 = exp(runif(n_sp,min(log(tol_range)), max(log(tol_range)))),
                tol_x_2 =  exp(runif(n_sp,min(log(tol_range)), max(log(tol_range)))))


# ggplot(cbind(opt_sp, tol_sp), aes(x = opt_x_1, y = opt_x_2, col = as.factor(1:n_sp))) +
#   theme(legend.position = '') +
#   coord_equal() +
#   geom_ellipse(aes(x0 = opt_x_1 , y0 = opt_x_2, a = 2*tol_x_1, b = 2*tol_x_2, angle = 0))
# 

# ---- 2: Sampling patterns ----

# Circular sampling
radius = 0.45
angles <- seq(0, 2 * pi, length.out = n_samples + 1)[-n_sp-1]

samples_circle <- data.frame(x_1 = 0.5 + radius * cos(angles) + rnorm(n_sp, 0, 0.03),
                            x_2 =  0.5 + radius * sin(angles) + rnorm(n_sp, 0, 0.03),
                            type = 'circle')
samples_circle[,1] <- ifelse(samples_circle[,1] >1 , 1, samples_circle[,1])
samples_circle[,1] <- ifelse(samples_circle[,1] <0 , 0, samples_circle[,1])
samples_circle[,2] <- ifelse(samples_circle[,2] >1 , 1, samples_circle[,2])
samples_circle[,2] <- ifelse(samples_circle[,2] <0 , 0, samples_circle[,2])

plot(samples_circle[,1:2], ylim = c(0,1), xlim = c(0,1))

# Normally distributed samples
samples_norm <-  data.frame(x_1 = rnorm(n_samples, 0.5, 0.2),
                            x_2 = rnorm(n_samples, 0.5, 0.2),
                            type = 'normal')
samples_norm[,1] <- ifelse(samples_norm[,1] >1 , 1, samples_norm[,1])
samples_norm[,1] <- ifelse(samples_norm[,1] <0 , 0, samples_norm[,1])
samples_norm[,2] <- ifelse(samples_norm[,2] >1 , 1, samples_norm[,2])
samples_norm[,2] <- ifelse(samples_norm[,2] <0 , 0, samples_norm[,2])

plot(samples_norm[,1:2], ylim = c(0,1), xlim = c(0,1))

# skewed distribution of samples
samples_skewed <- data.frame(x_1 = 1-(runif(n_samples, (0.0001)^(1/3), (1)^(1/3)))^3,
                             x_2 =  (runif(n_samples, (0.0001)^(1/3), (1)^(1/3)))^3,
                             type ='skewed')

plot(samples_skewed[,1:2],ylim = c(0,1), xlim = c(0,1))

# Grid sampling

samples_grid <- data.frame(x_1 = rep(seq(0.1, 0.9, length.out = sqrt(n_samples)), each = sqrt(n_samples)),
                           x_2 = rep(seq(0.1, 0.9, length.out = sqrt(n_samples)), sqrt(n_samples)), 
                           type ='grid')

# merge
samples <- rbind(samples_circle, samples_norm, samples_skewed, samples_grid)
rm(samples_circle, samples_skewed, samples_norm)

# ---- 3: Species composition ----
# Calculate probability of occurrence from optima tolerance
out.list <- list()
for (i in unique(samples$type)) {
  print(i)
  kk <- subset(samples, type == i)
  for (j in 1:nrow(kk)) {
    for (sp in 1:n_sp) {
      p_x_1 =  gaussfunc(kk[j,1], 
                           mu = opt_sp[sp,1],
                           sigma = tol_sp[sp,1])
      p_x_2 =  gaussfunc(kk[j,2], 
                           mu = opt_sp[sp,2],
                           sigma = tol_sp[sp,2])
      
      p_sp = p_x_1*p_x_2
      
      out.list[[length(out.list) + 1]] <- data.frame(type = i,
                                                   x_1 = kk[j,1],
                                                   x_2 = kk[j,2],
                                                   sp = sp,
                                                   p = p_sp)
    }
  }
}

out.data <- data.frame(do.call(rbind, out.list))

# Simulate draws from binomial distribution to obtain PA
out.data$PA <- vrbinom(1,1, out.data$p)
out.data <- out.data %>% pivot_wider(names_from = sp, values_from = PA, id_cols = c(x_1, x_2, type))

# ---- 4: Calculate LCBD ----

data_LCBD = NULL
for (i in unique(out.data$type)) {
  kk <- subset(out.data, type == i)
  LCBD <- LCBD.comp(vegdist(kk[,-c(1:3)]))$LCBD
  data_LCBD = rbind(data_LCBD, data.frame(type = i, x_1 = kk$x_1, x_2 = kk$x_2, LCBD = LCBD))
}

library(splines)
m1 <- glmmTMB(LCBD~ ns(x_1, 2) + ns(x_2,2), data = subset(data_LCBD,type == 'circle'))
m2 <- glmmTMB(LCBD~ ns(x_1, 2) + ns(x_2,2), data = subset(data_LCBD,type == 'normal'))
m3 <- glmmTMB(LCBD~ ns(x_1, 2) + ns(x_2,2), data = subset(data_LCBD,type == 'skewed'))
m4 <- glmmTMB(LCBD~ ns(x_1, 2) + ns(x_2,2), data = subset(data_LCBD,type == 'grid'))

pred.data <- rbind(
  cbind(data.frame(ggeffects::ggpredict(m1, terms = c('x_1[0:1, by=0.01]'), type = 'fixed')), var = 'x_1', type = 'circle'),
  cbind(data.frame(ggeffects::ggpredict(m1, terms = c('x_2[0:1, by=0.01]'))), var = 'x_2', type = 'circle'),
  cbind(data.frame(ggeffects::ggpredict(m2, terms = c('x_1[0:1, by=0.01]'))), var = 'x_1', type = 'normal'),
  cbind(data.frame(ggeffects::ggpredict(m2, terms = c('x_2[0:1, by=0.01]'))), var = 'x_2', type = 'normal'),
  cbind(data.frame(ggeffects::ggpredict(m3, terms = c('x_1[0:1, by=0.01]'))), var = 'x_1', type = 'skewed'),
  cbind(data.frame(ggeffects::ggpredict(m3, terms = c('x_2[0:1, by=0.01]'))), var = 'x_2', type = 'skewed'),
  cbind(data.frame(ggeffects::ggpredict(m4, terms = c('x_1[0:1, by=0.01]'))), var = 'x_1', type = 'grid'),
  cbind(data.frame(ggeffects::ggpredict(m4, terms = c('x_2[0:1, by=0.01]'))), var = 'x_2', type = 'grid'))


plot_grid(
ggplot(samples, aes( x = x_1, y= x_2)) +
  theme_bw()+
  geom_point(aes(col = type, fill=type), shape =21, alpha = 0.7) +
  scale_colour_manual(values = c('steelblue', 'darkorange', 'darkolivegreen', 'brown')) +
  scale_fill_manual(values = c('steelblue', 'darkorange', 'darkolivegreen', 'brown')) +
  facet_grid(type ~ 1) +
  theme(aspect.ratio = 1, 
        plot.title = element_text(hjust = 0.5), 
        strip.text = element_blank(),
        panel.grid = element_blank(),
        legend.position = ''),
ggplot(pred.data, aes(x = x, y = predicted))+
  geom_line(aes(col = type))+
  theme_bw()+
  theme(panel.border = element_rect(fill= 'transparent'), 
        strip.background = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1, 
        plot.title = element_text(hjust = 0.5), 
        legend.position = '') +
   geom_point(data = pivot_longer(samples, cols = c(x_1, x_2), names_to = 'var'), 
              aes(x = value, y = 0.007, col = type), shape = 3 )+
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high, fill = type), alpha = 0.2) +
  scale_x_continuous('environment', expand = c(0,0), limits = c(0,1), breaks = c(0, 0.2,0.4, 0.6, 0.8)) +
  scale_y_continuous('LCBD',expand = c(0,0), limits = c(0.007, 0.017)) + 
  scale_colour_manual(values = c('steelblue', 'darkorange', 'darkolivegreen', 'brown')) +
  scale_fill_manual(values = c('steelblue', 'darkorange', 'darkolivegreen', 'brown')) +
  facet_grid(type~var), align = 'hv', rel_widths = c(1,2), labels = c('A', 'B'))


ggsave('SI/sampling_patterns.svg', width = 16, height = 18, dpi = 600, units = 'cm')
