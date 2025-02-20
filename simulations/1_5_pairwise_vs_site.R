library(greta)
library(latex2exp)
library(ggplot2)
library(extrafont)
library(sysfonts)
library(showtext)
library(cowplot)

# ggplot settings
theme_set(theme_bw())
theme_update(panel.grid = element_blank())
set.seed(1)

## This allows us to annotate ggplots in specific facets
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}

n_pred = 100

# Load functions and models
dir <- './functions/'
files.sources = list.files(dir)
sapply(paste0(dir,files.sources), source)

m2_2 <- readRDS('models/m2_2.rds')
m4 <- readRDS('models/m4.rds')
m5 <- readRDS('models/m5.rds')
m6 <- readRDS('models/m6.rds')
m7 <- readRDS('models/m7.rds')

# Calculate expected values
pred_u_m4_exp <- with(m4, {
  pred_data <- data.frame(X_iso = mean(X_iso),
                          X_e_uni = seq(min(X_e_uni), max(X_e_uni), length.out = n_pred))
  u <- alpha/2 + pred_data$X_e_uni %*% beta_e_uni
  SS <- ilogit(u*2)
  pred <- data.frame(calculate(SS, values = draws, nsim = 10000))
  pred_sum <- t(apply(pred, 2, function(x)
    quantile(as.numeric(x), probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
  data.frame(cbind(pred_sum, pred_data, m = 'm4'))
})

pred_u_m5_exp <- with(m5, {
  pred_data <- data.frame(X_iso = mean(X_iso),
                          X_e_uni = seq(min(X_e_uni), max(X_e_uni), length.out = n_pred))
  u <- alpha/2 + pred_data$X_e_uni %*% beta_e_uni
  SS <- ilogit(u*2)
  pred <- data.frame(calculate(SS, values = draws, nsim = 10000))
  pred_sum <- t(apply(pred, 2, function(x)
    quantile(as.numeric(x), probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
  data.frame(cbind(pred_sum, pred_data, m = 'm5'))
})

pred_diss_m5_exp <- with(m5, {
  pred_data <- data.frame(X_iso = mean(X_iso),
                          X_e_dis = seq(min(X_e_dis), max(X_e_dis), length.out = n_pred))
  eta <- alpha + pred_data$X_e_dis %*% beta_e_dis
  mu <- ilogit(eta)
  pred <- data.frame(calculate(mu, values = draws, nsim = 10000))
  pred_sum <- t(apply(pred, 2, function(x)
    quantile(as.numeric(x), probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
  data.frame(cbind(pred_sum, pred_data, m = 'm5'))
})

pred_u_m6_exp <- with(m6, {
  pred_data <- data.frame(X_iso = mean(X_iso),
                          X_env = seq(min(scale(env_raw)), max(scale(env_raw)), length.out = n_pred))

  X_env2 <- predict(X_env, pred_data$X_env)
  u <- alpha/2 + X_env2 %*% beta_env
  SS <- ilogit(u*2)
  pred <- data.frame(calculate(SS, values = draws, nsim = 10000))
  pred_sum <- t(apply(pred, 2, function(x)
    quantile(as.numeric(x), probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
  data.frame(cbind(pred_sum, pred_data, m = 'm6'))
})

pred_u_m7_exp <- with(m7, {
  
  pred_data <- data.frame(X_iso = mean(X_iso),
                          X_env = seq(min(scale(env_raw)), max(scale(env_raw)), length.out = n_pred))
  
  X_env2 <- predict(X_env, pred_data$X_env)
  u <- alpha/2 + X_env2 %*% beta_env
  SS <- ilogit(u*2)
  pred <- data.frame(calculate(SS, values = draws, nsim = 10000))
  pred_sum <- t(apply(pred, 2, function(x)
    quantile(as.numeric(x), probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
  data.frame(cbind(pred_sum, pred_data, m = 'm7'))
})

pred_diss_m7_exp <- with(m7, {
  pred_data <- data.frame(X_iso = mean(X_iso),
                          X_e_dis = seq(min(X_e_dis), max(X_e_dis), length.out = n_pred))
  eta <- alpha + pred_data$X_e_dis %*% beta_e_dis 
  mu <- ilogit(eta)
  pred <- data.frame(calculate(mu, values = draws, nsim = 10000))
  pred_sum <- t(apply(pred, 2, function(x)
    quantile(as.numeric(x), probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
  data.frame(cbind(pred_sum, pred_data, m = 'm7'))
})

# Plots  
diss_m5_plot <- ggplotGrob(ggplot(pred_diss_m5_exp, aes(x = X_e_dis, y = `X50.`)) +
    theme(legend.position ='', 
          panel.background = element_blank(), 
          axis.ticks = element_line(linewidth = 0.2),
          plot.background = element_blank(),
          text = element_text(size = 6)) +
    geom_ribbon(aes(ymin = `X2.5.`, ymax = `X97.5.`, 
                    fill = '95% CI'), alpha = 0.5) + 
    geom_ribbon(aes(ymin = `X25.`, ymax = `X75.`, 
                    fill = '50% CI'), alpha = 0.5) + 
    geom_line(colour = 'darkblue', size = 0.5) + 
    scale_y_continuous(breaks = c(0.4, 0.6, 0.8))+
    scale_fill_manual('', values = c('steelblue','lightblue3')) + 
    ylab(expression(bold(E)*"(Z"["ij"]*")") ) +
    xlab( TeX("x") ))

(uniq_plot <- ggplot(rbind(pred_u_m4_exp, pred_u_m5_exp)) + 
    geom_ribbon(aes(ymin = `X2.5.`, ymax = `X97.5.`, x = X_e_uni, y = `X50.`,
                    fill = '95% CI'), alpha = 0.2) + 
    geom_ribbon(aes(ymin = `X25.`, ymax = `X75.`, x = X_e_uni, y = `X50.`,
                    fill = '50% CI'), alpha = 0.5) + 
    geom_line(colour = 'coral4', size = 1, aes(x = X_e_uni, y = `X50.`)) + 
    scale_fill_manual('', values = c('coral', 'coral')) +
    ylab(TeX("uniqueness ( $\\u_{i}$)")) +
    xlab(TeX('u(x) (sd units)')) +
    facet_wrap(~m) + 
    ylim(0.4,0.8)+
    theme(legend.position = '', strip.background = element_blank(),text = element_text(size = 7)) +
    annotation_custom2(grob = diss_m5_plot, data = data.frame(m="m5"), ymin = 0.57, ymax = 0.82, xmin = 0, xmax = 2.3))

diss_m7_plot <- ggplotGrob(ggplot(pred_diss_m7_exp, aes(x = X_e_dis, y = `X50.`)) +
                             theme(legend.position ='',
                                   panel.background = element_blank(),
                                   text = element_text(size = 6), 
                                   plot.background = element_blank(),
                                   axis.ticks = element_line(linewidth = 0.2)) +
                             geom_ribbon(aes(ymin = `X2.5.`, ymax = `X97.5.`, 
                                             fill = '95% CI'), alpha = 0.5) + 
                             geom_ribbon(aes(ymin = `X25.`, ymax = `X75.`, 
                                             fill = '50% CI'), alpha = 0.5) + 
                             geom_line(colour = 'darkblue', size = 0.5) + 
                             
                             scale_y_continuous(breaks = c(0.4, 0.6, 0.8))+
                             scale_fill_manual('', values = c('steelblue','lightblue3')) + 
                             ylab(expression(bold(E)*"(Z"["ij"]*")") ) +
                             xlab( TeX("x") ))

(env_plot <- ggplot(rbind(pred_u_m6_exp, pred_u_m7_exp)) + 
    geom_ribbon(aes(ymin = `X2.5.`, ymax = `X97.5.`, x = X_env, y = `X50.`,
                    fill = '95% CI'), alpha = 0.5) + 
    geom_ribbon(aes(ymin = `X25.`, ymax = `X75.`, x = X_env, y = `X50.`,
                    fill = '50% CI'), alpha = 0.5) + 
    geom_line(colour = 'darkorange3', size = 1, aes(x = X_env, y = `X50.`)) + 
    scale_fill_manual('', values = c('orange', 'wheat')) +
    ylab(TeX("uniqueness ( $\\u_{i}$)")) +
    xlab(TeX('x (sd units)')) +
    facet_wrap(~m) + 
    ylim(0.4,0.8)+
    theme(legend.position = '', strip.background = element_blank(),
          text = element_text(size = 7)) +
    annotation_custom2(grob = diss_m7_plot, data = data.frame(m="m7"), ymin = 0.57, ymax = 0.82, xmin = -0.4, xmax = 1.75))


## Summary of coeff
draws4_df <- data.frame(do.call(rbind, m4$draws))
draws4_df  <- draws4_df[,1:2]
draws4_df$m <- 'm4'

draws5_df <- data.frame(do.call(rbind, m5$draws))
draws5_df  <- draws5_df[,1:3]
draws5_df$m <- 'm5'

draws6_df <- data.frame(do.call(rbind, m6$draws))
draws6_df  <- draws6_df[,1:3]
draws6_df$m <- 'm6'

draws7_df <- data.frame(do.call(rbind, m7$draws))
draws7_df  <- draws7_df[,1:4]
draws7_df$m <- 'm7'

library(tidyverse)
library(ggridges)
draws_data <- bind_rows(draws4_df, draws5_df, draws6_df, draws7_df) %>% 
  pivot_longer(-c(m)) %>% na.omit()


(coef_plot_m6_7 <- ggplot(subset(draws_data, m %in% c('m6', 'm7')), aes(x = value, y = name )) + 
  facet_wrap(~m, ncol =2) +
  geom_vline(xintercept = 0, colour = 'grey50', size = 0.3, lty = 5) + 
  stat_summary(geom = 'linerange', 
               fun.min = function(x) quantile(na.omit(x), probs = 0.025),
               fun.max = function(x) quantile(na.omit(x), probs = 0.975),
               alpha = 0.5,
               aes(col = name),
               position = position_dodge(0.5)) +
  stat_summary(geom = 'point', 
                 fun = mean,
                 shape = 21,
                 fill = 'white',
                 size = 1.5,
                 position = position_dodge(0.5)) +
  stat_summary(geom = 'point', 
               fun = mean,
               shape = 21,
               aes(fill = name, colour = name), 
               size = 1.5,
               alpha = 0.5,
               position = position_dodge(0.5)) +
  theme_classic() +
    ylab('model coefficient') +
    xlab('effect') +
    xlim(NA, 0.7) + 
  theme(strip.placement = 'outside',
        panel.spacing = unit(0.5, "cm", data = NULL),
        strip.background = element_blank(),
        strip.text = element_blank(),
        text = element_text(size = 7),
        axis.line.y = element_blank(),
        legend.position = '',
        legend.position.inside = c(0.5,0.6)) + 
    coord_flip() +
    scale_y_discrete(labels =  c(
      TeX('$\\beta_{\\Delta x}$', italic = TRUE),
      TeX('$\\lambda_{x_1}$', italic = TRUE),
      TeX('$\\lambda_{x_2}$', italic = TRUE),
      TeX('$\\lambda_{w}$', italic = TRUE)
    )) +
  scale_colour_manual('', values = c('steelblue', 'darkorange', 'darkorange', 'grey20'), 
                      labels = c(TeX('$\\beta_{\\Delta x}$', italic = TRUE), 
                                 TeX('$\\lambda_{x_1}$', italic = TRUE),
                                 TeX('$\\lambda_{x_2}$', italic = TRUE),
                                 TeX('$\\lambda_{w}$', italic = TRUE))) +
  scale_fill_manual('', values = c('steelblue', 'darkorange', 'darkorange', 'grey20'), 
                      labels = c(TeX('$\\beta_{\\Delta x}$', italic = TRUE), 
                                 TeX('$\\lambda_{x_1}$', italic = TRUE),
                                 TeX('$\\lambda_{x_2}$', italic = TRUE),
                                 TeX('$\\lambda_{w}$', italic = TRUE))))

(coef_plot_m4_5 <- ggplot(subset(draws_data, m %in% c('m4', 'm5')), aes(x = value, y = name )) + 
    facet_wrap(~m, ncol =2) +
    geom_vline(xintercept = 0, colour = 'grey50', size = 0.3, lty = 5) + 
    stat_summary(geom = 'linerange', 
                 fun.min = function(x) quantile(na.omit(x), probs = 0.025),
                 fun.max = function(x) quantile(na.omit(x), probs = 0.975),
                 alpha = 0.5,
                 aes(col = name),
                 position = position_dodge(0.5)) +
    stat_summary(geom = 'point', 
                 fun = mean,
                 shape = 21,
                 fill = 'white',
                 size = 1.5,
                 position = position_dodge(0.5)) +
    stat_summary(geom = 'point', 
                 fun = mean,
                 shape = 21,
                 aes(fill = name, colour = name), 
                 size = 1.5,
                 alpha = 0.5,
                 position = position_dodge(0.5)) +
    theme_classic() +
    theme(strip.placement = 'outside',
          panel.spacing = unit(0.5, "cm", data = NULL),
          strip.background = element_blank(),
          strip.text = element_blank(),
          text = element_text(size = 7),
          axis.line.y = element_blank(),
          legend.position = '',
          legend.position.inside = c(0.5,0.6)) + 
    ylab('model coefficient') +
    xlab('effect') +
    xlim(NA, 0.7) + 
    coord_flip() +
    scale_y_discrete(labels = c(TeX('$\\beta_{\\Delta x}$', italic = TRUE), 
                                TeX('$\\lambda_{u}$', italic = TRUE),
                                TeX('$\\lambda_{w}$', italic = TRUE))) +
    scale_colour_manual('', values = c('steelblue','coral3','grey20'), 
                        labels = c(TeX('$\\beta_{\\Delta x}$', italic = TRUE), 
                                   TeX('$\\lambda_{u}$', italic = TRUE),
                                   TeX('$\\lambda_{w}$', italic = TRUE))) +
    scale_fill_manual('', values = c('steelblue','coral3','grey20'), 
                      labels = c(TeX('$\\beta_{\\Delta x}$', italic = TRUE), 
                                 TeX('$\\lambda_{u}$', italic = TRUE),
                                 TeX('$\\lambda_{w}$', italic = TRUE)))) 
  


pair_vs_site_fig <- plot_grid(plot_grid(uniq_plot, coef_plot_m4_5, rel_heights = c(2,1), ncol =1, align = 'v', labels = c('A', 'C')), 
          plot_grid(env_plot, coef_plot_m6_7 , rel_heights = c(2,1), ncol =1, align = 'v', labels = c('B', 'D')),
          ncol = 2)

ggsave( 'plots/pair_vs_site_plot.png',pair_vs_site_fig, width = 18, height = 8, dpi = 600, units = 'cm')


