# Load packages 
library(greta)
library(vegan)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(latex2exp)

# Load functions
dir <- './Functions/'
files.sources = list.files(dir)
sapply(paste0(dir,files.sources), source)

# ggplot settings
theme_set(theme_bw())
theme_update(panel.grid = element_blank())
set.seed(1)

# Load datasets
doubs.spe <- read.csv('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/DoubsSpe.csv', row.names = 1)
doubs.env <- read.csv('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/DoubsEnv.csv', row.names = 1)
doubs.spa <- read.csv('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/DoubsSpa.csv', row.names = 1)
doubs.spa$s <- 1:30

# Create response matrix
sor_comp <- make_y_df( doubs.spe, method = 'abcd')
Y_sha <- 2*sor_comp$a # n species shared 
Y_rep <- pmin(sor_comp$b, sor_comp$c) # Replacement component (t)
Y_sim <- abs(sor_comp$b - sor_comp$c) # Richness difference component (r)
Y <- cbind(Y_rep, Y_sim, Y_sha) # Multivariate response

# Predictor (distance along the river)
X  <- scale(make_x_df(data.frame(doubs.env$das))[,3])

# Design
row <- sor_comp[, 1]
col <- sor_comp[, 2]

# ---- M1 FULL ----
# Note that although the response has three categories, the liear predictor contains only 2 as the third can be derived as (1-p(k_1) - p(k_1))
# priors
beta <- normal(0, 1, dim = c(1,2),  truncation = c(0, Inf)) 
alpha <- normal(0, 1, dim = 2)
SD_s <- normal(0, 1, truncation = c(0, Inf), dim = 2)
e_s  <- sweep(normal(0, 1, dim = c(nrow(doubs.spe),2)),2 , SD_s, FUN = '*') 
S_s <- e_s # Site level component only contains RE
eta <- sweep( X %*% beta + S_s[row,] + S_s[col,], 2, alpha, FUN = '+' ) 

# model
mu <- imultilogit(eta) # inverse link

# dimensions
trials <- rowSums(Y) 
n_realisations <- nrow(Y)

# distribution
distribution(Y) <- greta::multinomial(trials, mu, n_realisations)

# Estimate
m1 <- model(beta, alpha, e_s, SD_s)
draws <- greta::mcmc(m1, hmc(Lmin = 15, Lmax = 20), 
                     warmup = 5000,
                     n_samples = 10000)

 # check chains
gelman.diag(draws) #Rhat
bayesplot::mcmc_trace(draws, regex_pars = c('beta', 'alpha', 'SD'))
bayesplot::mcmc_rank_overlay(draws, regex_pars = c('beta', 'alpha', 'SD'))
bayesplot::mcmc_intervals(draws, regex_pars = c('beta', 'alpha', 'SD'))
bayesplot::mcmc_rank_overlay(draws, regex_pars = c('e_s'))
bayesplot::mcmc_intervals(draws, regex_pars = 'e_s')

# Calculate expected site re (proportional to u in this case)
e_s_m1 <- data.frame(calculate(e_s, values = draws, nsim = 1000))
e_s_m1 <- data.frame(t(apply(e_s_m1, 2, function(x) quantile(x, prob = c(0.025, 0.5, 0.975)))))
e_s_m1$k <- rep(c('rep', 'sim'), each = 30)
e_s_m1$s <- rep(1:30, 2)
e_s_m1$sig <- as.factor((sign(e_s_m1$X2.5.) + sign(e_s_m1$X97.5.)) / 2)
e_s_m1$col <- paste0(e_s_m1$k, e_s_m1$sig)

plot.data <- e_s_m1 %>% 
  inner_join(doubs.spa) %>%
  mutate(x = ifelse(k == 'rep', x - 170, x),
         y = ifelse(k == 'rep', y + 100, y))


(map.plot <- ggplot(plot.data, aes(x = x, y = y)) +  
  theme_void() +
  geom_path(colour = 'grey', size = 2, alpha = 0.5, aes(group = k)) + 
  coord_equal() +
  theme(legend.position = '',
        strip.background = element_blank()) + 
  geom_point(aes(size = X50., colour = sig, fill = col), shape = 21, stroke = 1) + 
  geom_text(data = subset(plot.data, s %in% c(8, 23:25)), aes(label = s, x = x + 15, y = y+10), size = 2.5) + 
  geom_point(data = subset(plot.data, sig != 0), aes(size = X50., colour = sig, fill = col), shape = 21, stroke = 1) + 
  scale_fill_manual('k',values = c( 'lightblue', 'lightblue3', 'steelblue3',
                                    'wheat', 'wheat3', 'darkorange')) + 
  scale_colour_manual('', values = c('black', 'transparent', 'black')) +
  geom_segment(aes(x = 120, y = 220, yend = 205, xend = 90),
               arrow = arrow(length = unit(0.1, 'inch')),
               colour = 'grey') + 
  geom_segment(aes(x = 120 - 170, y = 220 + 100, yend = 205 + 100, xend = 90 - 170),
               arrow = arrow(length = unit(0.1, 'inch')),
               colour = 'grey') +
  geom_segment(aes(x = 120, y = 20, yend = 50, xend = 130),
               arrow = arrow(length = unit(0.1, 'inch')),
               colour = 'grey') +
  geom_segment(aes(x = 120 - 170, y = 20 + 100, yend = 50 + 100, xend = 130 - 170),
               arrow = arrow(length = unit(0.1, 'inch')),
               colour = 'grey') +
  scale_size_continuous(range = c(0.05, 4)))

# Predictions for river-distance gradient
mu_1 <- imultilogit(
  sweep(
    seq(min(X), max(X), length.out = 100) %*% beta, 2, alpha, FUN = '+'))
pred_m1 <- data.frame(calculate(mu_1,
  values = draws, nsim = 1000))
pred_m1 <- data.frame(t(data.frame(apply(pred_m1,2 , function(x) quantile(x, prob = c(0.025, 0.25, 0.5, 0.75, 0.975))))))

pred_m1$X <- rep(seq(min(X), max(X), length.out = 100),3)
pred_m1$k <- rep(c('$P(k = t)$ \n $(Z_{t})$', '$P(k = r)$ \n $(Z_{r})$',' $P(k = s)$ \n$ (1- Z_{sor})$'), each = 100)


(diss.plot <- ggplot(pred_m1, aes(x = X, y = `X50.`, col = k, fill = k))  + 
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.2, col = 'transparent') +
  geom_ribbon(aes(ymin = X25., ymax = X75.), alpha = 0.5, col = 'transparent') +
  geom_line() + 
  xlab(TeX('$|x_i -x_j|$ (sd units)')) + 
  facet_wrap( ~ TeX(k, output = "character"), labeller = label_parsed, ncol = 1,
              strip.position = 'left') + 
  theme(legend.position = '',
        axis.title.y = element_blank(),
        strip.placement = 'outside',
        strip.background = element_blank()) + 
  scale_fill_manual('k',values = c( 'grey50',  'darkorange', 'steelblue4')) + 
  scale_colour_manual('k',values = c(  'grey20', 'darkorange3', 'darkblue'))) 


# Calculate u (x constant)
S_s_m1 <-  data.frame(calculate(sweep(e_s,2,alpha/2,FUN = '+' ),
                                values = draws, nsim = 1000))
S_s_m1 <- data.frame(data.frame(apply(S_s_m1,2 , function(x) quantile(x, prob = c(0.5)))))


# ---- M2 NO DIS ----
# priors
alpha <- normal(0, 1, dim = 2)
SD_s <- normal(0, 1, truncation = c(0, Inf), dim = 2)
e_s  <- sweep(normal(0, 1, dim = c(nrow(doubs.spe),2)),2 , SD_s, FUN = '*') 
S_s <- e_s 
eta <- sweep(S_s[row,] + S_s[col,], 2, alpha, FUN = '+' )

# model
mu <- imultilogit(eta)
trials <- rowSums(Y)
n_realisations <- nrow(Y)

distribution(Y) <- greta::multinomial(trials, mu, n_realisations)

m2 <- model(alpha, e_s, SD_s)
draws2 <- greta::mcmc(m2, hmc(Lmin = 15, Lmax = 20), 
                     warmup = 5000,
                     n_samples = 10000)

# check chains
gelman.diag(draws2) #Rhat
bayesplot::mcmc_trace(draws2, regex_pars = c('beta', 'alpha', 'SD'))
bayesplot::mcmc_rank_overlay(draws2, regex_pars = c('beta', 'alpha', 'SD'))
bayesplot::mcmc_intervals(draws2, regex_pars = c('beta', 'alpha', 'SD'))
bayesplot::mcmc_rank_overlay(draws2, regex_pars = c('e_s'))
bayesplot::mcmc_intervals(draws2, regex_pars = 'e_s')


# Calculate u 
S_s_m2 <-  data.frame(calculate(sweep(e_s,2,alpha/2,FUN = '+' ),
                                values = draws2, nsim = 1000))

S_s_m2 <- data.frame(data.frame(apply(S_s_m2,2 , function(x) quantile(x, prob = c(0.5)))))


# Change in u when accounting for environmental gradient
lcbd_comp <- data.frame(dif_m = (S_s_m1 - S_s_m2), k = rep(c('rep', 'sim'), each = 30), s = rep(1:30, 2))
names(lcbd_comp)[1] <- 'dif_m'

(change.plot <- ggplot(lcbd_comp, aes(x = as.factor(s), y = dif_m )) + 
  theme_classic() +
  geom_segment(aes(yend = 0, y = dif_m, colour = k), position = position_dodge(width = 0.7)) + 
  geom_hline(yintercept = 0) + 
  theme(legend.position = '',
        aspect.ratio = 0.2,
        axis.text.x = element_text(size = 6)) + 
  geom_point(aes(colour = k), position = position_dodge(width = 0.7)) + geom_hline(yintercept = 0) +
  scale_colour_manual('', values = c('steelblue4', 'darkorange')) + 
  scale_y_continuous(TeX("$\\Delta u_i$") , breaks = c(-1.5, -1, -0.5, 0, 0.5, 1)) + 
  xlab('site'))

# prop explained
plot_doubs <- plot_grid(diss.plot,
                        plot_grid(map.plot,change.plot, rel_heights = c(2,1), ncol = 1, labels = c(' ', 'C')),
                        rel_widths = c(1,2), 
                        labels = c('A', 'B'))



ggsave('plots/plot_doubs.png', plot_doubs, width = 18, height = 12, units = 'cm', dpi = 600)


summary(draws)[[2]]
summary(draws2)[[2]]
