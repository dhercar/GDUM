library(ggplot2)
theme_set(theme_bw())
theme_update(strip.background = element_blank(),
             panel.grid = element_blank())
# Load data
data('doubs')
env <- doubs$env
com <- doubs$fish
env$site <- rownames(env)

m <- gdmm(Y = com, 
          X = env, 
          diss_formula = ~ isp(dfs),
          uniq_formula = ~ (1|site),
          method = 'bray', 
          binary = TRUE, 
          family = 'binomial', 
          link = 'logit',
          mono = T,
          bboot = F) 

summary(m)

n = 30
nit = sort(env$nit)
bdo = sort(env$bdo)
mean_nit <- rep(mean(env$nit), n)
mean_bdo <- rep(mean(env$bdo), n)
pred_lcbd <- list()

# TOTAL effect of nitrogen (total)
newdata_W =  data.frame(bdo = mean_bdo, nit = nit)
newdata_X =  data.frame(bdo = mean_bdo, nit = nit)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
          new_X = newdata_X,
          new_W = newdata_W,
          component = 'uniqueness',
          scale_uniq = F,
          type = 'response'),
        var = 'nit',
        value = nit,
        which = ' full')

# TOTAL effect of bdo (total)
newdata_W =  data.frame(bdo = bdo, nit = mean_nit)
newdata_X =  data.frame(bdo = bdo, nit = mean_nit)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                new_X = newdata_X,
                new_W = newdata_W,
                component = 'uniqueness',
                scale_uniq = F,
                type = 'response'),
        var = 'bdo',
        value = bdo,
        which = ' full')

# diss effect of nit (total)
newdata_W =  data.frame(bdo = mean_bdo, nit = mean_nit)
newdata_X =  data.frame(bdo = mean_bdo, nit = nit)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                new_X = newdata_X,
                new_W = newdata_W,
                component = 'uniqueness',
                scale_uniq = F,
                type = 'response'),
        var = 'nit',
        value = nit,
        which = 'diss. gradient only')

# diss effect of bdo (total)
newdata_W =  data.frame(bdo = mean_bdo, nit = mean_nit)
newdata_X =  data.frame(bdo = bdo, nit = mean_nit)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = F,
                     type = 'response'),
             var = 'bdo',
             value = bdo,
             which = 'diss. gradient only')

# uniq effect nit (total)
newdata_W =  data.frame(bdo = mean_bdo, nit = nit)
newdata_X =  data.frame(bdo = mean_bdo, nit = mean_nit)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = F,
                     type = 'response'),
             var = 'nit',
             value = nit,
             which = 'direct only')

# uniq effect bdo
newdata_W =  data.frame(bdo = bdo, nit = mean_nit)
newdata_X =  data.frame(bdo = mean_bdo, nit = mean_nit)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = F,
                     type = 'response'),
             var = 'bdo',
             value = bdo,
             which = 'direct only')


### PLOT PARTIAL EFFECTS ####
pred_lcbd <- do.call(rbind, pred_lcbd)

partial_eff_plot <- ggplot(pred_lcbd, aes(x = value, y = mean)) + 
  geom_point(aes(col = which), shape = 21) +
  geom_line(aes(col = which)) +
  theme(strip.placement = 'outside',
        legend.position = '',
        axis.title.x = element_blank()) +
  geom_ribbon(aes(ymin = CI.2.5., ymax = CI.97.5., col = which), fill = 'transparent', lty = 5) +
  scale_colour_manual(values = c('black', 'coral', 'steelblue')) + 
  facet_grid(which~var, scales = 'free_x', switch = 'x') + 
  ylab('Expected uniqueness (SS)')


### DISSIMILARITY GRADIENTS #### 
diss_gradients <- do.call(rbind, diss_gradient(m, CI_quant = c(0.95, 0.5)))

(fx_plot <- ggplot(diss_gradients, aes(x = x, y = f_x)) + 
  ylab('f(x)') +
  theme(strip.placement = 'outside', 
        axis.title.x = element_blank()) +
  geom_ribbon(aes(ymin = `CI 2.5%`, ymax = `CI 97.5%`), fill = 'steelblue', alpha = 0.1) +
  geom_ribbon(aes(ymin = `CI 25%`, ymax = `CI 75%`), fill = 'steelblue', alpha = 0.2) +
  geom_line(col = 'steelblue4') + 
  facet_wrap(~var, scales = 'free_x', switch = 'x'))

cowplot::plot_grid(partial_eff_plot, fx_plot, ncol = 1, rel_heights = c(3,1.5), labels = c('C', 'D'), align = 'v')

ggsave('Proposal/case_study.png', dpi = 600, height = 16, width = 10, units = 'cm')
