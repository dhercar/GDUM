library(tidyverse)
library(vegan)

sp <- read.table("./Case-study/Schneck_2022_RiverCatchmentLCBD/inverts20_5.txt",
                header = TRUE)
env <-  read.table("./Case-study/Schneck_2022_RiverCatchmentLCBD/env20_5.txt",
                        header = TRUE)

env$site <- as.factor(paste0(env$catchment, '_', env$stream))
env$catchment <- as.factor(env$catchment)
env$soft <- (env$mud + env$sand + env$gravel)


m <- gdmm(Y = sp[,-c(1:2)],
          X = env, 
          diss_formula = ~ isp(cover_site, degree = 2) + isp(log(cond), degree = 2) + isp(temp, degree = 2), # In the future: s( bdo ) + s(nit)
          uniq_formula = ~ scale(cover_site) + scale(log(cond)) + scale(temp)  + (1|catchment) + (1|site),
          binary = T,
          family = 'binomial',
          method = 'bray',
          mono = T)

summary(m)


n = 100
cover_site = sort(env$cover_site)
temp = sort(env$temp)
cond = sort(env$cond)

mean_cov <- rep(median(env$cover_site), n)
mean_tem <- rep(median(env$temp), n)
mean_con <- rep(median(env$cond), n)

pred_lcbd <- list()

# TOTAL effect of cover (total)
newdata_W =  data.frame(cover_site = cover_site, temp = mean_tem, cond = mean_con)
newdata_X =  data.frame(cover_site = cover_site, temp = mean_tem,  cond = mean_con)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = F,
                     type = 'response'),
             var = 'cover_site',
             value = cover_site,
             which = ' full')

# TOTAL effect of temp (total)
newdata_W =  data.frame(cover_site = mean_cov, temp = temp, cond = mean_con)
newdata_X =  data.frame(cover_site = mean_cov, temp = temp,  cond = mean_con)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = F,
                     type = 'response'),
             var = 'temp',
             value = temp,
             which = ' full')

# TOTAL effect of cond (total)
newdata_W =  data.frame(cover_site = mean_cov, temp = mean_tem, cond = cond)
newdata_X =  data.frame(cover_site = mean_cov, temp = mean_tem,  cond = cond)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = F,
                     type = 'response'),
             var = 'cond',
             value = cond,
             which = ' full')

# diss effect of cover (total)
newdata_W =  data.frame(cover_site = mean_cov, temp = mean_tem, cond = mean_con)
newdata_X =  data.frame(cover_site = cover_site, temp = mean_tem,  cond = mean_con)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = F,
                     type = 'response'),
             var = 'cover_site',
             value = cover_site,
             which = 'diss. gradient only')

# diss effect of temp (total)
newdata_W =  data.frame(cover_site = mean_cov, temp = mean_tem, cond = mean_con)
newdata_X =  data.frame(cover_site = mean_cov, temp = temp,  cond = mean_con)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = F,
                     type = 'response'),
             var = 'temp',
             value = temp,
             which = 'diss. gradient only')

# diss effect of cond (total)
newdata_W =  data.frame(cover_site = mean_cov, temp = mean_tem, cond = mean_con)
newdata_X =  data.frame(cover_site = mean_cov, temp = mean_tem,  cond = cond)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = F,
                     type = 'response'),
             var = 'cond',
             value = cond,
             which = 'diss. gradient only')

# uniq effect cover (total)
newdata_W =  data.frame(cover_site = cover_site, temp = mean_tem, cond = mean_con)
newdata_X =  data.frame(cover_site = mean_cov, temp = mean_tem,  cond = mean_con)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = F,
                     type = 'response'),
             var = 'cover_site',
             value = cover_site,
             which = 'direct only')

# uniq effect temp  (total)
newdata_W =  data.frame(cover_site = mean_cov, temp = temp, cond = mean_con)
newdata_X =  data.frame(cover_site = mean_cov, temp = mean_tem,  cond = mean_con)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = F,
                     type = 'response'),
             var = 'temp',
             value = temp,
             which = 'direct only')

# uniq effect cond  (total)
newdata_W =  data.frame(cover_site = mean_cov, temp = mean_tem, cond = cond)
newdata_X =  data.frame(cover_site = mean_cov, temp = mean_tem,  cond = mean_con)

pred_lcbd[[length(pred_lcbd) + 1]] <- 
  data.frame(predict(m,
                     new_X = newdata_X,
                     new_W = newdata_W,
                     component = 'uniqueness',
                     scale_uniq = F,
                     type = 'response'),
             var = 'cond',
             value = cond,
             which = 'direct only')


### PLOT PARTIAL EFFECTS ####
pred_lcbd <- do.call(rbind, pred_lcbd)

(partial_eff_plot <- ggplot(pred_lcbd, aes(x = value, y = mean)) + 
  geom_point(aes(col = which), shape = 21) +
  geom_line(aes(col = which)) +
  theme(strip.placement = 'outside',
        legend.position = '',
        axis.title.x = element_blank()) +
  geom_ribbon(aes(ymin = CI.2.5., ymax = CI.97.5., col = which), fill = 'transparent', lty = 5) +
  scale_colour_manual(values = c('black', 'coral', 'steelblue')) + 
  facet_grid(which~var, scales = 'free_x', switch = 'x') + 
  ylab('Expected uniqueness (SS)'))


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


