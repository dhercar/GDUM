library(tidyverse)
library(gdmmTMB)

sp <- gdm::southwest %>% select(site, species) %>% 
  group_by(site, species) %>%
  summarise(n = 1) %>%
  pivot_wider(id_cols = site, names_from = species,values_from = n, values_fill = 0 ) %>%
  arrange(site)

env <- gdm::southwest %>% select(-c(species)) %>%
  distinct() %>% arrange(site)


D <- t(combn(nrow(env), 2))
dist = data.frame(dist = scale(dist(env[,c('bio19')]))[D])

m2 <- gdmm(Y = sp[,-1],
          X = env,
          X_pair = dist,
          D = D,
          diss_formula = ~
            isp(bio5) +
            isp(bio6) +
            isp(bio15) +
            isp(log(bio19)),
          uniq_formula = ~
            scale(bio5) +
            scale(bio6) +
            scale(bio15) +
            scale(log(bio19)),
          pair_formula = ~ dist,
          mono = T,
          mono_pair = T,
          family = 'binomial',
          n_boot = 1000,
          bboot = T)

summary(m2)
sum <- summary(m2)
sum$CI
boxplot(m$boot_samples[,1:17])
abline(h=0)
f_x <- diss_gradient(m2, CI_quant = c(0.5, 0.95))

library(ggplot2)

do.call(rbind, f_x) %>% ggplot(aes(fill = var)) +
  geom_line(aes(x = x, y = f_x)) +
  geom_ribbon(aes(x = x, ymax = `CI 2.5%`, ymin = `CI 97.5%`), alpha = 0.1) +
  geom_ribbon(aes(x = x, ymax = `CI 25%`, ymin = `CI 75%`), alpha = 0.3) +
  facet_wrap(~var, scales = 'free_x') 




new_X_mean <- data.frame(
  awcA = mean(env$awcA),
  phTotal = mean(env$phTotal),
  sandA = mean(env$sandA),
  shcA = mean(env$shcA),
  solumDepth = mean(env$solumDepth),
  bio5 = mean(env$bio5),
  bio6 = mean(env$bio6),
  bio15 = mean(env$bio15),
  bio19 = rep(mean(env$bio19), nrow(envDF)))  

plot(predict(m2, CI = F), m2$Y_diss)
abline(0,1)


ggplot(pres, aes(x = x/100000, y = y/100000, fill = mean)) + 
  geom_tile() +
  coord_equal() + 
  scale_fill_viridis_c()


pres <- cbind(envDF, predict(m, new_X = envDF, new_W = envDF, component = 'uniqueness', scale_uniq = F))

ggplot(pres, aes(x = x/100000, y = y/100000, fill = mean)) + 
  geom_tile() +
  coord_equal() + 
  scale_fill_viridis_c()

pres <- cbind(env, predict(m, new_X = env, new_W = new_X_mean, scale_uniq = F, CI = F))

ggplot(pres, aes(x = x/100000, y = y/100000, fill = mean)) + 
  geom_tile() +
  coord_equal() + 
  scale_fill_viridis_c()
