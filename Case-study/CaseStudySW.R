library(tidyverse)
library(gdmmTMB)

data <- gdm::southwest


sp <- data %>% select(site, species) %>% 
  group_by(site, species) %>%
  summarise(n = 1) %>%
  pivot_wider(id_cols = site, names_from = species,values_from = n, values_fill = 0 ) %>%
  arrange(site)

env <- data %>%  select(-species) %>%
  distinct() %>% arrange(site)


m <- gdmm(Y = sp[,-1],
          X = env,
          diss_formula = ~ 
            isp(awcA) + 
            isp(phTotal) + 
            isp(sandA) + 
            isp(shcA) +
            isp(solumDepth) + 
            isp(bio5) + 
            isp(bio6) + 
            isp(bio15) + 
            isp(log(bio19)), 
          uniq_formula = ~ 
            scale(awcA) + 
            scale(phTotal) +
            scale(sandA) + 
            scale(shcA) + 
            scale(solumDepth) +
            scale(bio5) + 
            scale(bio6) + 
            scale(bio15) + 
            scale(log(bio19)) + (1|site),
          mono= F,
          family = 'beta',
          replace_01  = c(0.001, 0.999),
          n_boot = 500,
          bboot = F)

summary(m)

f_x <- diss_gradient(m)

library(ggplot2)

do.call(rbind, f_x) %>% ggplot() +
  geom_line(aes(x = x, y = f_x)) +
  geom_ribbon(aes(x = x, ymax = `CI 2.5%`, ymin = `CI 97.5%`), alpha = 0.1) +
  facet_wrap(~var, scales = 'free_x') 
