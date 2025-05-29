library(ade4)
library(splines2)
library(hardhat)
library(Matrix)
library(lme4)
library(TMB)

compile('TMB/gdum.cpp')
dyn.load(dynlib('TMB/gdum')) # load TMB model

dir <- './functions/'
files.sources = list.files(dir)
sapply(paste0(dir,files.sources), source)
set.seed(1)

# Load data
data("microbialdata", package = 'gllvm')
env <- microbialdata$Xenv
com <- microbialdata$Y
env$site <- rownames(env)
hist(log(env$Phosp))

m <- gdmm(Y = com, 
          X = env, 
          diss_formula = ~isp(pH,4),
          uniq_formula = NULL,
          binary = TRUE,
          family = 'normal',
          method = 'bray',
          link = 'identity',
          bboot = F,
          n_boot = 1000,
          mono = T)

summary(m)

LCBD_list <- adespatial::LCBD.comp(vegan::vegdist(com, binary = TRUE))
LCBD <- LCBD_list$LCBD*LCBD_list$beta['SStotal']
summary(lm(LCBD~splines::ns(env$pH,2)))

round(m$sdrep$cov.fixed, 5)
n = 100
newdata_W =  data.frame(pH = seq(min(env$pH), max(env$pH), length.out = n),
                        SOM = median(env$SOM),
                        Phosp = mean(env$Phosp))

newdata_X =  data.frame(pH = rep(mean(env$pH), n),
                        SOM = rep(mean(env$SOM), n),
                        Phosp = rep(mean(env$Phosp), n))

LCBD <- adespatial::LCBD.comp(vegan::vegdist(com, binary = T), sqrt.D = T)
obs <- LCBD$LCBD*LCBD$beta[1]

pred1 <- predict_gdmm(m, type = 'link', component = 'uniqueness', conditional = TRUE) 

plot(x = obs,
     y = pred1)
abline(0,1)

pred <- predict_gdmm(m, type = 'link', component = 'dissimilarity', conditional = TRUE)

plot(x = m$Y_diss,
     y = pred)

plot(as.matrix(vegan::vegdist(com, binary = TRUE, method = 'bray'))[cbind(m$D[,1], m$D[,2])],
     pred)
abline(0,1)

