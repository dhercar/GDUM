##  Load libraries
library(coda)
library(greta)
library(tidyverse)
library(vegan)
library(bayesplot)
library(ggplot2)
library(ggnewscale)

# functions
dir <- './functions/'
files.sources = list.files(dir)
sapply(paste0(dir,files.sources), source)


##  Read in data
##  File descriptions for this project
##  diatoms20_5.txt - diatom abundances in streams and catchments
DiatAbund <- read.table("./case_study_atlanticforest/Schneck_2022_RiverCatchmentLCBD/diatoms20_5.txt",
                        header = TRUE)

DiattEnv <-  read.table("./case_study_atlanticforest/Schneck_2022_RiverCatchmentLCBD/env20_5.txt",
                        header = TRUE)

##  Transform to presence/absence data
DiatPA <- cbind(DiatAbund[,c(1,2)],
                decostand(DiatAbund[3:ncol(DiatAbund)],
                          method = "pa"))

##  Rename stream IDs
DiatPA$stream <- 1:100

##  Running model for diatoms
DiatDiss <- make_y_df(DiatPA[, 3:ncol(DiatPA)],
                      id = DiatPA[,2],
                      method = "sorensen",
                      num_den = TRUE)


DiatDiss$s1 <- DiatPA$stream[DiatDiss$s1]
DiatDiss$s2 <- DiatPA$stream[DiatDiss$s2]
DiatDiss$c1 <- DiatPA$catchment[DiatDiss$s1]
DiatDiss$c2 <- DiatPA$catchment[DiatDiss$s2]

# Response
Y <- DiatDiss[, "num_sor"]
Trials <- DiatDiss[, "den_sor"]

# Dimensions
n_s <- length(unique(DiatPA$stream))
n_c <- length(unique(DiatPA$catchment))

##  Sampling design
s1 <- DiatDiss$s1
s2 <- DiatDiss$s2
c1 <- DiatDiss$c1
c2 <- DiatDiss$c2

##  MODEL
##  Setting uninformative priors
alpha <- normal(0, 0.5)

SD_s <- normal(0, 1,
               truncation = c(0, Inf))

SD_c <- normal(0, 1,
               truncation = c(0, Inf))

# Precision of beta-binomial distribution
phi <- exponential(0.01)

##  Random effects for catchment (e_c) and site (e_s) standardized by SD of each level
e_s <- normal(0, 1, dim = c(n_s, 1)) * SD_s
e_c <- normal(0, 1, dim = c(n_c, 1)) * SD_c

##  Building model
eta <- alpha + e_s[s1] + e_c[c1] + e_s[s2] + e_c[c2]
mu <- ilogit(eta)

# reparametrise mu phi to shape1 shape2
shape1 <- mu*phi
shape2 <- (1 - mu)*phi

##  Defining the distribution of the data 
distribution(Y) <- beta_binomial(size = Trials, alpha = shape1,  beta = shape2)
DiatMod <- model(alpha, e_s, SD_s, e_c, SD_c, phi)

# Sampling
DiatDraws <- greta::mcmc(DiatMod,
                         chains = 4,
                         warmup = 10000,
                         n_samples = 15000,
                         one_by_one = TRUE,
                         sampler = hmc(Lmin = 10,
                                       Lmax = 15))
# Checks
bayesplot::mcmc_trace(DiatDraws, regex_pars = c('alpha', 'SD', 'phi'))
bayesplot::mcmc_intervals(DiatDraws,regex_pars = c('e_c', 'e_s'))
bayesplot::mcmc_rank_overlay(DiatDraws,regex_pars = c('e_c', 'e_s'))
gelman.diag(DiatDraws, autoburnin = FALSE, multivariate = FALSE)

saveRDS(DiatDraws,
        "./case_study_atlanticforest/Diatom_mLCBD.rds")

##  Diatom site-level LCBD calculation
DiatLCBDSite <- as.data.frame(cbind(DiatPA[, 2],
                                    adespatial::beta.div(DiatPA[, 3:ncol(DiatPA)],
                                                         method = "sorensen")$LCBD))
saveRDS(DiatLCBDSite,
        "./case_study_atlanticforest/DiatomSiteLCBD.rds")

##  Diatom catchment-level LCBD calculation
DiatCatch <- DiatPA %>%
  group_by(catchment) %>%
  summarize(across(.fns = sum,
                   .cols = 3:ncol(.)-1))
DiatCatchPA <- cbind(DiatCatch[, 1], decostand(DiatCatch[, 2:ncol(DiatCatch)],
                                               method = "pa"))
DiatLCBDCatch <- as.data.frame(cbind(DiatCatchPA[, 1],
                                     adespatial::beta.div(DiatCatchPA[, 3:ncol(DiatCatchPA)],
                                                          method = "sorensen")$LCBD))
saveRDS(DiatLCBDCatch,
        "./case_study_atlanticforest/DiatomCatchmentLCBD.rds")


