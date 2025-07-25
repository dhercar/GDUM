library(ade4)
library(splines2)
library(hardhat)
library(Matrix)
library(lme4)
library(TMB)

library(cowplot)
library(ggplot2)

theme_set(theme_bw())
theme_update(text = element_text(size = 8))


compile('TMB/gdum.cpp')
dyn.load(dynlib('TMB/gdum')) # load TMB model

dir <- './TMB/functions/'
files.sources = list.files(dir)
sapply(paste0(dir,files.sources), source)
set.seed(1)

