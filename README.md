GDUM: Simultaneusly modelling community dissimilarity and uniqueness
================

This repository contains the code associated with the manuscript
Hernandez-Carrasco et al. (2025).

- `functions` contains custom helper functions used in other scripts.
- `simulations`, `case_study_doubs`, and `case_study_atlanticforest`
  contain scripts to replicate the analyses presented in the manuscript
- `SI` contains code to replicate the figures in Supplementary
  Information
- `models` contains the models fitted to simulated data
- `plots` contains the figures in the main text

Below we include a short tutorial showing how to fitting GDUMs using the
wrapping functions provided in this repository.

NOTE: For transparency, our scripts do not use some of the helper
functions described below. Instead, we provide the full code used to fit
each model and generate predictions.

# Fitting GDUMs in `greta`

Here we show how to fit and inspect Generalized Dissimilarity Uniqueness
Models (DGUM) using `greta`. We provide some helper functions to
construct the models and generate predictions.This scripts assumes that
a folder named `functions` exists in the working directory with the
helper functions provided
[here](https://github.com/dhercar/GDUM/tree/main/Functions).

``` r
# Load functions
dir <- './functions/'
files.sources = list.files(dir)
sapply(paste0(dir,files.sources), source)
```

Some other necessary packages:

``` r
library(greta)
library(coda)
library(bayesplot)
library(ggplot2)
library(tidyverse)
```

Note that the package `greta` requires some additional dependencies.
Please check [greta-stats.org/](https://greta-stats.org/) for more
information.

## Case-study: Fish from the doubs river

We will use the [doubs
dataset](https://www.davidzeleny.net/anadat-r/doku.php/en:data:doubs) as
an example. The dataset contains abundances of fish species along the
Doubs river.

``` r
# Load datasets
doubs.spe <- read.csv('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/DoubsSpe.csv', 
                      row.names = 1)
doubs.env <- read.csv('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/DoubsEnv.csv', 
                      row.names = 1)
doubs.spa <- read.csv('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/DoubsSpa.csv', 
                      row.names = 1)
```

doubs.spe

| CHA | TRU | VAI | LOC | OMB | BLA | HOT | TOX | VAN | CHE | BAR | SPI | GOU | BRO | PER | BOU | PSO | ROT | CAR | TAN | BCO | PCH | GRE | GAR | BBO | ABL | ANG |
|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|----:|
|   0 |   3 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |
|   0 |   5 |   4 |   3 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |
|   0 |   5 |   5 |   5 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   1 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |
|   0 |   4 |   5 |   5 |   0 |   0 |   0 |   0 |   0 |   1 |   0 |   0 |   1 |   2 |   2 |   0 |   0 |   0 |   0 |   1 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |
|   0 |   2 |   3 |   2 |   0 |   0 |   0 |   0 |   5 |   2 |   0 |   0 |   2 |   4 |   4 |   0 |   0 |   2 |   0 |   3 |   0 |   0 |   0 |   5 |   0 |   0 |   0 |
|   0 |   3 |   4 |   5 |   0 |   0 |   0 |   0 |   1 |   2 |   0 |   0 |   1 |   1 |   1 |   0 |   0 |   0 |   0 |   2 |   0 |   0 |   0 |   1 |   0 |   0 |   0 |

doubs.env

|  das | alt |  pen |  deb |  pH | dur |  pho |  nit |  amm |  oxy | dbo |
|-----:|----:|-----:|-----:|----:|----:|-----:|-----:|-----:|-----:|----:|
|  0.3 | 934 | 48.0 | 0.84 | 7.9 |  45 | 0.01 | 0.20 | 0.00 | 12.2 | 2.7 |
|  2.2 | 932 |  3.0 | 1.00 | 8.0 |  40 | 0.02 | 0.20 | 0.10 | 10.3 | 1.9 |
| 10.2 | 914 |  3.7 | 1.80 | 8.3 |  52 | 0.05 | 0.22 | 0.05 | 10.5 | 3.5 |
| 18.5 | 854 |  3.2 | 2.53 | 8.0 |  72 | 0.10 | 0.21 | 0.00 | 11.0 | 1.3 |
| 21.5 | 849 |  2.3 | 2.64 | 8.1 |  84 | 0.38 | 0.52 | 0.20 |  8.0 | 6.2 |
| 32.4 | 846 |  3.2 | 2.86 | 7.9 |  60 | 0.20 | 0.15 | 0.00 | 10.2 | 5.3 |

doubs.spa

|   x |   y |
|----:|----:|
|  88 |   7 |
|  94 |  14 |
| 102 |  18 |
| 100 |  28 |
| 106 |  39 |
| 112 |  51 |

## Preparing data

We provide two functions to help formatting the data: `make_x_df()` and
`make_y_df()`. `make_x_df()` generates a matrix of pairwise distances
for the variables of interest:

``` r
X <- make_x_df(env = doubs.env, 
               scale = FALSE, # Whether variables should be scaled prior to computing distances
               method = 'euclidean', # Method passed to `dist`
               collapse = FALSE)  # Whether environmental variables should be collapsed into a single environmental distance

kableExtra::kable(head(X))
```

|     |  s1 |  s2 | dist_das | dist_alt | dist_pen | dist_deb | dist_pH | dist_dur | dist_pho | dist_nit | dist_amm | dist_oxy | dist_dbo |
|:----|----:|----:|---------:|---------:|---------:|---------:|--------:|---------:|---------:|---------:|---------:|---------:|---------:|
| 2   |   2 |   1 |      1.9 |        2 |     45.0 |     0.16 |     0.1 |        5 |     0.01 |     0.00 |     0.10 |      1.9 |      0.8 |
| 3   |   3 |   1 |      9.9 |       20 |     44.3 |     0.96 |     0.4 |        7 |     0.04 |     0.02 |     0.05 |      1.7 |      0.8 |
| 4   |   4 |   1 |     18.2 |       80 |     44.8 |     1.69 |     0.1 |       27 |     0.09 |     0.01 |     0.00 |      1.2 |      1.4 |
| 5   |   5 |   1 |     21.2 |       85 |     45.7 |     1.80 |     0.2 |       39 |     0.37 |     0.32 |     0.20 |      4.2 |      3.5 |
| 6   |   6 |   1 |     32.1 |       88 |     44.8 |     2.02 |     0.0 |       15 |     0.19 |     0.05 |     0.00 |      2.0 |      2.6 |
| 7   |   7 |   1 |     36.5 |       93 |     41.4 |     3.16 |     0.2 |       43 |     0.06 |     0.05 |     0.00 |      1.1 |      0.5 |

`make_y_df()` computes pairwise dissimilarities and returns a data frame
in long format. Different options are available, including:

- Methods available in `vegan::vegdist()`:
  `make_y_df(data, method = 'bray')`
- Numerator-denominator versions of Sørensen-Dice and Jaccard indices:
  `make_y_df(data, method = 'jaccard', num_den = TRUE)`
- Full matrix of shared and unshared species between samples:
  `make_y_df(data, method = 'abcd')`

We generate different matrices to show how GDUM can handle different
response variables.

``` r
Y_bin <- make_y_df(com = doubs.spe, method = 'sorensen', num_den = TRUE)
kableExtra::kable(head(Y_bin))
```

|     |  s1 |  s2 | num_sor | den_sor |
|:----|----:|----:|--------:|--------:|
| 2   |   2 |   1 |       2 |       4 |
| 3   |   3 |   1 |       3 |       5 |
| 4   |   4 |   1 |       7 |       9 |
| 5   |   5 |   1 |      10 |      12 |
| 6   |   6 |   1 |       9 |      11 |
| 7   |   7 |   1 |       4 |       6 |

``` r
Y_bray <- make_y_df(com = doubs.spe, trans = log1p)
```

    ## Warning in vegan::vegdist(com, method = method): you have empty rows: their dissimilarities may be
    ##                  meaningless in method "bray"

``` r
kableExtra::kable(head(Y_bray))
```

|     |  s1 |  s2 |      diss |
|:----|----:|----:|----------:|
| 2   |   2 |   1 | 0.5509095 |
| 3   |   3 |   1 | 0.6280761 |
| 4   |   4 |   1 | 0.7446012 |
| 5   |   5 |   1 | 0.8664653 |
| 6   |   6 |   1 | 0.7657643 |
| 7   |   7 |   1 | 0.6519273 |

Outputs from `make_y_df()` and `make_x_df()` contain two columns that
indicate which pair of samples are compared in each row: `s1` and `s2`.
These columns should be consistent across datasets.

``` r
all.equal(X[,c('s1','s2')],Y_bin[,c('s1','s2')])
```

    ## [1] TRUE

## Fitting models

We provide a function `fit_gdum()` to help built GDUMs in `greta`. By
deafault, `fit_gdum()` uses uninformative priors that expect predictors
to be scaled. In this example, we use nitrogen (nit) concentration as a
pairwise predictor and pH as a site-level effect.

``` r
X$s_dist_nit <- scale(X$dist_nit)
doubs.env$s_pH <- scale(doubs.env$pH)
```

In addition, a design matrix `D` indicating which sites are compared in
each row is required to match dissimilarities with pairwise- and
site-level predictor.

``` r
D <- X[,c('s1','s2')]
```

### Gaussian distribution

``` r
m_gaus <- fit_gdum(Y = Y_bray[,3], # Dissimilarities
         W = doubs.env, # Site level predictors
         X = X, # Pairwise level predictors
         D = D, # Design matrix 
         family = 'gaussian',
         link = 'identity',
         diss_formula = ~ s_dist_nit, # pairwise level formula 
         site_formula = ~ s_pH, # site level formula
         warmup = 2000, # mcmc warm up
         n_samples = 4000, # samples taken after warm up
         Lmin = 15, 
         Lmax = 20)
```

`Lmin` and `Lmax` are parameters passed to `greta::hmc()`. See
[greta-stats.org/reference/samplers](https://greta-stats.org/reference/samplers)
for more information.

We can explore MCMC chains using functions in `coda` and `bayesplot`.
MCMC chains are in `m_gaus$draws`.

``` r
bayesplot::mcmc_trace(m_gaus$draws, regex_pars = c('alpha','beta','lambda'))
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

- `alpha` is the global intercept
- `beta_X` are the coefficients for pairwise-level predictors
- `labmda_X` are the coefficients for site-level predictors

In addition, `coda::gelman.diag()` can be used to compute the
[Gelman-Rubin](https://en.wikipedia.org/wiki/Gelman-Rubin_statistic)
statistic

``` r
coda::gelman.diag(m_gaus$draws)
```

    ## Potential scale reduction factors:
    ## 
    ##                 Point est. Upper C.I.
    ## alpha                    1       1.00
    ## SD_s                     1       1.00
    ## sigma                    1       1.01
    ## beta_s_dist_nit          1       1.00
    ## lambda_s_pH              1       1.00
    ## e_s[1,1]                 1       1.00
    ## e_s[2,1]                 1       1.00
    ## e_s[3,1]                 1       1.00
    ## e_s[4,1]                 1       1.00
    ## e_s[5,1]                 1       1.00
    ## e_s[6,1]                 1       1.00
    ## e_s[7,1]                 1       1.00
    ## e_s[8,1]                 1       1.00
    ## e_s[9,1]                 1       1.00
    ## e_s[10,1]                1       1.00
    ## e_s[11,1]                1       1.00
    ## e_s[12,1]                1       1.00
    ## e_s[13,1]                1       1.00
    ## e_s[14,1]                1       1.01
    ## e_s[15,1]                1       1.00
    ## e_s[16,1]                1       1.00
    ## e_s[17,1]                1       1.00
    ## e_s[18,1]                1       1.00
    ## e_s[19,1]                1       1.00
    ## e_s[20,1]                1       1.00
    ## e_s[21,1]                1       1.00
    ## e_s[22,1]                1       1.01
    ## e_s[23,1]                1       1.00
    ## e_s[24,1]                1       1.00
    ## e_s[25,1]                1       1.00
    ## e_s[26,1]                1       1.00
    ## e_s[27,1]                1       1.00
    ## e_s[28,1]                1       1.00
    ## e_s[29,1]                1       1.00
    ## e_s[30,1]                1       1.00
    ## 
    ## Multivariate psrf
    ## 
    ## 1.01

In our case, all values are relatively small ($\geq 1.01$) indicating
that chains have succesfully converged.

The `summary` of the draws provides the estimated value of each model
parameter and credibility intervals

``` r
summary(m_gaus$draws)[[1]][1:5,]
```

    ##                        Mean          SD     Naive SE Time-series SE
    ## alpha           0.629882142 0.051003526 4.032183e-04   5.968936e-04
    ## SD_s            0.135112214 0.020880990 1.650787e-04   2.676229e-04
    ## sigma           0.187578585 0.006761673 5.345572e-05   8.728534e-05
    ## beta_s_dist_nit 0.166954539 0.013105565 1.036086e-04   1.678528e-04
    ## lambda_s_pH     0.008071676 0.026123290 2.065227e-04   2.397642e-04

- `sigma` is the standard deviation of the normally distributed
  residuals
- `SD_s` is the standard deviation associated to site-level random
  effects

### Other distributions

#### Beta distribution

Note that the beta distribution does not allow for 0s or 1s. For this
example, we replace those values with `0.001` and `0.999` respectively.

``` r
Y_bray$diss2 <- ifelse(Y_bray$diss == 0, 0.001, Y_bray$diss)
Y_bray$diss2 <- ifelse(Y_bray$diss2 == 1, 0.999, Y_bray$diss2)

m_beta <- fit_gdum(Y = Y_bray$diss2, 
         W = doubs.env, 
         X = X, 
         D = D, 
         family = 'beta',
         link = 'logit',
         diss_formula = ~ s_dist_nit,
         site_formula = ~ s_pH,
         warmup = 2000,
         n_samples = 4000,
         Lmin = 15,
         Lmax = 20)
```

``` r
summary(m_beta$draws)[[1]][1:5,]
```

    ##                       Mean         SD     Naive SE Time-series SE
    ## alpha           0.65830955 0.28777637 0.0022750719   0.0032740946
    ## SD_s            0.81462620 0.12025762 0.0009507200   0.0015741483
    ## phi             4.97920309 0.35145430 0.0027784902   0.0081228368
    ## beta_s_dist_nit 0.81246401 0.06458139 0.0005105607   0.0006183887
    ## lambda_s_pH     0.04346403 0.15216096 0.0012029380   0.0014451408

Instead of `sigma`, `phi` is the scale parameter of the beta
distribution. The larger the value of `phi`, the narrower the
distribution gets around the estimated dissimilarity.

#### Binomial distribution

The binomial distribution can be used when the dissimilarity index
corresponds to a proportion of shared and unshared items, such as the
Sørensen–Dice, Jaccard or Bray-Curtis dissimilarity indices. Both the
numerator and the denominator need to be supplied independently:

``` r
m_bin <- fit_gdum(Y = Y_bin$num_sor, # Dissimilarity (numerator: b + c)
         Y_den = Y_bin$den_sor,     # Dissimilarity (denominator: 2a + b + c)
         W = doubs.env, # Site level predictors
         X = X, # Pairwise level predictors
         D = D, # Design matrix 
         family = 'binomial',
         link = 'logit',
         diss_formula = ~ s_dist_nit,
         site_formula = ~ s_pH,
         warmup = 2000,
         n_samples = 4000,
         Lmin = 10,
         Lmax = 15)
```

``` r
summary(m_bin$draws)[[1]][1:4,]
```

    ##                       Mean         SD     Naive SE Time-series SE
    ## alpha           0.57446869 0.41928481 0.0033147374   0.0256615531
    ## SD_s            1.32487658 0.20056531 0.0015856080   0.0077188083
    ## beta_s_dist_nit 0.95495893 0.03989445 0.0003153933   0.0006919637
    ## lambda_s_pH     0.06518421 0.24365061 0.0019262272   0.0116601073

The binomial distribution does not contain any scale parameter.

#### Beta binomial distribution

The beta-binomial distribution can be useful when the data is
overdispersed compared to the binomial distribution.

``` r
m_bbin <- fit_gdum(Y = Y_bin$num_sor, # Dissimilarity matrix (numerator: b + c)
         Y_den = Y_bin$den_sor,     
         W = doubs.env, # Site level predictors
         X = X, # Pairwise level predictors
         D = D, # Design matrix 
         family = 'betabinomial',
         link = 'logit',
         diss_formula = ~ s_dist_nit,
         site_formula = ~ s_pH,
         warmup = 2000,
         n_samples = 4000,
         Lmin = 10,
         Lmax = 15)
```

    ## running 4 chains simultaneously on up to 24 cores

    ## 

    ##     warmup                                           0/2000 | eta:  ?s              warmup =                                        50/2000 | eta:  1m              warmup ==                                      100/2000 | eta: 41s              warmup ===                                     150/2000 | eta: 36s              warmup ====                                    200/2000 | eta: 32s              warmup =====                                   250/2000 | eta: 29s              warmup ======                                  300/2000 | eta: 28s              warmup =======                                 350/2000 | eta: 26s              warmup ========                                400/2000 | eta: 25s              warmup =========                               450/2000 | eta: 24s              warmup ==========                              500/2000 | eta: 23s              warmup ==========                              550/2000 | eta: 22s              warmup ===========                             600/2000 | eta: 21s              warmup ============                            650/2000 | eta: 21s              warmup =============                           700/2000 | eta: 20s              warmup ==============                          750/2000 | eta: 19s              warmup ===============                         800/2000 | eta: 18s              warmup ================                        850/2000 | eta: 17s | <1% bad    warmup =================                       900/2000 | eta: 17s | <1% bad    warmup ==================                      950/2000 | eta: 16s | <1% bad    warmup ===================                    1000/2000 | eta: 15s | <1% bad    warmup ====================                   1050/2000 | eta: 14s | <1% bad    warmup =====================                  1100/2000 | eta: 13s | <1% bad    warmup ======================                 1150/2000 | eta: 13s | <1% bad    warmup =======================                1200/2000 | eta: 12s | <1% bad    warmup ========================               1250/2000 | eta: 11s | <1% bad    warmup =========================              1300/2000 | eta: 10s | <1% bad    warmup ==========================             1350/2000 | eta: 10s | <1% bad    warmup ===========================            1400/2000 | eta:  9s | <1% bad    warmup ============================           1450/2000 | eta:  8s | <1% bad    warmup ============================           1500/2000 | eta:  7s | <1% bad    warmup =============================          1550/2000 | eta:  7s | <1% bad    warmup ==============================         1600/2000 | eta:  6s | <1% bad    warmup ===============================        1650/2000 | eta:  5s | <1% bad    warmup ================================       1700/2000 | eta:  4s | <1% bad    warmup =================================      1750/2000 | eta:  4s | <1% bad    warmup ==================================     1800/2000 | eta:  3s | <1% bad    warmup ===================================    1850/2000 | eta:  2s | <1% bad    warmup ====================================   1900/2000 | eta:  1s | <1% bad    warmup =====================================  1950/2000 | eta:  1s | <1% bad    warmup ====================================== 2000/2000 | eta:  0s | <1% bad
    ##   sampling                                           0/4000 | eta:  ?s            sampling                                          50/4000 | eta:  1m            sampling =                                       100/4000 | eta:  1m            sampling =                                       150/4000 | eta:  1m            sampling ==                                      200/4000 | eta:  1m            sampling ==                                      250/4000 | eta:  1m            sampling ===                                     300/4000 | eta:  1m            sampling ===                                     350/4000 | eta:  1m            sampling ====                                    400/4000 | eta:  1m            sampling ====                                    450/4000 | eta:  1m            sampling =====                                   500/4000 | eta:  1m            sampling =====                                   550/4000 | eta:  1m            sampling ======                                  600/4000 | eta:  1m            sampling ======                                  650/4000 | eta:  1m            sampling =======                                 700/4000 | eta:  1m            sampling =======                                 750/4000 | eta:  1m            sampling ========                                800/4000 | eta:  1m            sampling ========                                850/4000 | eta: 50s            sampling =========                               900/4000 | eta: 49s            sampling =========                               950/4000 | eta: 49s            sampling ==========                             1000/4000 | eta: 48s            sampling ==========                             1050/4000 | eta: 47s            sampling ==========                             1100/4000 | eta: 47s            sampling ===========                            1150/4000 | eta: 46s            sampling ===========                            1200/4000 | eta: 45s            sampling ============                           1250/4000 | eta: 45s            sampling ============                           1300/4000 | eta: 44s            sampling =============                          1350/4000 | eta: 43s            sampling =============                          1400/4000 | eta: 42s            sampling ==============                         1450/4000 | eta: 42s            sampling ==============                         1500/4000 | eta: 41s            sampling ===============                        1550/4000 | eta: 40s            sampling ===============                        1600/4000 | eta: 39s            sampling ================                       1650/4000 | eta: 38s            sampling ================                       1700/4000 | eta: 37s            sampling =================                      1750/4000 | eta: 36s            sampling =================                      1800/4000 | eta: 36s            sampling ==================                     1850/4000 | eta: 35s            sampling ==================                     1900/4000 | eta: 34s            sampling ===================                    1950/4000 | eta: 33s            sampling ===================                    2000/4000 | eta: 32s            sampling ===================                    2050/4000 | eta: 31s            sampling ====================                   2100/4000 | eta: 31s            sampling ====================                   2150/4000 | eta: 30s            sampling =====================                  2200/4000 | eta: 29s            sampling =====================                  2250/4000 | eta: 28s            sampling ======================                 2300/4000 | eta: 27s            sampling ======================                 2350/4000 | eta: 26s            sampling =======================                2400/4000 | eta: 26s            sampling =======================                2450/4000 | eta: 25s            sampling ========================               2500/4000 | eta: 24s            sampling ========================               2550/4000 | eta: 23s            sampling =========================              2600/4000 | eta: 23s            sampling =========================              2650/4000 | eta: 22s            sampling ==========================             2700/4000 | eta: 21s            sampling ==========================             2750/4000 | eta: 20s            sampling ===========================            2800/4000 | eta: 19s            sampling ===========================            2850/4000 | eta: 19s            sampling ============================           2900/4000 | eta: 18s            sampling ============================           2950/4000 | eta: 17s            sampling ============================           3000/4000 | eta: 16s            sampling =============================          3050/4000 | eta: 15s            sampling =============================          3100/4000 | eta: 14s            sampling ==============================         3150/4000 | eta: 14s            sampling ==============================         3200/4000 | eta: 13s            sampling ===============================        3250/4000 | eta: 12s            sampling ===============================        3300/4000 | eta: 11s            sampling ================================       3350/4000 | eta: 10s            sampling ================================       3400/4000 | eta: 10s            sampling =================================      3450/4000 | eta:  9s            sampling =================================      3500/4000 | eta:  8s            sampling ==================================     3550/4000 | eta:  7s            sampling ==================================     3600/4000 | eta:  6s            sampling ===================================    3650/4000 | eta:  6s            sampling ===================================    3700/4000 | eta:  5s            sampling ====================================   3750/4000 | eta:  4s            sampling ====================================   3800/4000 | eta:  3s            sampling =====================================  3850/4000 | eta:  2s            sampling =====================================  3900/4000 | eta:  2s            sampling ====================================== 3950/4000 | eta:  1s            sampling ====================================== 4000/4000 | eta:  0s

``` r
bayesplot::mcmc_trace(m_bbin$draws, regex_pars = c('alpha','beta','lambda'))
```

![](README_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
summary(m_bbin$draws)[[1]][1:5,]
```

    ##                       Mean        SD     Naive SE Time-series SE
    ## alpha           0.43412739 0.4032296 0.0031878099   0.0119662633
    ## SD_s            1.13597682 0.1833565 0.0014495600   0.0051332733
    ## phi             7.18023169 0.7186433 0.0056813738   0.0070814105
    ## beta_s_dist_nit 0.90343739 0.0726593 0.0005744222   0.0006357859
    ## lambda_s_pH     0.06069787 0.2061461 0.0016297282   0.0035196688

## Predictions and partial effects

### Assessing compositional uniqueness

Community uniqueness (u) values, [equivalent to $SS_i$ in the
computation of LCBD](https://doi.org/10.1111/ele.12141), can be obtained
from the fitted model using the `predict_gdum()` function.

``` r
u_pred <- predict_gdum(m_gaus, response = 'uniqueness')
```

We can compare those with the ones obtained with `adespatial::` and map
the values across the river:

``` r
u_pred <- data.frame(u_pred, lcbd = adespatial::LCBD.comp(vegan::vegdist(log1p(doubs.spe)))$LCBD)
```

    ## Registered S3 methods overwritten by 'adegraphics':
    ##   method         from
    ##   biplot.dudi    ade4
    ##   kplot.foucart  ade4
    ##   kplot.mcoa     ade4
    ##   kplot.mfa      ade4
    ##   kplot.pta      ade4
    ##   kplot.sepan    ade4
    ##   kplot.statis   ade4
    ##   scatter.coa    ade4
    ##   scatter.dudi   ade4
    ##   scatter.nipals ade4
    ##   scatter.pco    ade4
    ##   score.acm      ade4
    ##   score.mix      ade4
    ##   score.pca      ade4
    ##   screeplot.dudi ade4

    ## Registered S3 method overwritten by 'spdep':
    ##   method   from
    ##   plot.mst ape

    ## Registered S3 methods overwritten by 'adespatial':
    ##   method             from       
    ##   plot.multispati    adegraphics
    ##   print.multispati   ade4       
    ##   summary.multispati ade4

    ## Warning in vegan::vegdist(log1p(doubs.spe)): you have empty rows: their dissimilarities may be
    ##                  meaningless in method "bray"

``` r
ggplot(u_pred, aes(x = lcbd, y = `X50.`/sum(`X50.`))) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  coord_equal() + 
  geom_abline(lty = 5) +
  ylab('u/sum(u) (GDUM)') +
  xlab('LCBD (adespatial)') +
  geom_linerange(aes(ymax = `X97.5.`/sum(`X50.`), ymin = `X2.5.`/sum(`X50.`)), linewidth = 1, col = 'lightblue') +
  geom_linerange(aes(ymax = `X75.`/sum(`X50.`), ymin = `X25.`/sum(`X50.`)), linewidth = 2, col = 'steelblue4') + 
  geom_point(shape = 21, fill = 'white', size = 2)
```

![](README_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
ggplot(cbind(u_pred, doubs.spa) ,aes( x = x, y = y)) + 
  theme_void() +
  coord_equal() +
  geom_path() + 
  geom_point(shape = 21, aes(size = X50.), fill = 'steelblue') +
  scale_size_continuous('u_i', range = c(1,7))
```

![](README_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

Similarly, we can obtain the expected $u$ values when the site-level
random effect is set to 0, indicating the expected site uniqueness given
the combination of predictors in each site.

``` r
u_pred <- predict_gdum(m_gaus, response = 'uniqueness', re = FALSE)
u_pred <- data.frame(u_pred)
ggplot(cbind(u_pred, doubs.spa) ,aes( x = x, y = y)) + 
  theme_void() +
  coord_equal() +
  geom_path() + 
  geom_point(shape = 21, aes(size = X50.), fill = 'darkorange') +
  scale_size_continuous('u_i', range = c(1,7))
```

![](README_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

We can also control for nitrogen concentration (`s_dist_nit`) to asses
the importance of this pairwise variable in determining uniqueness.

``` r
u_pred <- predict_gdum(m_gaus, response = 'uniqueness', D_new = D, W_new = doubs.env, re = TRUE) 
u_pred <- data.frame(u_pred)
ggplot(cbind(u_pred, doubs.spa) ,aes( x = x, y = y)) + 
  theme_void() +
  coord_equal() +
  geom_path() + 
  geom_point(shape = 21, aes(size = X50.), fill = 'darkorange') +
  scale_size_continuous('u_i', range = c(1,7))
```

![](README_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

### Predict with new data

It is often useful to visualize the effect of individual predictors.
This can also be achieved with `predict_gdum()`.

``` r
# generate new data
n = 50
D_new <- t(combn(n, 2)) # n pairwise combinations
X_new <- data.frame(s_dist_nit = seq(min(X$s_dist_nit), 
                                     max(X$s_dist_nit), 
                                     length.out = nrow(D_new))) # Simulate pairwise predictors
W_new <- data.frame(s_pH = rep(mean(doubs.env$s_pH), n)) # Simulate site-level predictor

# Make predictions
pred <- predict_gdum(fit = m_bbin, X_new = X_new, W_new = W_new, D_new = D_new, 
                     re =FALSE, samples = 5000)
ggplot(cbind(pred, X_new), aes(x = s_dist_nit,y = `50%`)) +
  theme_bw() + 
  ylab('expected dissimilarity') + 
  geom_ribbon(aes(ymin = `2.5%`,max = `97.5%`), fill = 'lightblue', alpha = 0.5) + 
  geom_ribbon(aes(ymin = `25%`,max = `75%`), fill = 'steelblue4', alpha = 0.5) + 
  geom_line(colour = 'white', linewidth = 1.5) +
  geom_line()
```

![](README_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
# generate new data
n = 50
D_new <- t(combn(n, 2)) # n pairwise combinations
X_new <- data.frame(s_dist_nit = rep(mean(doubs.env$s_pH),
                                     length.out = nrow(D_new))) # Simulate pairwise predictors
W_new <- data.frame(s_pH = seq(min(doubs.env$s_pH),max(doubs.env$s_pH), length.out = n)) # Simulate site-level predictor

# Make predictions
pred <- predict_gdum(fit = m_bbin, X_new = X_new, W_new = W_new, D_new = D_new, 
                     re = FALSE, samples = 5000, response = 'uniqueness')

ggplot(cbind(pred, W_new), aes(x = s_pH,y = `50%`)) +
  theme_bw() + 
  ylab('expected uniqueness (u_i)') + 
  geom_ribbon(aes(ymin = `2.5%`,max = `97.5%`), fill = 'lightblue', alpha = 0.5) + 
  geom_ribbon(aes(ymin = `25%`,max = `75%`), fill = 'steelblue4', alpha = 0.5) + 
  geom_line(colour = 'white', linewidth = 1.5) +
  geom_line()
```

![](README_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->
