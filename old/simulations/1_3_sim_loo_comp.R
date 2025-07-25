library(greta)

dir <- './functions/'
files.sources = list.files(dir)
sapply(paste0(dir,files.sources), source)

m0 <- readRDS('models/m0.rds')
m1 <- readRDS('models/m1.rds')
m2_1 <- readRDS('models/m2_1.rds')
m2_2 <- readRDS('models/m2_2.rds')
m3 <- readRDS('models/m3.rds')
m4 <- readRDS('models/m4.rds')
m5 <- readRDS('models/m5.rds')
m6 <- readRDS('models/m6.rds')
m7 <- readRDS('models/m7.rds')

# --- Predicted values ----
mu0 <- with(m0, {
  mu <- alpha + X_e_dis %*% beta_e_dis + e_s[row] + e_s[col]
  pred_m1 <- calculate(mu, sd,  values = draws)})
mu0_no_re <- with(m0, {
  mu <- alpha + X_e_dis %*% beta_e_dis
  pred_m1 <- calculate(mu, sd,  values = draws)})

mu1 <- with(m1, {
  mu <- ilogit(alpha + X_e_dis %*% beta_e_dis+ e_s[row] + e_s[col])
  pred_m1 <- calculate(mu, values = draws)})
mu1_no_re <- with(m1, {
  mu <- ilogit(alpha + X_e_dis %*% beta_e_dis)
  pred_m1 <- calculate(mu, values = draws)})

mu2_1 <-with(m2_1, {
  mu <- ilogit(alpha + X_iso[row] %*% beta_iso + X_iso[col] %*% beta_iso + e_s[row] + e_s[col]) 
  pred_m2_1 <- calculate(mu, values = draws)})
mu2_1_no_re <-with(m2_1, {
  mu <- ilogit(alpha + X_iso[row] %*% beta_iso + X_iso[col] %*% beta_iso) 
  pred_m2_1 <- calculate(mu, values = draws)})

mu2_2_no_re <- with(m2_2, {
  mu <- alpha + X_iso %*% beta_iso
  pred_m2_2 <- calculate(mu, sd, values = draws)})

mu3 <- with(m3, {
  mu <- ilogit(alpha + X_e_dis %*% beta_e_dis + X_iso[row] %*% beta_iso + X_iso[col] %*% beta_iso + e_s[row] + e_s[col]) 
  pred_m3 <- calculate(mu, values = draws)})
mu3_no_re <- with(m3, {
  mu <- ilogit(alpha + X_e_dis %*% beta_e_dis + X_iso[row] %*% beta_iso + X_iso[col] %*% beta_iso) 
  pred_m3 <- calculate(mu, values = draws)})

mu4 <- with(m4, {
  S_s <-  X_iso %*% beta_iso + X_e_uni %*% beta_e_uni
  mu <- ilogit(alpha + S_s[col] + S_s[row] + e_s[row] + e_s[col]) 
  pred_m4 <- calculate(mu, values = draws)})
mu4_no_re  <- with(m4, {
  S_s <-  X_iso %*% beta_iso + X_e_uni %*% beta_e_uni
  mu <- ilogit(alpha + S_s[col] + S_s[row]) 
  pred_m4 <- calculate(mu, values = draws)})

mu5 <- with(m5, {
  S_s <-  X_iso %*% beta_iso + X_e_uni %*% beta_e_uni
  mu <- ilogit(alpha + S_s[col] + S_s[row] + X_e_dis %*%beta_e_dis + e_s[row] + e_s[col]) 
  pred_m5 <- calculate(mu, values = draws)})
mu5_no_re  <- with(m5, {
  S_s <-  X_iso %*% beta_iso + X_e_uni %*% beta_e_uni
  mu <- ilogit(alpha + S_s[col] + S_s[row] + X_e_dis %*%beta_e_dis) 
  pred_m5 <- calculate(mu, values = draws)})

mu6 <- with(m6, {
  S_s <-  X_iso %*% beta_iso + X_env %*% beta_env
  mu <- ilogit(alpha + S_s[col] + S_s[row] + e_s[row] + e_s[col]) 
  pred_m6 <- calculate(mu, values = draws)})
mu6_no_re <- with(m6, {
  S_s <-  X_iso %*% beta_iso + X_env %*% beta_env
  mu <- ilogit(alpha + S_s[col] + S_s[row]) 
  pred_m6 <- calculate(mu, values = draws)})

mu7 <- with(m7, {
  S_s <-  X_iso %*% beta_iso + X_env %*% beta_env
  mu <- ilogit(alpha + S_s[col] + S_s[row] + X_e_dis %*%beta_e_dis + e_s[row] + e_s[col]) 
  pred_m4 <- calculate(mu, values = draws)})
mu7_no_re  <- with(m7, {
  S_s <-  X_iso %*% beta_iso + X_env %*% beta_env
  mu <- ilogit(alpha + S_s[col] + S_s[row] + X_e_dis %*% beta_e_dis ) 
  pred_m7 <- calculate(mu, values = draws)})


# ---- R2 Z (dissimilarity) ----
(r2_m0 <- r2_greta(observed = m0$Y,
                  estimated = mu0_no_re,
                  mean = "mu",
                  summary = TRUE))

(r2_m1 <- r2_greta(observed = m1$Y[,1],
                  estimated = mu1_no_re,
                  mean = "mu",
                  trials = m1$Y[,2],
                  summary = TRUE))

(r2_m2_2 <- r2_greta(observed = m2_2$Y,
                    estimated = mu2_2_no_re,
                    mean = "mu",
                    summary = TRUE))

(r2_m3 <- r2_greta(observed = m3$Y[,1],
                  estimated = mu3_no_re,
                  mean = "mu",
                  trials = m3$Y[,2],
                  summary = TRUE))

(r2_m4 <- r2_greta(observed = m4$Y[,1],
                  estimated = mu4_no_re,
                  trials = m4$Y[,2],
                  mean = "mu",
                  summary = TRUE))

(r2_m5 <- r2_greta(observed = m5$Y[,1],
                  estimated = mu5_no_re,
                  trials = m5$Y[,2],
                  mean = "mu",
                  summary = TRUE))

(r2_m6 <- r2_greta(observed = m6$Y[,1],
                  estimated = mu6_no_re,
                  trials = m6$Y[,2],
                  mean = "mu",
                  summary = TRUE))

(r2_m7 <- r2_greta(observed = m7$Y[,1],
                  estimated = mu7_no_re,
                  trials = m7$Y[,2],
                  mean = "mu",
                  summary = TRUE))

# ---- R2 U (uniqueness) ----
(r2_u_m0 <- with(m0, {
  D_s <- X_e_dis %*% beta_e_dis
  D_s_m <- zeros(length(unique(c(row,col))), length(unique(c(row,col))))
  D_s_m[cbind(row,col)]  <-  D_s
  D_s_m[cbind(col,row)]  <-  D_s
  D_s_s <- greta::rowSums(D_s_m) / nrow(D_s_m)
  eta1 <- alpha/2 + D_s_s  + X_iso %*% beta_iso
  eta2 <- alpha/2 + D_s_s +X_iso %*% beta_iso + e_s
  u <- cbind(eta2,eta1) 
  pred_u <- data.frame(calculate(u, values = draws, nsim = 1000))
  R2 <- apply(pred_u, 1, function(x){
    ypred <- x[51:100]  
    error <- (x[1:50] - x[51:100])
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
  quantile(R2, c(0.025, 0.5, 0.975))
}))

(r2_u_m1 <- with(m1, {
  D_s <- X_e_dis %*% beta_e_dis
  D_s_m <- zeros(length(unique(c(row,col))), length(unique(c(row,col))))
  D_s_m[cbind(row,col)]  <-  D_s
  D_s_m[cbind(col,row)]  <-  D_s
  D_s_s <- greta::rowSums(D_s_m) / nrow(D_s_m)
  eta1 <- alpha/2 + D_s_s 
  eta2 <- alpha/2 + D_s_s + e_s
  u <- cbind(eta2,eta1) 
  pred_u <- data.frame(calculate(u, values = draws, nsim = 1000))
  R2 <- apply(pred_u, 1, function(x){
    ypred <- x[51:100]  
    error <- (x[1:50] - x[51:100])
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
  
  quantile(R2, c(0.025, 0.5, 0.975))
}))

(r2_u_m2_1 <- with(m2_1, {
  D_s <- 0
  D_s_m <- zeros(length(unique(c(row,col))), length(unique(c(row,col))))
  D_s_m[cbind(row,col)]  <-  D_s
  D_s_m[cbind(col,row)]  <-  D_s
  D_s_s <- greta::rowSums(D_s_m) / nrow(D_s_m)
  eta1 <- alpha/2 + D_s_s + X_iso %*% beta_iso
  eta2 <- alpha/2 + D_s_s + X_iso %*% beta_iso + e_s
  u <- cbind(eta2,eta1) 
  pred_u <- data.frame(calculate(u, values = draws, nsim = 1000))
  R2 <- apply(pred_u, 1, function(x){
    ypred <- x[51:100]  
    error <- (x[1:50] - x[51:100])
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
  quantile(R2, c(0.025, 0.5, 0.975))
}))

(r2_u_m2_2 <- with(m2_2, {
  LCBD <-  cbind(Y, alpha + X_iso %*% beta_iso)
  # calculate
  pred_u <- data.frame(calculate(LCBD, values = draws, nsim = 1000))
  R2 <- apply(pred_u, 1, function(x){
    ypred <- x[51:100]  
    error <- (x[1:50] - x[51:100])
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
  
  quantile(R2, c(0.025, 0.5, 0.975))
}))

(r2_u_m3 <- with(m3, {
  D_s <- X_e_dis %*% beta_e_dis
  D_s_m <- zeros(length(unique(c(row,col))), length(unique(c(row,col))))
  D_s_m[cbind(row,col)]  <-  D_s
  D_s_m[cbind(col,row)]  <-  D_s
  D_s_s <- greta::rowSums(D_s_m) / nrow(D_s_m)
  eta1 <- alpha/2 + D_s_s  + X_iso %*% beta_iso
  eta2 <- alpha/2 + D_s_s +X_iso %*% beta_iso + e_s
  u <- cbind(eta2,eta1) 
  pred_u <- data.frame(calculate(u, values = draws, nsim = 1000))
  R2 <- apply(pred_u, 1, function(x){
    ypred <- x[51:100]  
    error <- (x[1:50] - x[51:100])
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
  quantile(R2, c(0.025, 0.5, 0.975))
}))

(r2_u_m4 <- with(m4, {
  D_s <- 0
  D_s_m <- zeros(length(unique(c(row,col))), length(unique(c(row,col))))
  D_s_m[cbind(row,col)]  <-  D_s
  D_s_m[cbind(col,row)]  <-  D_s
  D_s_s <- greta::rowSums(D_s_m) / nrow(D_s_m)
  eta1 <- alpha/2 + D_s_s  + X_iso %*% beta_iso + X_e_uni%*%beta_e_uni
  eta2 <- alpha/2 + D_s_s +X_iso %*% beta_iso + X_e_uni%*%beta_e_uni + e_s
  u <- cbind(eta2,eta1) 
  pred_u <- data.frame(calculate(u, values = draws, nsim = 1000))
  R2 <- apply(pred_u, 1, function(x){
    ypred <- x[51:100]  
    error <- (x[1:50] - x[51:100])
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
  quantile(R2, c(0.025, 0.5, 0.975))
}))

(r2_u_m5 <- with(m5, {
  D_s <- X_e_dis %*% beta_e_dis
  D_s_m <- zeros(length(unique(c(row,col))), length(unique(c(row,col))))
  D_s_m[cbind(row,col)]  <-  D_s
  D_s_m[cbind(col,row)]  <-  D_s
  D_s_s <- greta::rowSums(D_s_m) / nrow(D_s_m)
  eta1 <- alpha/2 + D_s_s  + X_iso %*% beta_iso + X_e_uni%*%beta_e_uni
  eta2 <- alpha/2 + D_s_s +X_iso %*% beta_iso + X_e_uni%*%beta_e_uni + e_s
  u <- cbind(eta2,eta1) 
  pred_u <- data.frame(calculate(u, values = draws, nsim = 1000))
  R2 <- apply(pred_u, 1, function(x){
    ypred <- x[51:100]  
    error <- (x[1:50] - x[51:100])
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
  quantile(R2, c(0.025, 0.5, 0.975))
}))

(r2_u_m6 <- with(m6, {
  D_s <- 0
  D_s_m <- zeros(length(unique(c(row,col))), length(unique(c(row,col))))
  D_s_m[cbind(row,col)]  <-  D_s
  D_s_m[cbind(col,row)]  <-  D_s
  D_s_s <- greta::rowSums(D_s_m) / nrow(D_s_m)
  eta1 <- alpha/2 + D_s_s  + X_iso %*% beta_iso + X_env%*%beta_env
  eta2 <- alpha/2 + D_s_s +X_iso %*% beta_iso + X_env%*%beta_env + e_s
  u <- cbind(eta2,eta1) 
  pred_u <- data.frame(calculate(u, values = draws, nsim = 1000))
  R2 <- apply(pred_u, 1, function(x){
    ypred <- x[51:100]  
    error <- (x[1:50] - x[51:100])
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
  quantile(R2, c(0.025, 0.5, 0.975))
}))

(r2_u_m7 <- with(m7, {
  D_s <- X_e_dis %*% beta_e_dis
  D_s_m <- zeros(length(unique(c(row,col))), length(unique(c(row,col))))
  D_s_m[cbind(row,col)]  <-  D_s
  D_s_m[cbind(col,row)]  <-  D_s
  D_s_s <- greta::rowSums(D_s_m) / nrow(D_s_m)
  eta1 <- alpha/2 + D_s_s  + X_iso %*% beta_iso + X_env%*%beta_env
  eta2 <- alpha/2 + D_s_s +X_iso %*% beta_iso + X_env%*%beta_env + e_s
  u <- cbind(eta2,eta1) 
  pred_u <- data.frame(calculate(u, values = draws, nsim = 1000))
  R2 <- apply(pred_u, 1, function(x){
    ypred <- x[51:100]  
    error <- (x[1:50] - x[51:100])
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
  quantile(R2, c(0.025, 0.5, 0.975))
}))

# ---- loo ----
(loo_m0 <- loo_greta(observed = m0$Y, 
                    posterior = mu0_no_re, 
                    mean = "mu",
                    scale = 'sd',
                    chains = 4,
                    moment_match = FALSE,
                    family = "normal", 
                    method = "loo"))

(loo_m1 <- loo_greta(observed = m1$Y[,1], 
          posterior = mu1_no_re, 
          trials = m1$Y[,2],
          mean = "mu",
          chains = 4,
          moment_match = FALSE,
          family = "binomial", 
          method = "loo"))

(loo_m2_1 <- loo_greta(observed = m2_1$Y[,1], 
                      posterior = mu2_1_no_re, 
                      trials = m2_1$Y[,2],
                      mean = "mu",
                      chains = 4,
                      moment_match = FALSE,
                      family = "binomial", 
                      method = "loo"))

(loo_m2_2 <- loo_greta(observed = m2_2$Y, 
                    posterior = mu2_2_no_re,
                    mean = "mu",
                    chains = 4,
                    moment_match = FALSE,
                    scale = 'sd',
                    method = "loo"))

(loo_m3 <- loo_greta(observed = m3$Y[,1], 
                    posterior = mu3_no_re, 
                    trials = m3$Y[,2],
                    mean = "mu",
                    chains = 4,
                    moment_match = FALSE,
                    family = "binomial", 
                    method = "loo"))

(loo_m4 <- loo_greta(observed = m4$Y[,1], 
                    posterior = mu4_no_re, 
                    trials = m4$Y[,2],
                    mean = "mu",
                    chains = 4,
                    moment_match = FALSE,
                    family = "binomial", 
                    method = "loo"))

(loo_m5 <- loo_greta(observed = m5$Y[,1], 
                    posterior = mu5_no_re, 
                    trials = m5$Y[,2],
                    mean = "mu",
                    chains = 4,
                    moment_match = FALSE,
                    family = "binomial", 
                    method = "loo"))

(loo_m6 <- loo_greta(observed = m6$Y[,1], 
                    posterior = mu6_no_re, 
                    trials = m6$Y[,2],
                    mean = "mu",
                    chains = 4,
                    moment_match = FALSE,
                    family = "binomial", 
                    method = "loo"))

(loo_m7 <- loo_greta(observed = m7$Y[,1], 
                    posterior = mu7_no_re, 
                    trials = m7$Y[,2],
                    mean = "mu",
                    chains = 4,
                    moment_match = FALSE,
                    family = "binomial", 
                    method = "loo"))

# ---- RMSE ----
RMSE_m0 <- RMSE_greta(observed = m0$Y, 
                       posterior = mu0_no_re, 
                       mean = "mu", 
                      summary = TRUE)

RMSE_m1 <- RMSE_greta(observed = m1$Y[,1], 
                      posterior = mu1_no_re, 
                      trials = m1$Y[,2],
                      mean = "mu",
                      summary = TRUE)

RMSE_m2_1 <- RMSE_greta(observed = m2_1$Y[,1], 
                      posterior = mu2_1_no_re, 
                      trials = m2_1$Y[,2],
                      mean = "mu",
                      summary = TRUE)

RMSE_m2_2 <- RMSE_greta(observed = m2_2$Y, 
                        posterior = mu2_2_no_re, 
                        mean = "mu",
                        summary = TRUE)

RMSE_m3 <- RMSE_greta(observed = m3$Y[,1], 
                      posterior = mu3_no_re, 
                      trials = m3$Y[,2],
                      mean = "mu",
                      summary = TRUE)

RMSE_m4 <- RMSE_greta(observed = m4$Y[,1], 
                      posterior = mu4_no_re, 
                      trials = m4$Y[,2],
                      mean = "mu",
                      summary = TRUE)


RMSE_m5 <- RMSE_greta(observed = m5$Y[,1], 
                      posterior = mu5_no_re, 
                      trials = m5$Y[,2],
                      mean = "mu",
                      summary = TRUE)


RMSE_m6 <- RMSE_greta(observed = m6$Y[,1], 
                      posterior = mu6_no_re, 
                      trials = m6$Y[,2],
                      mean = "mu",
                      summary = TRUE)


RMSE_m7 <- RMSE_greta(observed = m7$Y[,1], 
                      posterior = mu7_no_re, 
                      trials = m7$Y[,2],
                      mean = "mu",
                      summary = TRUE)
