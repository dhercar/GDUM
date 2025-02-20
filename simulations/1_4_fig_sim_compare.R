library(greta)
library(latex2exp)
library(ggplot2)
library(extrafont)
library(sysfonts)
library(showtext)
library(cowplot)

# ggplot settings
theme_set(theme_bw())
theme_update(panel.grid = element_blank())
set.seed(1)

# Load functions
dir <- './functions/'
files.sources = list.files(dir)
sapply(paste0(dir,files.sources), source)

# Load models 
m1 <- readRDS('models/m1.rds')
m2_1 <- readRDS('models/m2_1.rds')
m2_2 <- readRDS('models/m2_2.rds')
m3 <- readRDS('models/m3.rds')

# ----Dissimilarity plots ----

# Expected dissimilarity at increasing |x_i - x_j| in each model

pred_m1 <- with(m1, {
  pred_data <- data.frame(X_e_dis = seq(min(X_e_dis), max(X_e_dis), length.out = n_pred))
  
  # Linear predictor
  e_s <- cbind(normal(0, SD_s, dim = n_pred), normal(0, SD_s, dim = n_pred)) # Simulate new site combinations
  D_s <- pred_data$X_e_dis %*% beta_e_dis  # Dissimilarity component
  eta <- alpha + D_s + e_s[, 1] + e_s[, 2] # Linear predictor
  mu <- ilogit(eta)
  # calculate
  pred <- data.frame(calculate(mu, values = draws, nsim = 10000))
  
  # quantiles
  pred_sum <- t(apply(pred, 2, function(x)
    quantile(
      as.numeric(x), probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
    )))
  
  cbind(pred_sum, pred_data)
})

pred_m3 <- with(m3, {
  pred_data <- data.frame(X_e_dis = seq(min(X_e_dis), max(X_e_dis), length.out = n_pred),
                          X_iso = mean(X_iso))
  # Linear predictor
  e_s <- cbind(normal(0, SD_s, dim = n_pred), normal(0, SD_s, dim = n_pred)) # Simulate new site combination
  S_s <- pred_data$X_iso # Simulated (unobserved) site effects
  D_s <- pred_data$X_e_dis %*% beta_e_dis # Linear predictor # Dissimilarity component
  eta <- alpha + D_s + e_s[, 1] + e_s[, 2] # Linear predictor
  mu <- ilogit(eta) # link
  
  # calculate
  pred <- data.frame(calculate(mu, values = draws, nsim = 10000))
  
  # quantiles
  pred_sum <- t(apply(pred, 2, function(x)
    quantile(
      as.numeric(x), probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
    )))
  
  
  
  cbind(pred_sum, pred_data)
})

# Points
pred_m3_pairs <- with(m3, {
  # PRED 1: With effect of isolation
  # Linear predictor
  D_s <- X_e_dis %*% beta_e_dis  # Dissimilarity component
  S_s <- X_iso %*% beta_iso
  eta <- alpha + D_s + S_s[row] + S_s[col] # Linear predictor
  mu <- ilogit(eta) # link
  
  # calculate
  pred1 <- data.frame(calculate(mu, values = draws, nsim = 10000))
  
  # quantiles
  pred1_sum <- t(apply(pred1, 2, function(x)
    quantile(
      as.numeric(x), probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
    )))
  
  # PRED 2: Without effect of isolation
  # Linear predictor
  D_s <- X_e_dis %*% beta_e_dis  # Dissimilarity component
  eta <- alpha + D_s 
  mu2 <- ilogit(eta) # link
  
  # calculate
  pred2 <- data.frame(calculate(mu2, values = draws, nsim = 10000))
  
  # quantiles
  pred2_sum <- t(apply(pred2, 2, function(x)
    quantile(
      as.numeric(x), probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
    )))
  
  data.frame(pred2_sum, X_e_dis = X_e_dis, Y = Y[,1]/Y[,2], y_marginal = Y[,1]/Y[,2] - (pred1_sum[,'50%'] - pred2_sum[,'50%']))
})


(m1_diss_plot <- ggplot(pred_m1, aes(x = X_e_dis, y = `50%`)) +
    geom_point(data = data.frame(X_e_dis = m1$X_e_dis, Y = m1$Y[,1]/m1$Y[,2]), aes(y = Y),
               shape = 21, colour = 'grey60') +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`,
                    fill = '95% CI'), alpha = 0.5) +
    geom_ribbon(aes(ymin = `25%`, ymax = `75%`,
                    fill = '50% CI'), alpha = 0.5) +
    scale_fill_manual('', values = c('steelblue', 'lightblue')) +
    geom_line(colour = 'darkblue', size = 1) + 
    annotate("text", -Inf, Inf, label = "m1", colour = 'steelblue', hjust = -0.5, vjust = 2, size = 3) +
    ylab(expression(bold(E)*"(Z"["ij"]*")") ) +
    xlab( TeX("$|x_i - x_j|$ (sd units)") ) + 
    ylim(c(0.185,1)) +
    theme(legend.position = '')) 

(m3_diss_plot <- ggplot(pred_m3, aes(x = X_e_dis, y = `50%`)) + 
    geom_point(data = pred_m3_pairs, aes(y = y_marginal),
               shape = 21, colour = 'grey60') +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`,
                    fill = '95% CI'), alpha = 0.5) +
    geom_ribbon(aes(ymin = `25%`, ymax = `75%`,
                    fill = '50% CI'), alpha = 0.5) +
    scale_fill_manual('', values = c('steelblue', 'lightblue')) +
    geom_line(colour = 'darkblue', size = 1) + 
    annotate("text", -Inf, Inf, label = "m3", colour = 'steelblue', hjust = -0.5, vjust = 2, size = 3) +
    ylab(expression(bold(E)*"(Z"["ij"]*")") ) +
    xlab( TeX("$|x_i - x_j|$ (sd units)") ) + 
    ylim(c(0.185,1)) +
    theme(legend.position = '')) 

# ---- Uniqueness plots ----

## m2_1 ##
# Expected uniqueness (unobserved site)
pred_u_m2_1_exp <- with(m2_1, {
  pred_data <- data.frame(X_iso = seq(min(X_iso), max(X_iso), length.out = n_pred))
  # Linear predictor
  u <- alpha/2 + pred_data$X_iso %*% beta_iso + normal(0, SD_s) # Simulate new site combinations
  # Draws
  pred <- data.frame(calculate(u, values = draws, nsim = 10000))
  # quantiles
  pred_sum <- t(apply(pred, 2, function(x)
    quantile(as.numeric(x), probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
  cbind(pred_sum, pred_data)
})

# Expected uniqueness (sites)
pred_u_m2_1 <- with(m2_1, {
  # Linear predictor
  u <- alpha/2 + X_iso %*% beta_iso + e_s # Simulate new site combinations
  # Draws
  pred <- data.frame(calculate(u, values = draws, nsim = 10000))
  # quantiles
  pred_sum <- t(apply(pred, 2, function(x)
    quantile(as.numeric(x), probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
    )
  
  data.frame(pred_sum, X_iso)
})

## m2_2 ##
# Expected uniqueness (unobserved site)
pred_u_m2_2_exp <- with(m2_2, {
  pred_data <- data.frame(X_iso = seq(min(X_iso), max(X_iso), length.out = n_pred))
  SS        <- alpha + pred_data$X_iso %*% beta_iso + normal(0, sd)
  pred      <- data.frame(calculate(SS, values = draws, nsim = 10000))
  pred_sum  <- t(apply(pred, 2, function(x)  quantile(as.numeric(x), probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
  
  cbind(pred_sum, pred_data)
})

## m3
pred_u_m3_exp <- with(m3, {
  pred_data    <- data.frame(X_e_dis = mean(X_e_dis),
                          X_iso = seq(min(X_iso), max(X_iso), length.out = n_pred))
  D_s          <- mean(X_e_dis) %*% beta_e_dis # Dissimilariy component (long)
  D_s_m        <- matrix(D_s, ncol = n_pred + 1, nrow = n_pred + 1)
  diag(D_s_m)  <-  0
  D_s_s        <- greta::rowSums(D_s_m) / nrow(D_s_m)
  u          <- alpha/2 + pred_data$X_iso %*% beta_iso + normal(0,SD_s)
  pred         <- data.frame(calculate(u, values = draws, nsim = 1000))
  pred_sum     <- t(apply(pred, 2, function(x)
    quantile(as.numeric(x), probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
  cbind(pred_sum, pred_data)
})

# Expected uniqueness (sites)
pred_u_m3 <- with(m3, {
  D_s         <- mean(X_e_dis) %*% beta_e_dis # Dissimilarity component (long)
  # Long dissimilarity component to matrix
  D_s_m       <- zeros(length(unique(row)) + 1, length(unique(col)) + 1)
  D_s_m[cbind(row,col)]  <-  D_s
  D_s_m[cbind(col,row)]  <-  D_s
  D_s_s       <- greta::rowSums(D_s_m) / nrow(D_s_m)
  u           <- alpha/2 + X_iso %*% beta_iso + D_s_s + e_s
  pred        <- data.frame(calculate(u, values = draws, nsim = 1000))
  pred_sum    <- t(apply(pred, 2, function(x)
    quantile(as.numeric(x), probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  )
  
  data.frame(pred_sum, X_iso)
})

# Plot
(m2_1_uniq_plot <- ggplot(pred_u_m2_1_exp, aes(x = X_iso, y = `50%`)) + 
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, 
                    fill = '95% CI'), alpha = 0.5) + 
    geom_ribbon(aes(ymin = `25%`, ymax = `75%`, 
                    fill = '50% CI'), alpha = 0.5) + 
    scale_fill_manual('', values = c('orange', 'wheat')) +
    geom_line(colour = 'darkorange3', size = 1) + 
    ylab(TeX("uniqueness ($u_i$)")) +
    xlab('w (sd units)') +
    ylim(-0.5, 0.9)+
    geom_point(data = pred_u_m2_1,
               aes(x = X_iso, y = `X50.`), shape = 21, col = 'grey30') + 
    annotate("text", -Inf, Inf, label = "m2.1", colour = 'darkorange', hjust = -0.5, vjust = 2, size = 3) +
    theme(legend.position = ''))

# plot
(m2_2_uniq_plot <- ggplot(pred_u_m2_2_exp, aes(x = X_iso, y = `50%`)) + 
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, 
                    fill = '95% CI'), alpha = 0.5) + 
    geom_ribbon(aes(ymin = `25%`, ymax = `75%`, 
                    fill = '50% CI'), alpha = 0.5) + 
    scale_fill_manual('', values = c('orange', 'wheat')) +
    geom_line(colour = 'darkorange3', size = 1) + 
    ylab('predicted') + 
    geom_point(data = data.frame(X_iso = m2_2$X_iso, Y= m2_2$Y), aes(x = X_iso, y = Y), shape = 21, col = 'grey30') +
    ylab(TeX("uniqueness ($SS_i$)")) +
    xlab('w (sd units)') +    
    ylim(0.135, 0.45)+
    annotate("text", -Inf, Inf, label = "m2.2", colour = 'darkorange', hjust = -0.5, vjust = 2, size = 3) +
    theme(legend.position = ''))


# Plot
(m3_uniq_plot <- ggplot(pred_u_m3_exp, aes(x = X_iso, y = `50%`)) + 
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, 
                    fill = '95% CI'), alpha = 0.5) + 
    geom_ribbon(aes(ymin = `25%`, ymax = `75%`, 
                    fill = '50% CI'), alpha = 0.5) + 
    scale_fill_manual('', values = c('orange', 'wheat')) +
    geom_line(colour = 'darkorange3', size = 1) + 
    ylab(TeX("uniqueness ( $u_{i}$)")) +
    xlab('w (sd units)') +
    ylim(-0.5, 0.9)+
    geom_point(data = pred_u_m3,
               aes(x = X_iso, y = `X50.`), shape = 21, col = 'grey30') + 
    annotate("text", -Inf, Inf, label = "m3", colour = 'darkorange', hjust = -0.5, vjust = 2, size = 3) +
    theme(legend.position = ''))


# ----  R2 dissimilarity ----
r2_m1 <-  with(m1, {
  D_s <- X_e_dis %*% beta_e_dis
  eta <- alpha + D_s
  mu <- ilogit(eta) # link
  pred <- data.frame(calculate(mu, values = draws, nsim = 1000))
  
  apply(pred, 1, function(x) {
    ypred <- x
    error <- Y[, 1] / Y[, 2] - x
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
})

r2_m2_1 <-  with(m2_1, {
  S_s <- X_iso %*% beta_iso
  eta <- alpha + S_s[row] + S_s[col]
  mu <- ilogit(eta)
  pred <- data.frame(calculate(mu, values = draws, nsim = 1000))
  
  apply(pred, 1, function(x) {
    ypred <- x
    error <- Y[, 1] / Y[, 2] - x
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
  
})

r2_m3 <-  with(m3, {
  S_s <- X_iso %*% beta_iso
  D_s <- X_e_dis %*% beta_e_dis 
  eta <- alpha + D_s + S_s[row] + S_s[col]
  mu <- ilogit(eta) 
  pred <- data.frame(calculate(mu, values = draws, nsim = 1000))
  
  apply(pred, 1, function(x) {
    ypred <- x
    error <- Y[, 1] / Y[, 2] - x
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
})


# ---- R2 uniqueness plots ----
r2_u_m1 <- with(m1, {
  D_s <- X_e_dis %*% beta_e_dis 
  D_s_m <- zeros(length(unique(row)) + 1, length(unique(col)) + 1)
  D_s_m[cbind(row,col)]  <-  D_s
  D_s_m[cbind(col,row)]  <-  D_s
  D_s_s <- greta::rowSums(D_s_m) / nrow(D_s_m)
  
  eta1 <- alpha/2 + D_s_s 
  
  eta2 <- alpha/2 + D_s_s + e_s
  
  LCBD <- cbind(eta2,
                eta1) # Predictor alone
  # calculate
  pred_u <- data.frame(calculate(LCBD, values = draws, nsim = 1000))
  apply(pred_u, 1, function(x){
    ypred <- x[51:100]  
    error <- (x[1:50] - x[51:100])
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
})

r2_u_m2_1 <- with(m2_1, {
  
  eta1 <- alpha/2 + X_iso %*% beta_iso 
  eta2 <- alpha/2 + X_iso %*% beta_iso + e_s
  
  
  LCBD <- cbind(eta2,
                eta1)

  pred_u <- data.frame(calculate(LCBD, values = draws, nsim = 1000))
  apply(pred_u, 1, function(x){
    ypred <- x[51:100]  
    error <- (x[1:50] - x[51:100])
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
})


r2_u_m2_2 <- with(m2_2, {
  LCBD <-  cbind(Y, alpha + X_iso %*% beta_iso)

  pred_u <- data.frame(calculate(LCBD, values = draws, nsim = 1000))
  apply(pred_u, 1, function(x){
    ypred <- x[51:100]  
    error <- (x[1:50] - x[51:100])
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
})

r2_u_m1 <- with(m1, {
  D_s <- X_e_dis %*% beta_e_dis
  D_s_m <- zeros(length(unique(row)) + 1, length(unique(col)) + 1)
  D_s_m[cbind(row,col)]  <-  D_s
  D_s_m[cbind(col,row)]  <-  D_s
  D_s_s <- greta::rowSums(D_s_m) / nrow(D_s_m)
  
  eta1 <- alpha/2 + D_s_s 
  
  eta2 <- alpha/2 + D_s_s + e_s
  
  LCBD <- cbind(eta2,
                eta1) 
  pred_u <- data.frame(calculate(LCBD, values = draws, nsim = 1000))
  apply(pred_u, 1, function(x){
    ypred <- x[51:100]  
    error <- (x[1:50] - x[51:100])
    var_ypred <- var(ypred)
    var_e <- var(error)
    var_ypred / (var_ypred + var_e)
  })
})

r2_data <- rbind(data.frame( r2 = r2_m1, m = 'm1', comp = 'pairwise dissimilarity'),
                 data.frame( r2 = r2_m2_1, m = 'm2.2', comp = 'pairwise dissimilarity'),
                 data.frame( r2 = r2_m3, m = 'm3', comp = 'pairwise dissimilarity'),
                 data.frame( r2 = r2_u_m1, m = 'm1', comp = 'site uniqueness'),
                 data.frame( r2 = r2_m3, m = 'm3', comp = 'site uniqueness'),
                 data.frame( r2 = r2_u_m2_1, m = 'm2.1', comp ='site uniqueness'),
                 data.frame( r2 = r2_u_m2_2, m = 'm2.2', comp ='site uniqueness'))

# ---- VARPART ----
pred_m3_full <- with(m3, {

  D_s <- X_e_dis %*% beta_e_dis 
  D_s_m <- zeros(length(unique(row)) + 1, length(unique(col)) + 1)
  D_s_m[cbind(row,col)]  <-  D_s
  D_s_m[cbind(col,row)]  <-  D_s
  D_s_s <- greta::rowSums(D_s_m) / nrow(D_s_m)
  u <- alpha/2 + X_iso %*% beta_iso + D_s_s + e_s

  pred <- data.frame(calculate(u, values = draws, nsim = 1000))

  pred_sum <- apply(pred, 2, function(x){
    quantile(as.numeric(x), probs = c(0.5))
    })
  
  pred_sum
})

pred_m3_no_x <- with(m3, {
  u <- alpha/2 + X_iso %*% beta_iso + e_s
  pred <- data.frame(calculate(u, values = draws, nsim = 1000))
  pred_sum <- apply(pred, 2, function(x){
    quantile(as.numeric(x), probs = c(0.5))
  })
  
  pred_sum
})

pred_m3_no_w <- with(m3, {
  D_s <- X_e_dis %*% beta_e_dis
  D_s_m <- zeros(length(unique(row)) + 1, length(unique(col)) + 1)
  D_s_m[cbind(row,col)]  <-  D_s
  D_s_m[cbind(col,row)]  <-  D_s
  D_s_s <- greta::rowSums(D_s_m) / nrow(D_s_m)
  u <- alpha/2 + D_s_s + e_s
  pred <- data.frame(calculate(u, values = draws, nsim = 1000))
  pred_sum <- apply(pred, 2, function(x){
    quantile(as.numeric(x), probs = c(0.5))
  })
  pred_sum
})

pred_m3_no_xw <- with(m3, {
  u <- alpha/2 + e_s
  pred <- data.frame(calculate(u, values = draws, nsim = 1000))
  pred_sum <- apply(pred, 2, function(x){
    quantile(as.numeric(x), probs = c(0.5))
  })
  pred_sum
})

var_exp = (var(pred_m3_full)) 
var_xw = (var_exp - var(pred_m3_no_xw)) / var_exp
var_x = (var_exp - var(pred_m3_no_x)) / var_exp
var_w = (var_exp - var(pred_m3_no_w)) / var_exp
var_shared = var_xw - var_x - var_w

var_comp = data.frame(component = c('01_pairwise', '03_site'),
                      name = c('x \n (pairwise)', 'w \n (site-level)'),
                      prop = c(var_x, var_w),
                      prop2 = c(paste0(round(var_x*100,1), ' %'), paste0(round(var_w*100, 1), ' %')))

# ---- FINAL PLOT ----
(sim_plot <- plot_grid(
  plot_grid(
    plot_grid(
    m1_diss_plot + theme(legend.position = '', text = element_text(size = 7), aspect.ratio = 1),
    m3_diss_plot + theme(legend.position = '', text = element_text(size = 7), aspect.ratio = 1), ncol = 2),
    
    plot_grid(
    m2_1_uniq_plot + theme(legend.position = '', text = element_text(size = 7), aspect.ratio = 1),
    m2_2_uniq_plot + theme(legend.position = '', text = element_text(size = 7), aspect.ratio = 1),
    m3_uniq_plot + theme(legend.position = '', text = element_text(size = 7), aspect.ratio = 1), ncol = 3
    ), ncol = 1, rel_heights = c(2,2), labels = c('A', 'B'), label_size = 9
  ),
  plot_grid(
  ggplot(data = var_comp, aes(x = 0, y = prop, fill = component)) +
    geom_col() + 
    theme_minimal() +
    theme(
      text = element_text(size = 7),
      aspect.ratio = 1.4,
      legend.position = '',
      panel.grid = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.y = element_line(colour = 'grey20')) +
    geom_text(aes(label = name, x = 0.85), 
              position = position_stack(vjust = 0.5), size = 2) +
    geom_text(aes(label = prop2, x =0),
               alpha = 0.7,
               col = 'grey10',
              fontface = 2,
              position = position_stack(vjust = 0.5), size =3) +
    scale_fill_manual(values = c('lightblue3', 'orange')) +
    ylab(TeX("% var ($u_i$)")) +
    scale_y_continuous(limits = c(0,01), expand = c(0,0), labels = scales::percent_format()) +
    scale_x_continuous(limits = c(-0.5,1.1)) ,
  ggplot(r2_data, aes(x = m, y = r2)) + 
    geom_violin(aes(fill = comp), col = 'transparent', alpha = 0.4, scale = 'width') + 
    facet_wrap(~ comp, scales = 'free', ncol = 2) + 
    stat_summary(fun.data = function(x){ 
      data.frame(y = quantile(x, 0.5),
                 ymin = quantile(x, 0.025),
                 ymax = quantile(x, 0.975))},
      size = 0.2, aes(col = comp)) +
    xlab('model') + 
    ylab(TeX("$R^{2}$")) +
    scale_fill_manual(values = c( 'lightblue3', 'orange')) +
    scale_colour_manual(values = c( 'darkblue', 'darkorange')) +
    theme(legend.position = '',
          aspect.ratio = 1,
          text = element_text(size = 7),
          axis.title.x = element_blank(),
          strip.background = element_blank()), ncol =2, rel_widths = c(2,3),
  labels = c('C', 'D'), label_size = 9), ncol = 1, label_size = 9, rel_heights = c(2.5,1)))

ggsave('plots/simulations.png', sim_plot, width = 15, height = 15, dpi = 600, units = 'cm')
