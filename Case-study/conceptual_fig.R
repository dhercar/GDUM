library(ggplot2)
library(dplyr)
library(ggforce)
library(cowplot)
library(patchwork)
library(adespatial)

set.seed(123)

# plotting preferences
theme_set(theme_void(base_size = 8))
theme_update(strip.text = element_text(size = 7, face = 2,  hjust = 0.5))

n_points <- 100
x_vals <- runif(n_points, 0, 1)  # continuous variable

# simulate points given a continuous x
make_points_cont <- function(x_vals, mean_fun, sd_fun, pattern) {
  tibble(
    x_var = x_vals,
    mean_x = mean_fun(x_vals),
    mean_y = 0,
    sd = sd_fun(x_vals)) %>%
    rowwise() %>%
    mutate(
      x = rnorm(1, mean_x, sd),
      y = rnorm(1, mean_y, sd),
      pattern = pattern) %>%
    ungroup()
}

# 1. Centroid shift only
df1 <- make_points_cont(
  x_vals,
  mean_fun = function(x) x*4 - 2,  # centroid shifts
  sd_fun = function(x) 0.5,        # constant dispersion
  pattern = "A1. directional change\n(centroid shifts)"
)
df1$LCBD <- adespatial::LCBD.comp(dist(df1[,c('x','y')]))$LCBD

# 2. Dispersion change only
df2 <- make_points_cont(
  x_vals,
  mean_fun = function(x) 0,
  sd_fun = function(x) 0.3 + 0.7*x, # sd increases with x
  pattern = "A2. non-directional change\n(dispersion increases)"
)

df2$LCBD <- adespatial::LCBD.comp(dist(df2[,c('x','y')]))$LCBD

# 3. Both centroid and dispersion
df3 <- make_points_cont(
  x_vals,
  mean_fun = function(x) x*4 - 2,
  sd_fun = function(x) 0.3 + 0.7*x,
  pattern = "A3. directional and\nnon-directional change"
)

df3$LCBD <- adespatial::LCBD.comp(dist(df3[,c('x','y')]))$LCBD

df <- bind_rows(df1, df2, df3)

# Select three representative x values per pattern: min, mean, max
circles <- df %>%
  group_by(pattern) %>%
  summarise(
    x_min = quantile(x_var, 0.05),
    x_mean = mean(x_var),
    x_max =  quantile(x_var, 0.95),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c(x_min, x_mean, x_max), names_to = "pos", values_to = "x_var") %>%
  rowwise() %>%
  mutate(
    mean_x = case_when(
      pattern == "A1. directional change\n(centroid shifts)" ~ x_var*4 - 2,
      pattern == "A2. non-directional change\n(dispersion increases)" ~ 0,
      pattern == "A3. directional and\nnon-directional change" ~ x_var*3 - 2
    ),
    sd = case_when(
      pattern == "A1. directional change\n(centroid shifts)" ~ 0.5,
      pattern == "A2. non-directional change\n(dispersion increases)" ~ 0.3 + 0.7*x_var,
      pattern == "A3. directional and\nnon-directional change" ~ 0.3 + 0.7*x_var
    ),
    mean_y = 0
  )

components_plot <- cowplot::plot_grid(
ggplot(df, aes(x, y, color = x_var)) +
  geom_point(alpha = 1, size = 1, shape = 21, col = 'black',stroke = 0.2, aes(fill= x_var)) +
  ggforce::geom_circle(
    data = circles,
    aes(x0 = mean_x, y0 = mean_y, r = 1.96*sd, color = x_var),
    inherit.aes = FALSE,
    linetype = 2,
    size = 0.5,
    alpha = 0.3
  ) +
  facet_wrap(~pattern, ncol = 1, strip.position = 'left') + 
  scale_color_gradientn(colours = c("#2a0000", "#a30000", "#ffb3c6")) +
  scale_fill_gradientn(colours = c("#2a0000", "#a30000", "#ffb3c6")) +
  theme(legend.position = '') +
  labs(color = "environemnt") +
  coord_equal(),

# Plot
ggplot(df, aes(x_var, LCBD, colour = x_var)) +
  facet_wrap(~pattern, ncol = 1, scales = 'free_x') +
  scale_color_gradientn(colours = c("#2a0000", "#a30000", "#ffb3c6")) +
  scale_fill_gradientn(colours = c("#2a0000", "#a30000", "#ffb3c6")) +
  theme_classic(base_size = 8) +
  theme(strip.text = element_blank(), aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(1, "lines"),
        legend.position = '') +
  labs(color = "enviroment") +
  xlab('environment') + 
  geom_smooth(method = 'lm', formula = y~poly(x,2), col = 'black', fill = 'grey90', alpha = 1, linewidth = 0.5) + 
  geom_point(alpha = 1, size = 1, shape = 21, col = 'black', aes(fill= x_var), stroke = 0.2),
ncol = 2, align  = 'h', rel_widths = c(2,1))



#########
set.seed(123)

n_points <- 200
x_vals <- runif(n_points, 0, 1)  # continuous variable
x_vals <- (x_vals - min(x_vals))/diff(range(x_vals)) *1.2
x_vals2 <- exp(x_vals*2)/max(exp(x_vals*2))
x_vals2 <- (x_vals2 - min(x_vals2))/diff(range(x_vals2)) *1.2
x_vals3 <- 1.2-x_vals2

# 1. x_vals# 1. # 1. # 1. Centroid shift only
df4 <- make_points_cont(
  x_vals,
  mean_fun = function(x) x*4 - 2,  # centroid shifts
  sd_fun = function(x) 0.5,        # constant dispersion
  pattern = "B1.even sampling"
)
df4$LCBD <- adespatial::LCBD.comp(dist(df4[,c('x','y')]))$LCBD

# 2. Dispersion change only
df5 <- make_points_cont(
  x_vals2,
  mean_fun = function(x) x*4 - 2,  # centroid shifts
  sd_fun = function(x) 0.5,        # constant dispersion
  pattern = "B2.skewed sampling"
)

df5$LCBD <- adespatial::LCBD.comp(dist(df5[,c('x','y')]))$LCBD


df <- bind_rows(df4, df5)

hist_df <- df %>%
  group_by(pattern) %>%
  do({
    h <- hist(.$x_var, breaks = 15, plot = FALSE)
    tibble(
      x = h$mids,
      count = h$counts
    )
  })

# use same custom gradient
my_cols <- c("#2a0000", "#a30000", "#ffb3c6")

p_hist <- ggplot(hist_df, aes(x, count, fill = x)) +
  geom_col(color = "white", alpha = 1, size = 0.5) +
  facet_wrap(~pattern, ncol = 2, strip.position = "top", scales = "free") +
  scale_fill_gradientn(colours = my_cols) +
  xlab("") +
  theme(
    legend.position = "",
    aspect.ratio = 0.5,
    panel.spacing = unit(1, "lines"),
    axis.title.x = element_text(),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(margin = margin(b = 4)),
    strip.background = element_rect(fill = 'transparent', colour = 'transparent'))

p_scatter <- ggplot(df, aes(x_var, LCBD, fill = x_var)) +
  facet_wrap(~pattern, ncol = 2, scales = "free") +
  scale_fill_gradientn(colours = my_cols) +
  theme_classic(base_size = 8) +
  theme(
    strip.text = element_blank(),
    aspect.ratio = 1,
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.position = ""
  ) +
  xlab("environment") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2),
              alpha = 1,
              color = "black", fill = "grey90", linewidth = 0.5) + 
  geom_point(alpha = 1, size = 1, shape = 21, stroke = 0.2)

library(cowplot)
titleA <- ggdraw() + draw_label("Directional and non-directional environmental effects\non community composition shape environment-LCBD\nrelationships simulataneusly", size = 9)
titleB <- ggdraw() + draw_label("Sampling patterns influence\nenvironment-LCBD relationship under\ndirectional community change", size = 9)


# stack histogram + scatter 
right_panel <- plot_grid(
  plot_spacer(),
  p_hist, p_scatter,
  plot_spacer(),
  ncol = 1,
  align = 'v',
  axis = "lr",
  rel_heights = c(0.7,1, 2, 0.7)  
)

# add title
right_panel <- plot_grid(
  titleB, right_panel,
  ncol = 1,
  rel_heights = c(0.15, 1)
)

# left panel (unchanged)
left_panel <- plot_grid(
  titleA, components_plot,
  ncol = 1,
  rel_heights = c(0.15, 1)
)

(final_plot <- plot_grid(
  left_panel,
  plot_spacer(),
  right_panel,
  ncol = 3,
  labels = c('A', 'B'),
  rel_widths = c(1.7, 0.1, 1)  # middle controls spacing
))

ggsave('figs/conceptual.eps', device = cairo_ps, width = 19, height = 12, units = 'cm', dpi = 1200, fallback_resolution = 600)
ggsave('figs/conceptual.png', width = 19, height = 12, units = 'cm', dpi = 1200)

