##  Load libraries
library(greta)
library(tidyverse)
library(vegan)
library(bayesplot)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(latex2exp)
library(sf)

# Convenience functions
dir <- './Manuscript/Functions/'
files.sources = list.files(dir)
sapply(paste0(dir,files.sources), source)

theme_set(theme_bw())
theme_update(panel.grid = element_blank())
set.seed(1)



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

mLCBDrds <-  tools::file_path_sans_ext(list.files("./case_study_atlanticforest",
                                                    pattern = ".rds",
                                                    full.names = TRUE))

mLCBD <- lapply(paste0(mLCBDrds, '.rds'), readRDS)
DiatDraws <- mLCBD[[1]]

u_site <- summary(mLCBD[[1]])$statistics[2:101,'Mean']
u_catch <- summary(mLCBD[[1]])$statistics[103:122,'Mean']

map_data <- data.frame(catchment = rep(1:20, each = 5), 
                       site = 1:100,
                       u_site = u_site,
                       u_catch = u_catch[rep(1:20, each = 5)],
                       lat = DiattEnv$lat,
                       lon = DiattEnv$lon)

diat <- data.frame(summary(DiatDraws)$quantiles)
sum_diat_r <- diat[substr(row.names(diat), 1,3) == 'e_s',]
sum_diat_r$r <- 1:100
sum_diat_r$c <- rep(1:20, each = 5)

sum_diat_c <- diat[substr(row.names(diat), 1,3) == 'e_c',]
sum_diat_c$c <- 1:20

#  No Hydroshed added
diat_u_plot <-
  plot_grid(
    plot_grid(
      ggplot(map_data, aes(x = lon, y = lat)) +
        theme_minimal() +
        theme(axis.text = element_blank()) +
        coord_fixed() +
        xlab('longitude') +
        ylab('latitude') +
        theme(legend.position = '') +
        ggforce::geom_mark_hull(
          aes(group = as.factor(catchment), fill = u_catch),
          concavity = 1,
          alpha = 1,
          label.fontsize = 6,
          expand = unit(2, "mm"),
          radius = unit(2, "mm")
        ) +
        scale_fill_gradientn(
          colours = c('coral4','coral4', 'coral', 'white', 'lightblue3', 'steelblue4', 'steelblue4'),
          breaks = c(-2, -1, 0, 1, 2),
          limits = c(-2, 2)
        ) +
        geom_point(aes(fill = u_site), shape = 21, size = 2) +
        scale_y_continuous(expand = c(0.05,0.05)),

      plot_grid(
        map_data %>%
          ggplot(aes(x = u_catch, y = 1)) +
          theme_classic() +
          theme(
            axis.line.y = element_blank(),
            legend.position = '',
            plot.title = element_text(hjust = 0.5, size = 7),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()
          ) +
          ggridges::geom_density_ridges_gradient(aes(fill = after_stat(x))) +
          scale_fill_gradientn(
            colours = c('coral4','coral4', 'coral', 'white', 'lightblue3', 'steelblue4', 'steelblue4'),
            breaks = c(-2, -1, 0, 1, 2),
            limits = c(-2, 2)
          ) +
          ggtitle('Catchmnet') +
          scale_x_continuous(TeX("$\\alpha_{c}$", italic = T), limits = c(-1.5, 2.5)) +
          scale_y_continuous('', expand = c(0, 0)),
        map_data %>%
          ggplot(aes(x = u_site, y = 1)) +
          theme_classic() +
          theme(
            axis.line.y = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 7),
            axis.ticks.y = element_blank(),
            legend.position = '',
            axis.text.y = element_blank()
          ) +
          ggridges::geom_density_ridges_gradient(aes(fill = after_stat(x))) +
          scale_fill_gradientn(
            colours = c('coral4','coral4', 'coral', 'white', 'lightblue3', 'steelblue4', 'steelblue4'),
            breaks = c(-2, -1, 0, 1, 2),
            limits = c(-2, 2)
          ) +
          ggtitle('Stream') +
          scale_x_continuous(TeX("$\\alpha_{r}$", italic = T), limits = c(-1.5, 2.5)) +
          scale_y_continuous('', expand = c(0, 0)),
        ncol = 2
      ),
      rel_heights = c(3, 2),
      labels = c('A', 'B'),
      ncol = 1
    ),
    plot_grid(
      ggplot(sum_diat_c, aes(x = X50., y = c)) +
        geom_vline(xintercept = 0,lty = 5,size = 0.3) +
        theme(text = element_text(size = 7), plot.title = element_text(hjust = 0.5)) +
        geom_linerange(aes(xmin = X2.5., xmax = X97.5.), colour = 'grey60') +
        geom_linerange(aes(xmin = X25., xmax = X75.), colour = 'black') +
        geom_point(size = 0.7) +
        scale_y_continuous('Catchemnt', breaks = 1:20) +
        xlab(TeX("$\\alpha_{c}$", italic = T))+
        ggtitle('Catchment'),

      ggplot(sum_diat_r, aes(x = X50., y = r)) +
        geom_vline(xintercept = 0, lty = 5) +
        theme(text = element_text(size = 7), plot.title = element_text( hjust = 0.5)) +
        geom_linerange(aes(xmin = X2.5., xmax = X97.5.), colour = 'grey60') +
        geom_linerange(aes(xmin = X25., xmax = X75.), colour = 'black') +
        geom_point(size = 0.7) +
        scale_y_continuous('Stream', breaks = (1:20)*5,minor_breaks = 1:100) +
        xlab(TeX("$\\alpha_{r}$", italic = T)) +
        ggtitle('Stream')
    ),

    rel_widths = c(3, 2),
    labels = c('', 'C'),
    ncol = 2
  )



ggsave('plots/diat_plot.png', diat_u_plot,width = 18, height = 10, dpi = 600, units = 'cm')  

# Streams and catchments with CI above or below 0
sum(sign(sum_diat_r$X2.5.)*sign(sum_diat_r$X97.5.) == 1)
sum(sign(sum_diat_c$X2.5.)*sign(sum_diat_c$X97.5.) == 1)


