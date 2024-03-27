##------------------------------------------------------------------------------
## Imports libraries

library(ggplot2)
library(patchwork)
source("specifications_for_plots.R")

##------------------------------------------------------------------------------
## Imports values 

# farthest
load("RData/known_dep_psigma_farthest_K20_g1_delta0.RData")
farthest = vec_p
#closest
load("RData/known_dep_psigma_closest_K20_g1_delta0.RData")
closest = vec_p

##------------------------------------------------------------------------------
## Creates plots 

df_f = data.frame(p = farthest)
plot_f = ggplot(df_f, aes(sample = p)) + 
  stat_qq(distribution = qunif, size = 1, color = "#ef476f") + 
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linewidth = 2, linetype = "dashed") + 
  labs(title = "Farthest", 
       x = "Theoretical Quantiles", 
       y = "Empirical Quantiles") + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  plot_details 

df_c = data.frame(p = closest)
plot_c = ggplot(df_c, aes(sample = p)) + 
  stat_qq(distribution = qunif, size = 1, color = "#ef476f") + 
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linewidth = 2, linetype = "dashed") + 
  labs(title = "Closest", 
       x = "Theoretical Quantiles", 
       y = "Empirical Quantiles") + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  plot_details 

##------------------------------------------------------------------------------
## Combines and saves plots

plot_comb = (plot_f | plot_c) &
  theme(plot.title = element_text(size = 40)) 

to_save = paste0("Figures/figure_2.png")
if (dir.exists("Figures") == 0) {
  dir.create("Figures")
  ggsave(plot_comb, filename = to_save, width = 15.5, height = 7)
} else {
  ggsave(plot_comb, filename = to_save, width = 15.5, height = 7)
}

##------------------------------------------------------------------------------