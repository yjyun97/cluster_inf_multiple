##------------------------------------------------------------------------------
## Imports libraries

library(ggplot2)
library(patchwork)
source("specifications_for_plots.R")

##------------------------------------------------------------------------------
## Imports values 

# q = 5
load("RData/unknown_pre_pstar_q5_K3_delta0.RData")
q5_K3 = vec_p
load("RData/unknown_pre_pstar_q5_K5_delta0.RData")
q5_K5 = vec_p
load("RData/unknown_pre_pstar_q5_K10_delta0.RData")
q5_K10 = vec_p
# q = 10
load("RData/unknown_pre_pstar_q10_K3_delta0.RData")
q10_K3 = vec_p
load("RData/unknown_pre_pstar_q10_K5_delta0.RData")
q10_K5 = vec_p
load("RData/unknown_pre_pstar_q10_K10_delta0.RData")
q10_K10 = vec_p
# q = 20
load("RData/unknown_pre_pstar_q20_K3_delta0.RData")
q20_K3 = vec_p
load("RData/unknown_pre_pstar_q20_K5_delta0.RData")
q20_K5 = vec_p
load("RData/unknown_pre_pstar_q20_K10_delta0.RData")
q20_K10 = vec_p

##------------------------------------------------------------------------------
## Creates plot for q = 5

df_q5 = data.frame(p = c(q5_K3, q5_K5, q5_K10), 
                   case = as.factor(rep(c("3", "5", "10"), each = num_trial)))
df_q5$case = factor(df_q5$case, levels = c("3", "5", "10"))

plot_q5 = ggplot(df_q5, aes(sample = p, color = case)) + 
  stat_qq(distribution = qunif, size = 1) + 
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linewidth = 2, linetype = "dashed") + 
  labs(title = "q = 5",
       x = "Theoretical Quantiles", 
       y = "Empirical Quantiles",
       color = "K") + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  scale_color_manual(values = c("3" = "#ef476f",
                                "5" = "#ffd166", 
                                "10" = "#26547c")) + 
  plot_details + 
  theme(legend.position = "none") 

##------------------------------------------------------------------------------
## Creates plot for q = 10

df_q10 = data.frame(p = c(q10_K3, q10_K5, q10_K10), 
                    case = as.factor(rep(c("3", "5", "10"), each = num_trial)))
df_q10$case = factor(df_q10$case, levels = c("3", "5", "10"))

plot_q10 = ggplot(df_q10, aes(sample = p, color = case)) + 
  stat_qq(distribution = qunif, size = 1) + 
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linewidth = 2, linetype = "dashed") + 
  labs(title = "q = 10",
       x = "Theoretical Quantiles", 
       y = "Empirical Quantiles",
       color = "K") + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  scale_color_manual(values = c("3" = "#ef476f",
                                "5" = "#ffd166", 
                                "10" = "#26547c")) + 
  plot_details + 
  theme(legend.position = "none") 

##------------------------------------------------------------------------------
## Creates plot for q = 20

df_q20 = data.frame(p = c(q20_K5, q20_K10, q20_K3), 
                    case = as.factor(rep(c("3", "5", "10"), each = num_trial)))
df_q20$case = factor(df_q20$case, levels = c("3", "5", "10"))

plot_q20 = ggplot(df_q20, aes(sample = p, color = case)) + 
  stat_qq(distribution = qunif, size = 1) + 
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linewidth = 2, linetype = "dashed") + 
  labs(title = "q = 20",
       x = "Theoretical Quantiles", 
       y = "Empirical Quantiles",
       color = "K") + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  scale_color_manual(values = c("3" = "#ef476f",
                                "5" = "#ffd166", 
                                "10" = "#26547c")) + 
  plot_details + 
  theme(legend.position = "none") 

##------------------------------------------------------------------------------
## Combines and saves plots

plot_comb = ((plot_q5 | plot_q10 | plot_q20) & 
  theme(legend.position = "right") &
  #---
  # Citations: the line below for changing the size of points in the legend has 
  # been taken from https://stackoverflow.com/questions/20415963/how-to-
  # increase-the-size-of-points-in-legend-of-ggplot2
  guides(color = guide_legend(override.aes = list(size = 7)))) +
  #---
  plot_layout(guides = "collect")

to_save = paste0("Figures/figure_9.png")
if (dir.exists("Figures") == 0) {
  dir.create("Figures")
  ggsave(plot_comb, filename = to_save, width = 24.5, height = 7)
} else {
  ggsave(plot_comb, filename = to_save, width = 24.5, height = 7)
}

##------------------------------------------------------------------------------