##------------------------------------------------------------------------------
## Imports libraries

library(ggplot2)
library(patchwork)
source("specifications_for_plots.R")

##------------------------------------------------------------------------------
## Imports values 

# q = 5
load("RData/unknown_pre_psample_q20_K3_delta0.RData")
ps_K3 = vec_p
load("RData/unknown_pre_psample_q20_K5_delta0.RData")
ps_K5 = vec_p
load("RData/unknown_pre_psample_q20_K10_delta0.RData")
ps_K10 = vec_p
# q = 5
load("RData/unknown_pre_pmed_q20_K3_delta0.RData")
pm_K3 = vec_p
load("RData/unknown_pre_pmed_q20_K5_delta0.RData")
pm_K5 = vec_p
load("RData/unknown_pre_pmed_q20_K10_delta0.RData")
pm_K10 = vec_p

##------------------------------------------------------------------------------
## Creates plot for p_{sample}

df_ps = data.frame(p = c(ps_K3, ps_K5, ps_K10), 
                   case = as.factor(rep(c("3", "5", "10"), each = num_trial)))
df_ps$case = factor(df_ps$case, levels = c("3", "5", "10"))

plot_ps = ggplot(df_ps, aes(sample = p, color = case)) + 
  stat_qq(distribution = qunif, size = 1) + 
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linewidth = 2, linetype = "dashed") + 
  labs(title = expression(p[hat(sigma)["sample"]]),
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
## Creates plot for p_{med}

df_pm = data.frame(p = c(pm_K3, pm_K5, pm_K10), 
                   case = as.factor(rep(c("3", "5", "10"), each = num_trial)))
df_pm$case = factor(df_pm$case, levels = c("3", "5", "10"))

plot_pm = ggplot(df_pm, aes(sample = p, color = case)) + 
  stat_qq(distribution = qunif, size = 1) + 
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linewidth = 2, linetype = "dashed") + 
  labs(title = expression(p[hat(sigma)["med"]]),
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

plot_comb = ((plot_ps | plot_pm) & 
  theme(legend.position = "right", 
        axis.title = element_text(size = 30, color = "black")) &
  #---
  # Citations: the line below for changing the size of points in the legend has 
  # been taken from https://stackoverflow.com/questions/20415963/how-to-
  # increase-the-size-of-points-in-legend-of-ggplot2
  guides(color = guide_legend(override.aes = list(size = 7)))) +
  #---
  plot_layout(guides = "collect")

to_save = paste0("Figures/figure_10.png")
if (dir.exists("Figures") == 0) {
  dir.create("Figures")
  ggsave(plot_comb, filename = to_save, width = 16, height = 7)
} else {
  ggsave(plot_comb, filename = to_save, width = 16, height = 7)
}

##------------------------------------------------------------------------------