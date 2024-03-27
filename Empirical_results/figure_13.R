##------------------------------------------------------------------------------
## Imports libraries

library(ggplot2)
library(patchwork)
source("specifications_for_plots.R")
source("../Proposed_methods/functions_helper.R")

##------------------------------------------------------------------------------
## Specifies settings

num_sim = 1500
vec_del = seq(0, 6, 0.5)
ld = length(vec_del)
thresh = 0.05

##------------------------------------------------------------------------------
## Imports values 

load("RData/unknown_dep_pstarJ_closest_K20_g3_horizontal_delta0.RData")
null = vec_p

# horizontal
list_h = list()
for (which_test_ in c("pstarJ", "psigmaJ")) {
  for (K_ in c(20)) {
    mat_p_h = matrix(0, num_sim, ld)
    for (i_ in 1:ld) {
      if (which_test_ == "pstarJ") {
        load(paste0("RData/unknown_dep_", which_test_, "_closest_K", K_, 
                    "_g3_horizontal_delta", vec_del[i_], ".RData"))
        mat_p_h[, i_] = vec_p
      } else {
        load(paste0("RData/known_dep_", which_test_, "_closest_q20_K", K_, 
                    "_g3_horizontal_delta", vec_del[i_], ".RData"))
        mat_p_h[, i_] = vec_p
      }
    }
    list_h[[paste0(which_test_, "_K", K_)]] = mat_p_h
  }
}

# Kgon
list_K = list()
for (which_test_ in c("pstarJ", "psigmaJ")) {
  for (K_ in c(20)) {
    mat_p_K = matrix(0, num_sim, ld)
    for (i_ in 1:ld) {
      if (which_test_ == "pstarJ") {
        load(paste0("RData/unknown_dep_", which_test_, "_closest_K", K_, 
                    "_g3_Kgon_delta", vec_del[i_], ".RData"))
        mat_p_K[, i_] = vec_p
      } else {
        load(paste0("RData/known_dep_", which_test_, "_closest_q20_K", K_, 
                    "_g3_Kgon_delta", vec_del[i_], ".RData"))
        mat_p_K[, i_] = vec_p
      }
    }
    list_K[[paste0(which_test_, "_K", K_)]] = mat_p_K
  }
}

##------------------------------------------------------------------------------
## Computes empirical power

# horizontal
list_power_h = list()
list_eb_h = list()
for (which_test_ in c("pstarJ", "psigmaJ")) {
  for (K_ in c(20)) {
    res_h = apply(list_h[[paste0(which_test_, "_K", K_)]], 2, 
                  function(x) { fun_power(x, thresh) })
    list_power_h[[paste0(which_test_, "_K", K_)]] = res_h[1, ]
    list_eb_h[[paste0(which_test_, "_K", K_)]] = res_h[2, ]
  }
}

# Kgon
list_power_K = list()
list_eb_K = list()
for (which_test_ in c("pstarJ", "psigmaJ")) {
  for (K_ in c(20)) {
    res_K = apply(list_K[[paste0(which_test_, "_K", K_)]], 2, 
                  function(x) { fun_power(x, thresh) })
    list_power_K[[paste0(which_test_, "_K", K_)]] = res_K[1, ]
    list_eb_K[[paste0(which_test_, "_K", K_)]] = res_K[2, ]
  }
}

##------------------------------------------------------------------------------
## Creates plot for the null hypothesis

df = data.frame(p = null)

plot_ = ggplot(df, aes(sample = p)) + 
  stat_qq(distribution = qunif, size = 1, color = "#ef476f") + 
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linewidth = 2, linetype = "dashed") + 
  labs(title = expression(paste("QQ plot (", delta, " = 0)")), 
       x = "Theoretical Quantiles", 
       y = "Empirical Quantiles",
       color = "K") + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  plot_details + 
  guides(color = guide_legend(override.aes = list(size = 5))) 

##------------------------------------------------------------------------------
## Creates plot for horizontal

power_h = c(list_power_h$"pstarJ_K20", list_power_h$"psigmaJ_K20")
eb_h = c(list_eb_h$"pstarJ_K20", list_eb_h$"psigmaJ_K20")
df_h = data.frame(del = rep(vec_del, 2), 
                  power = power_h, 
                  case = factor(rep(c("pstarJ", "psigmaJ"), each = ld)), 
                  leb = power_h - eb_h, 
                  ueb = power_h + eb_h)
df_h$case = factor(df_h$case, levels = c("pstarJ", "psigmaJ"))
plot_h = ggplot(df_h, aes(del, power, color = case)) + 
  geom_errorbar(aes(ymin = leb, ymax = ueb, color = case), 
                width = 0.3, linewidth = 1) + 
  geom_line(aes(color = case), linewidth = 1, alpha = 0.5) + 
  geom_point(aes(color = case), size = 2) + 
  labs(title = "Horizontal", 
       x = expression(paste(delta)), 
       y = "Empirical Power", 
       color = "Test") +
  scale_color_manual(values = c("pstarJ" = "#ED2939",
                                "psigmaJ" = "black"), 
                     labels = expression(
                       "pstarJ" = paste(p, "*")[J], 
                       "psigmaJ" = p[paste(sigma, J)])) + 
  plot_details + 
  theme(legend.position = "none") 

##------------------------------------------------------------------------------
## Creates plot for Kgon

power_K = c(list_power_K$"pstarJ_K20", list_power_K$"psigmaJ_K20")
eb_K = c(list_eb_K$"pstarJ_K20", list_eb_K$"psigmaJ_K20")
df_K = data.frame(del = rep(vec_del, 2), 
                  power = power_K, 
                  case = factor(rep(c("pstarJ", "psigmaJ"), each = ld)), 
                  leb = power_K - eb_K, 
                  ueb = power_K + eb_K)
df_K$case = factor(df_K$case, levels = c("pstarJ", "psigmaJ"))
plot_K = ggplot(df_K, aes(del, power, color = case)) + 
  geom_errorbar(aes(ymin = leb, ymax = ueb, color = case), 
                width = 0.3, linewidth = 1) + 
  geom_line(aes(color = case), linewidth = 1, alpha = 0.5) + 
  geom_point(aes(color = case), size = 2) + 
  labs(title = "K-gon", 
       x = expression(paste(delta)), 
       y = "Empirical Power", 
       color = "Test") +
  scale_color_manual(values = c("pstarJ" = "#ED2939",
                                "psigmaJ" = "black"), 
                     labels = expression(
                       "pstarJ" = paste(p, "*")[J], 
                       "psigmaJ" = p[paste(sigma, J)])) + 
  plot_details + 
  theme(legend.position = "none") 

##------------------------------------------------------------------------------
## Combines and saves plots

plot_comb = ((plot_ | plot_h | plot_K) &
  theme(legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        legend.key.size = unit(2, "cm"), 
        legend.text.align = 0, 
        plot.title = element_text(hjust = 0.5, size = 40), 
        axis.title = element_text(size = 30, color = "black")) &
  #---
  # Citations: the line below for changing the size of points in the legend has 
  # been taken from https://stackoverflow.com/questions/20415963/how-to-
  # increase-the-size-of-points-in-legend-of-ggplot2
  guides(color = guide_legend(override.aes = list(size = 7)))) +
  #---
  plot_layout(guides = "collect")

to_save = paste0("Figures/figure_13.png")
if (dir.exists("Figures") == 0) {
  dir.create("Figures")
  ggsave(plot_comb, filename = to_save, width = 23.5, height = 7)
} else {
  ggsave(plot_comb, filename = to_save, width = 23.5, height = 7)
}

##------------------------------------------------------------------------------