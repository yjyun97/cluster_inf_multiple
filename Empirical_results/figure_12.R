##------------------------------------------------------------------------------
## Imports libraries

library(ggplot2)
library(patchwork)
source("specifications_for_plots.R")
source("../proposed_methods/functions_helper.R")

##------------------------------------------------------------------------------
## Specifies settings

num_sim = 1500
vec_del = seq(0, 6, 0.5)
ld = length(vec_del)
thresh = 0.05

##------------------------------------------------------------------------------
## Imports values 

# horizontal
list_h = list()
for (which_test_ in c("pstar", "psigma")) {
  for (K_ in c(3, 5, 10)) {
    mat_p_h = matrix(0, num_sim, ld)
    for (i_ in 1:ld) {
      if (which_test_ == "pstar") {
        load(paste0("RData/unknown_pre_", which_test_, "_K", K_, 
                    "_horizontal_delta", vec_del[i_], ".RData"))
        mat_p_h[, i_] = vec_p
      } else {
        load(paste0("RData/known_pre_", which_test_, "_q20_K", K_, 
                    "_horizontal_delta", vec_del[i_], ".RData"))
        mat_p_h[, i_] = vec_p
      }
    }
    list_h[[paste0(which_test_, "_K", K_)]] = mat_p_h
  }
}

# Kgon
list_K = list()
for (which_test_ in c("pstar", "psigma")) {
  for (K_ in c(3, 5, 10)) {
    mat_p_K = matrix(0, num_sim, ld)
    for (i_ in 1:ld) {
      if (which_test_ == "pstar") {
        load(paste0("RData/unknown_pre_", which_test_, "_K", K_, 
                    "_Kgon_delta", vec_del[i_], ".RData"))
        mat_p_K[, i_] = vec_p
      } else {
        load(paste0("RData/known_pre_", which_test_, "_q20_K", K_, 
                    "_Kgon_delta", vec_del[i_], ".RData"))
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
for (which_test_ in c("pstar", "psigma")) {
  for (K_ in c(3, 5, 10)) {
    res_h = apply(list_h[[paste0(which_test_, "_K", K_)]], 2, 
                  function(x) { fun_power(x, thresh) })
    list_power_h[[paste0(which_test_, "_K", K_)]] = res_h[1, ]
    list_eb_h[[paste0(which_test_, "_K", K_)]] = res_h[2, ]
  }
}

# Kgon
list_power_K = list()
list_eb_K = list()
for (which_test_ in c("pstar", "psigma")) {
  for (K_ in c(3, 5, 10)) {
    res_K = apply(list_K[[paste0(which_test_, "_K", K_)]], 2, 
                  function(x) { fun_power(x, thresh) })
    list_power_K[[paste0(which_test_, "_K", K_)]] = res_K[1, ]
    list_eb_K[[paste0(which_test_, "_K", K_)]] = res_K[2, ]
  }
}

##------------------------------------------------------------------------------
## Creates plot for K = 3 (horizontal)

power_K3_h = c(list_power_h$"pstar_K3", list_power_h$"psigma_K3")
eb_K3_h = c(list_eb_h$"pstar_K3", list_eb_h$"psigma_K3")
df_K3_h = data.frame(del = rep(vec_del, 2), 
                     power = power_K3_h, 
                     case = factor(rep(c("pstar", "psigma"), 
                                       each = ld)), 
                     leb = power_K3_h - eb_K3_h, 
                     ueb = power_K3_h + eb_K3_h)
df_K3_h$case = factor(df_K3_h$case, 
                      levels = c("pstar", "psigma"))
#---
# Citations: the codes below, as well as the lines followed by #---
# for adding two vertical titles for the y-axis is adapted from 
# https://stackoverflow.com/questions/74534630/ggplot2-two-different-font-
# sizes-on-y-axis 
yaxis_K3_h = list("names" = c("Horizontal", "Empirical power"))
yaxis_K3_h = lapply(yaxis_K3_h, function(x) {
  x[[1]] = paste0("<span style='font-size: 40pt'>", x[[1]], "</span>")
  return(x)
})
#---
plot_K3_h = ggplot(df_K3_h, 
                   aes(del, power, color = case)) + 
  geom_errorbar(aes(ymin = leb, ymax = ueb, color = case), 
                width = 0.3, linewidth = 1) + 
  geom_line(aes(color = case), linewidth = 1, alpha = 0.5) + 
  geom_point(aes(color = case), size = 2) + 
  labs(title = "K = 3", 
       x = expression(paste(delta)), 
       y = paste(yaxis_K3_h[["names"]], collapse = "<br>"), #---
       color = "Test") +
  ylim(0, 1) + 
  scale_color_manual(values = c("pstar" = "#ED2939",
                                "psigma" = "black"), 
                     labels = expression(
                       "pstar" = paste(p, "*"), 
                       "psigma" = p[sigma])) + 
  plot_details + 
  theme(legend.position = "none") + 
  theme(axis.title.y = ggtext::element_markdown()) #---

##------------------------------------------------------------------------------
## Creates plot for K = 5 (horizontal)

power_K5_h = c(list_power_h$"pstar_K5", list_power_h$"psigma_K5")
eb_K5_h = c(list_eb_h$"pstar_K5", list_eb_h$"psigma_K5")
df_K5_h = data.frame(del = rep(vec_del, 2), 
                     power = power_K5_h, 
                     case = factor(rep(c("pstar", "psigma"), 
                                       each = ld)), 
                     leb = power_K5_h - eb_K5_h, 
                     ueb = power_K5_h + eb_K5_h)
df_K5_h$case = factor(df_K5_h$case, 
                      levels = c("pstar", "psigma"))
plot_K5_h = ggplot(df_K5_h, 
                   aes(del, power, color = case)) + 
  geom_errorbar(aes(ymin = leb, ymax = ueb, color = case), 
                width = 0.3, linewidth = 1) + 
  geom_line(aes(color = case), linewidth = 1, alpha = 0.5) + 
  geom_point(aes(color = case), size = 2) + 
  labs(title = "K = 5", 
       x = expression(paste(delta)), 
       y = "Empirical Power", 
       color = "Test") +
  ylim(0, 1) + 
  scale_color_manual(values = c("pstar" = "#ED2939",
                                "psigma" = "black"), 
                     labels = expression(
                       "pstar" = paste(p, "*"), 
                       "psigma" = p[sigma])) + 
  plot_details + 
  theme(legend.position = "none") 

##------------------------------------------------------------------------------
## Creates plot for K = 10 (horizontal)

power_K10_h = c(list_power_h$"pstar_K10", list_power_h$"psigma_K10")
eb_K10_h = c(list_eb_h$"pstar_K10", list_eb_h$"psigma_K10")
df_K10_h = data.frame(del = rep(vec_del, 2), 
                      power = power_K10_h, 
                      case = factor(rep(c("pstar", "psigma"), 
                                        each = ld)), 
                      leb = power_K10_h - eb_K10_h, 
                      ueb = power_K10_h + eb_K10_h)
df_K10_h$case = factor(df_K10_h$case, 
                       levels = c("pstar", "psigma"))
plot_K10_h = ggplot(df_K10_h, 
                    aes(del, power, color = case)) + 
  geom_errorbar(aes(ymin = leb, ymax = ueb, color = case), 
                width = 0.3, linewidth = 1) + 
  geom_line(aes(color = case), linewidth = 1, alpha = 0.5) + 
  geom_point(aes(color = case), size = 2) + 
  labs(title = "K = 10", 
       x = expression(paste(delta)), 
       y = "Empirical Power", 
       color = "Test") +
  ylim(0, 1) + 
  scale_color_manual(values = c("pstar" = "#ED2939",
                                "psigma" = "black"), 
                     labels = expression(
                       "pstar" = paste(p, "*"), 
                       "psigma" = p[sigma])) + 
  plot_details + 
  theme(legend.position = "none") 

##------------------------------------------------------------------------------
## Creates plot for K = 3 (Kgon)

power_K3_K = c(list_power_K$"pstar_K3", list_power_K$"psigma_K3")
eb_K3_K = c(list_eb_K$"pstar_K3", list_eb_K$"psigma_K3")
df_K3_K = data.frame(del = rep(vec_del, 2), 
                     power = power_K3_K, 
                     case = factor(rep(c("pstar", "psigma"), 
                                       each = ld)), 
                     leb = power_K3_K - eb_K3_K, 
                     ueb = power_K3_K + eb_K3_K)
df_K3_K$case = factor(df_K3_K$case, 
                      levels = c("pstar", "psigma"))
#---
# Citations: the codes below, as well as the lines followed by #---
# for adding two vertical titles for the y-axis is adapted from 
# https://stackoverflow.com/questions/74534630/ggplot2-two-different-font-
# sizes-on-y-axis 
yaxis_K3_K = list("names" = c("K-gon", "Empirical power"))
yaxis_K3_K = lapply(yaxis_K3_K, function(x) {
  x[[1]] = paste0("<span style='font-size: 40pt'>", x[[1]], "</span>")
  return(x)
})
#---
plot_K3_K = ggplot(df_K3_K, 
                   aes(del, power, color = case)) + 
  geom_errorbar(aes(ymin = leb, ymax = ueb, color = case), 
                width = 0.3, linewidth = 1) + 
  geom_line(aes(color = case), linewidth = 1, alpha = 0.5) + 
  geom_point(aes(color = case), size = 2) + 
  labs(title = NULL, 
       x = expression(paste(delta)), 
       y = paste(yaxis_K3_K[["names"]], collapse = "<br>"), #---
       color = "Test") +
  ylim(0, 1) + 
  scale_color_manual(values = c("pstar" = "#ED2939",
                                "psigma" = "black"), 
                     labels = expression(
                       "pstar" = paste(p, "*"), 
                       "psigma" = p[sigma])) + 
  plot_details + 
  theme(legend.position = "none") + 
  theme(axis.title.y = ggtext::element_markdown()) #---

##------------------------------------------------------------------------------
## Creates plot for K = 5 (Kgon)

power_K5_K = c(list_power_K$"pstar_K5", list_power_K$"psigma_K5")
eb_K5_K = c(list_eb_K$"pstar_K5", list_eb_K$"psigma_K5")
df_K5_K = data.frame(del = rep(vec_del, 2), 
                     power = power_K5_K, 
                     case = factor(rep(c("pstar", "psigma"), 
                                       each = ld)), 
                     leb = power_K5_K - eb_K5_K, 
                     ueb = power_K5_K + eb_K5_K)
df_K5_K$case = factor(df_K5_K$case, 
                      levels = c("pstar", "psigma"))
plot_K5_K = ggplot(df_K5_K, 
                   aes(del, power, color = case)) + 
  geom_errorbar(aes(ymin = leb, ymax = ueb, color = case), 
                width = 0.3, linewidth = 1) + 
  geom_line(aes(color = case), linewidth = 1, alpha = 0.5) + 
  geom_point(aes(color = case), size = 2) + 
  labs(title = NULL, 
       x = expression(paste(delta)), 
       y = "Empirical Power", 
       color = "Test") +
  ylim(0, 1) + 
  scale_color_manual(values = c("pstar" = "#ED2939",
                                "psigma" = "black"), 
                     labels = expression(
                       "pstar" = paste(p, "*"), 
                       "psigma" = p[sigma])) + 
  plot_details + 
  theme(legend.position = "none") 

##------------------------------------------------------------------------------
## Creates plot for K = 10 (Kgon)

power_K10_K = c(list_power_K$"pstar_K10", list_power_K$"psigma_K10")
eb_K10_K = c(list_eb_K$"pstar_K10", list_eb_K$"psigma_K10")
df_K10_K = data.frame(del = rep(vec_del, 2), 
                      power = power_K10_K, 
                      case = factor(rep(c("pstar", "psigma"), 
                                        each = ld)), 
                      leb = power_K10_K - eb_K10_K, 
                      ueb = power_K10_K + eb_K10_K)
df_K10_K$case = factor(df_K10_K$case, 
                       levels = c("pstar", "psigma"))
plot_K10_K = ggplot(df_K10_K, 
                    aes(del, power, color = case)) + 
  geom_errorbar(aes(ymin = leb, ymax = ueb, color = case), 
                width = 0.3, linewidth = 1) + 
  geom_line(aes(color = case), linewidth = 1, alpha = 0.5) + 
  geom_point(aes(color = case), size = 2) + 
  labs(title = NULL, 
       x = expression(paste(delta)), 
       y = "Empirical Power", 
       color = "Test") +
  ylim(0, 1) + 
  scale_color_manual(values = c("pstar" = "#ED2939",
                                "psigma" = "black"), 
                     labels = expression(
                       "pstar" = paste(p, "*"), 
                       "psigma" = p[sigma])) + 
  plot_details + 
  theme(legend.position = "none") 

##------------------------------------------------------------------------------
## Combines and saves plots

plot_comb = ((plot_K3_h | plot_K5_h | plot_K10_h) /
  (plot_K3_K | plot_K5_K | plot_K10_K) &
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

to_save = paste0("Figures/figure_12.png")
if (dir.exists("Figures") == 0) {
  dir.create("Figures")
  ggsave(plot_comb, filename = to_save, width = 21, height = 11)
} else {
  ggsave(plot_comb, filename = to_save, width = 21, height = 11)
}

##------------------------------------------------------------------------------