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

# horizontal
list_h = list()
for (which_test_ in c("psigmaJ", "psigma")) {
  for (choice_of_V_ in c("farthest", "closest")) {
    mat_p_h = matrix(0, num_sim, ld)
    for (i_ in 1:ld) {
      load(paste0("RData/known_dep_", which_test_, "_", choice_of_V_, "_K20_g3", 
                  "_horizontal_delta", vec_del[i_], ".RData"))
      mat_p_h[, i_] = vec_p
    }
    list_h[[paste0(which_test_, "_", choice_of_V)]] = mat_p_h
  }
}

# Kgon
list_K = list()
for (which_test_ in c("psigmaJ", "psigma")) {
  for (choice_of_V_ in c("farthest", "closest")) {
    mat_p_K = matrix(0, num_sim, ld)
    for (i_ in 1:ld) {
      load(paste0("RData/known_dep_", which_test_, "_", choice_of_V_, "_K20_g3", 
                  "_Kgon_delta", vec_del[i_], ".RData"))
      mat_p_K[, i_] = vec_p
    }
    list_K[[paste0(which_test_, "_", choice_of_V_)]] = mat_p_K
  }
}

##------------------------------------------------------------------------------
## Computes empirical power

# horizontal
list_power_h = list()
list_eb_h = list()
for (which_test_ in c("psigmaJ", "psigma")) {
  for (choice_of_V_ in c("farthest", "closest")) {
    res_h = apply(list_h[[paste0(which_test_, "_", choice_of_V_)]], 2, 
                  function(x) { fun_power(x, thresh) })
    list_power_h[[paste0(which_test_, "_", choice_of_V_)]] = res_h[1, ]
    list_eb_h[[paste0(which_test_, "_", choice_of_V_)]] = res_h[2, ]
  }
}

# Kgon
list_power_K = list()
list_eb_K = list()
for (which_test_ in c("psigmaJ", "psigma")) {
  for (choice_of_V_ in c("farthest", "closest")) {
    res_K = apply(list_K[[paste0(which_test_, "_", choice_of_V_)]], 2,
                  function(x) { fun_power(x, thresh) })
    list_power_K[[paste0(which_test_, "_", choice_of_V_)]] = res_K[1, ]
    list_eb_K[[paste0(which_test_, "_", choice_of_V_)]] = res_K[2, ]
  }
}

##------------------------------------------------------------------------------
## Creates plot for farthest (horizontal)

power_f_h = c(list_power_h$"psigmaJ_farthest", list_power_h$"psigma_farthest")
eb_f_h = c(list_eb_h$"psigmaJ_farthest", list_eb_h$"psigma_farthest")
df_f_h = data.frame(del = rep(vec_del, 2), power = power_f_h, 
                    case = factor(rep(c("psigmaJ", "psigma"), each = ld)), 
                    leb = power_f_h - eb_f_h, ueb = power_f_h + eb_f_h)
df_f_h$case = factor(df_f_h$case, levels = c("psigmaJ", "psigma"))
#---
# Citations: the codes below, as well as the lines followed by #---
# for adding two vertical titles for the y-axis is adapted from 
# https://stackoverflow.com/questions/74534630/ggplot2-two-different-font-
# sizes-on-y-axis 
yaxis_f_h = list("names" = c("Horizontal", "Empirical power"))
yaxis_f_h = lapply(yaxis_f_h, function(x) {
  x[[1]] = paste0("<span style='font-size: 45pt'>", x[[1]], "</span>")
  return(x)
})
#---
plot_f_h = ggplot(df_f_h, aes(del, power, color = case)) + 
  geom_errorbar(aes(ymin = leb, ymax = ueb, color = case), 
                width = 0.3, linewidth = 1) + 
  geom_line(aes(color = case), linewidth = 1) + 
  geom_point(aes(color = case), size = 2) + 
  labs(title = "Setting 1", 
       x = expression(paste(delta)), 
       y = paste(yaxis_f_h[["names"]], collapse = "<br>"), #---
       color = "Test") +
  scale_color_manual(values = c("psigmaJ" = "#ED2939",
                                "psigma" = "black"), 
                     labels = expression(
                       "psigmaJ" = p[paste(sigma, ",J")], 
                       "psigma" = p[sigma])) + 
  plot_details + 
  theme(legend.position = "none") + 
  theme(axis.title.y = ggtext::element_markdown()) #---

##------------------------------------------------------------------------------
## Creates plot for closest (horizontal)

power_c_h = c(list_power_h$"psigmaJ_closest", list_power_h$"psigma_closest")
eb_c_h = c(list_eb_h$"psigmaJ_closest", list_eb_h$"psigma_closest")
df_c_h = data.frame(del = rep(vec_del, 2), power = power_c_h, 
                    case = factor(rep(c("psigmaJ", "psigma"), each = ld)), 
                    leb = power_c_h - eb_c_h, ueb = power_c_h + eb_c_h)
df_c_h$case = factor(df_c_h$case, levels = c("psigmaJ", "psigma"))
plot_c_h = ggplot(df_c_h, aes(del, power, color = case)) + 
  geom_errorbar(aes(ymin = leb, ymax = ueb, color = case), 
                width = 0.3, linewidth = 1) + 
  geom_line(aes(color = case), linewidth = 1) + 
  geom_point(aes(color = case), size = 2) + 
  labs(title = "Setting 2", 
       x = expression(paste(delta)), 
       y = "Empirical Power",
       color = "Test") +
  scale_color_manual(values = c("psigmaJ" = "#ED2939",
                                "psigma" = "black"), 
                     labels = expression(
                       "psigmaJ" = p[paste(sigma, ",J")], 
                       "psigma" = p[sigma])) + 
  plot_details + 
  theme(legend.position = "none") + 
  theme(axis.title.y = ggtext::element_markdown()) #---

##------------------------------------------------------------------------------
## Creates plot for farthest (Kgon)

power_f_K = c(list_power_K$"psigmaJ_farthest", list_power_K$"psigma_farthest")
eb_f_K = c(list_eb_K$"psigmaJ_farthest", list_eb_K$"psigma_farthest")
df_f_K = data.frame(del = rep(vec_del, 2), power = power_f_K, 
                    case = factor(rep(c("psigmaJ", "psigma"), each = ld)), 
                    leb = power_f_K - eb_f_K, ueb = power_f_K + eb_f_K)
df_f_K$case = factor(df_f_K$case, levels = c("psigmaJ", "psigma"))
#---
# Citations: the codes below, as well as the lines followed by #---
# for adding two vertical titles for the y-axis is adapted from 
# https://stackoverflow.com/questions/74534630/ggplot2-two-different-font-
# sizes-on-y-axis 
yaxis_f_K = list("names" = c("K-gon", "Empirical power"))
yaxis_f_K = lapply(yaxis_f_K, function(x) {
  x[[1]] = paste0("<span style='font-size: 45pt'>", x[[1]], "</span>")
  return(x)
})
#---
plot_f_K = ggplot(df_f_K, aes(del, power, color = case)) + 
  geom_errorbar(aes(ymin = leb, ymax = ueb, color = case), 
                width = 0.3, linewidth = 1) + 
  geom_line(aes(color = case), linewidth = 1) + 
  geom_point(aes(color = case), size = 2) + 
  labs(title = NULL, 
       x = expression(paste(delta)), 
       y = paste(yaxis_f_K[["names"]], collapse = "<br>"), #---
       color = "Test") +
  scale_color_manual(values = c("psigmaJ" = "#ED2939",
                                "psigma" = "black"), 
                     labels = expression(
                       "psigmaJ" = p[paste(sigma, ",J")], 
                       "psigma" = p[sigma])) + 
  plot_details + 
  theme(legend.position = "none") + 
  theme(axis.title.y = ggtext::element_markdown()) #---

##------------------------------------------------------------------------------
## Creates plot for closest (Kgon)

power_c_K = c(list_power_K$"psigmaJ_closest", list_power_K$"psigma_closest")
eb_c_K = c(list_eb_K$"psigmaJ_closest", list_eb_K$"psigma_closest")
df_c_K = data.frame(del = rep(vec_del, 2), power = power_c_K, 
                    case = factor(rep(c("psigmaJ", "psigma"), each = ld)), 
                    leb = power_c_K - eb_c_K, ueb = power_c_K + eb_c_K)
df_c_K$case = factor(df_c_K$case, levels = c("psigmaJ", "psigma"))
plot_c_K = ggplot(df_c_K, aes(del, power, color = case)) + 
  geom_errorbar(aes(ymin = leb, ymax = ueb, color = case), 
                width = 0.3, linewidth = 1) + 
  geom_line(aes(color = case), linewidth = 1) + 
  geom_point(aes(color = case), size = 2) + 
  labs(title = NULL, 
       x = expression(paste(delta)), 
       y = "Empirical Power", 
       color = "Test") +
  scale_color_manual(values = c("psigmaJ" = "#ED2939",
                                "psigma" = "black"), 
                     labels = expression(
                       "psigmaJ" = p[paste(sigma, ",J")], 
                       "psigma" = p[sigma])) + 
  plot_details + 
  theme(legend.position = "none") 

##------------------------------------------------------------------------------
## Combines and saves plots

plot_comb = ((plot_f_h | plot_c_h) / 
  (plot_f_K | plot_c_K)  &
  theme(legend.position = "right", 
        legend.key.size = unit(2, "cm"), 
        legend.text.align = 0) &
  #---
  # Citations: the line below for changing the size of points in the legend has 
  # been taken from https://stackoverflow.com/questions/20415963/how-to-
  # increase-the-size-of-points-in-legend-of-ggplot2
  guides(color = guide_legend(override.aes = list(size = 7)))) +
  #---
  plot_layout(guides = "collect")

to_save = paste0("Figures/figure_8.png")
if (dir.exists("Figures") == 0) {
  dir.create("Figures")
  ggsave(plot_comb, filename = to_save, width = 18, height = 13)
} else {
  ggsave(plot_comb, filename = to_save, width = 18, height = 13)
}

##------------------------------------------------------------------------------