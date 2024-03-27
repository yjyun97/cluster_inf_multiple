##------------------------------------------------------------------------------
## Imports libraries

library(ggplot2)
library(patchwork)
library(ggtext)
source("specifications_for_plots.R")

##------------------------------------------------------------------------------
## Imports values 

# setting 1, p_{sigma, J}
load("RData/known_dep_psigmaJ_farthest_K5_g3_delta0.RData")
s1_psigmaJ_K5 = vec_p
load("RData/known_dep_psigmaJ_farthest_K10_g3_delta0.RData")
s1_psigmaJ_K10 = vec_p
load("RData/known_dep_psigmaJ_farthest_K20_g3_delta0.RData")
s1_psigmaJ_K20 = vec_p
# setting 1, p_{sigma}
load("RData/known_dep_psigma_farthest_K5_g3_delta0.RData")
s1_psigma_K5 = vec_p
load("RData/known_dep_psigma_farthest_K10_g3_delta0.RData")
s1_psigma_K10 = vec_p
load("RData/known_dep_psigma_farthest_K20_g3_delta0.RData")
s1_psigma_K20 = vec_p
# setting 2, p_{sigma, J}
load("RData/known_dep_psigmaJ_closest_K5_g3_delta0.RData")
s2_psigmaJ_K5 = vec_p
load("RData/known_dep_psigmaJ_closest_K10_g3_delta0.RData")
s2_psigmaJ_K10 = vec_p
load("RData/known_dep_psigmaJ_closest_K20_g3_delta0.RData")
s2_psigmaJ_K20 = vec_p
# setting 2, p_{sigma}
load("RData/known_dep_psigma_closest_K5_g3_delta0.RData")
s2_psigma_K5 = vec_p
load("RData/known_dep_psigma_closest_K10_g3_delta0.RData")
s2_psigma_K10 = vec_p
load("RData/known_dep_psigma_closest_K20_g3_delta0.RData")
s2_psigma_K20 = vec_p

##------------------------------------------------------------------------------
## Creates plot setting 1, p_{sigma, J}

df_s1_psigmaJ = data.frame(p = c(s1_psigmaJ_K5, s1_psigmaJ_K10, s1_psigmaJ_K20), 
                           case = as.factor(rep(c("5", "10", "20"), 
                                                each = num_trial)))
df_s1_psigmaJ$case = factor(df_s1_psigmaJ$case, levels = c("5", "10", "20"))
#---
# Citations: the codes below, as well as the lines followed by #---
# for adding two vertical titles for the y-axis is adapted from 
# https://stackoverflow.com/questions/74534630/ggplot2-two-different-font-
# sizes-on-y-axis 
yaxis_s1_psigmaJ = list("names" = c("Setting 1", "Empirical Quantiles"))
yaxis_s1_psigmaJ = lapply(yaxis_s1_psigmaJ, function(x) {
  x[[1]] = paste0("<span style='font-size: 45pt'>", x[[1]], "</span>")
  return(x)
})
#---
plot_s1_psigmaJ = ggplot(df_s1_psigmaJ, aes(sample = p, color = case)) + 
  stat_qq(distribution = qunif, size = 1) + 
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linewidth = 2, linetype = "dashed") + 
  labs(title = expression(p[paste(sigma, ",J")]),
       x = "Theoretical Quantiles", 
       y = paste(yaxis_s1_psigmaJ[["names"]], collapse = "<br>"), #---
       color = "K") + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  scale_color_manual(values = c("5" = "#ef476f",
                                "10" = "#ffd166", 
                                "20" = "#26547c")) + 
  plot_details + 
  theme(legend.position = "none") + 
  theme(axis.title.y = ggtext::element_markdown()) #---

##------------------------------------------------------------------------------
## Creates plot setting 1, p_{sigma}

df_s1_psigma = data.frame(p = c(s1_psigma_K5, s1_psigma_K10, s1_psigma_K20), 
                          case = as.factor(rep(c("5", "10", "20"), 
                                               each = num_trial)))
df_s1_psigma$case = factor(df_s1_psigma$case, levels = c("5", "10", "20"))
plot_s1_psigma = ggplot(df_s1_psigma, aes(sample = p, color = case)) + 
  stat_qq(distribution = qunif, size = 1) + 
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linewidth = 2, linetype = "dashed") + 
  labs(title = expression(p[sigma]),
       x = "Theoretical Quantiles", 
       y = "Empirical Quantiles",
       color = "K") + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  scale_color_manual(values = c("5" = "#ef476f",
                                "10" = "#ffd166", 
                                "20" = "#26547c")) + 
  plot_details + 
  theme(legend.position = "none") 

##------------------------------------------------------------------------------
## Creates plot setting 2, p_{sigma, J}

df_s2_psigmaJ = data.frame(p = c(s2_psigmaJ_K5, s2_psigmaJ_K10, s2_psigmaJ_K20), 
                           case = as.factor(rep(c("5", "10", "20"), 
                                                each = num_trial)))
df_s2_psigmaJ$case = factor(df_s2_psigmaJ$case, levels = c("5", "10", "20"))
#---
# Citations: the codes below, as well as the lines followed by #---
# for adding two vertical titles for the y-axis is adapted from 
# https://stackoverflow.com/questions/74534630/ggplot2-two-different-font-
# sizes-on-y-axis 
yaxis_s2_psigmaJ = list("names" = c("Setting 2", "Empirical Quantiles"))
yaxis_s2_psigmaJ = lapply(yaxis_s2_psigmaJ, function(x) {
  x[[1]] = paste0("<span style='font-size: 45pt'>", x[[1]], "</span>")
  return(x)
})
#---
plot_s2_psigmaJ = ggplot(df_s2_psigmaJ, aes(sample = p, color = case)) + 
  stat_qq(distribution = qunif, size = 1) + 
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linewidth = 2, linetype = "dashed") + 
  labs(title = NULL,
       x = "Theoretical Quantiles", 
       y = paste(yaxis_s2_psigmaJ[["names"]], collapse = "<br>"), #---
       color = "K") + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  scale_color_manual(values = c("5" = "#ef476f",
                                "10" = "#ffd166", 
                                "20" = "#26547c")) + 
  plot_details + 
  theme(legend.position = "none") + 
  theme(axis.title.y = ggtext::element_markdown()) #---

##------------------------------------------------------------------------------
## Creates plot setting 2, p_{sigma}

df_s2_psigma = data.frame(p = c(s2_psigma_K5, s2_psigma_K10, s2_psigma_K20), 
                          case = as.factor(rep(c("5", "10", "20"), 
                                               each = num_trial)))
df_s2_psigma$case = factor(df_s2_psigma$case, levels = c("5", "10", "20"))
plot_s2_psigma = ggplot(df_s2_psigma, aes(sample = p, color = case)) + 
  stat_qq(distribution = qunif, size = 1) + 
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linewidth = 2, linetype = "dashed") + 
  labs(title = NULL,
       x = "Theoretical Quantiles", 
       y = "Empirical Quantiles",
       color = "K") + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  scale_color_manual(values = c("5" = "#ef476f",
                                "10" = "#ffd166", 
                                "20" = "#26547c")) + 
  plot_details + 
  theme(legend.position = "none") 

##------------------------------------------------------------------------------
## Combines and saves plots

plot_comb = ((plot_s1_psigmaJ | plot_s1_psigma) / 
  (plot_s2_psigmaJ | plot_s2_psigma) & 
  theme(legend.position = "right", 
        legend.direction = "vertical", 
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 40), 
        plot.title = element_text(hjust = 0.5, size = 45)) &
  #---
  # Citations: the line below for changing the size of points in the legend has 
  # been taken from https://stackoverflow.com/questions/20415963/how-to-
  # increase-the-size-of-points-in-legend-of-ggplot2
  guides(color = guide_legend(override.aes = list(size = 7)))) +
  #---
  plot_layout(guides = "collect")

to_save = paste0("Figures/figure_7.png")
if (dir.exists("Figures") == 0) {
  dir.create("Figures")
  ggsave(plot_comb, filename = to_save, width = 17, height = 13)
} else {
  ggsave(plot_comb, filename = to_save, width = 17, height = 13)
}

##------------------------------------------------------------------------------