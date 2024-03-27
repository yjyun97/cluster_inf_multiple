##------------------------------------------------------------------------------
## Imports libraries

library(ggplot2)
library(patchwork)
source("specifications_for_plots.R")

##------------------------------------------------------------------------------
## Imports values 

load("RData/known_pre_psigmaBonV1_q2_K3_delta0.RData")
K3 = vec_p
load("RData/known_pre_psigmaBonV1_q2_K5_delta0.RData")
K5 = vec_p
load("RData/known_pre_psigmaBonV1_q2_K7_delta0.RData")
K7 = vec_p

##------------------------------------------------------------------------------
## Creates plot 

df = data.frame(p = c(K3, K5, K7), 
                case = as.factor(rep(c("3", "5", "7"), each = num_trial)))
df$case = factor(df$case, levels = c("3", "5", "7"))

plot_ = ggplot(df, aes(sample = p, color = case)) + 
  stat_qq(distribution = qunif, size = 1) + 
  geom_abline(intercept = 0, slope = 1, color = "black", 
              linewidth = 2, linetype = "dashed") + 
  labs(title = expression("QQ Plot for" ~ p[paste(sigma, ",Bon")]), 
       x = "Theoretical Quantiles", 
       y = "Empirical Quantiles",
       color = "K") + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  scale_color_manual(values = c("3" = "#ef476f",
                                "5" = "#ffd166", 
                                "7" = "#26547c")) + 
  plot_details + 
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  theme(axis.title = element_text(size = 30, 
                                  color = "black"),
        plot.title = element_text(hjust = 0.5, size = 35)) 

##------------------------------------------------------------------------------
## Saves plot

to_save = paste0("Figures/figure_1.png")
if (dir.exists("Figures") == 0) {
  dir.create("Figures")
  ggsave(plot_, filename = to_save, width = 9, height = 7)
} else {
  ggsave(plot_, filename = to_save, width = 9, height = 7)
}

##------------------------------------------------------------------------------