##------------------------------------------------------------------------------
## Imports libraries 

library(ggplot2)
library(patchwork)
library(xtable)
library(palmerpenguins)
source("specifications_for_plots.R")
source("../Proposed_methods/functions_helper.R")
source("../Proposed_methods/functions_pval.R")
source("../Proposed_methods/functions_trunc.R")
source("../Proposed_methods/functions_settings.R")

##------------------------------------------------------------------------------
## Specifies settings

K = 4

##------------------------------------------------------------------------------
## Subsets the penguins dataset to bill depth and bill length of female penguins 

penguins = na.omit(penguins)
penguins = penguins[penguins$sex == "female", ]
X_orig = as.matrix(penguins[, c(3, 4)])
true_y_orig = as.vector(as.matrix(penguins[, 1]))

##------------------------------------------------------------------------------
## Subsets further to the Adelie species and standardizes

ind_null = which(true_y_orig == "Adelie")
X_null = X_orig[ind_null, ]
X_null = apply(X_null, 2, function(x){(x - mean(x)) / sd(x)})
true_y_null = true_y_orig[ind_null]

##------------------------------------------------------------------------------
## Subsets further to the Adelie and Gentoo species and standardizes

ind_alter = which(true_y_orig %in% c("Adelie", "Gentoo"))
X_alter = X_orig[ind_alter, ]
X_alter = apply(X_alter, 2, function(x){(x - mean(x)) / sd(x)})
true_y_alter = true_y_orig[ind_alter]

##------------------------------------------------------------------------------
## Visualizes the data set consistent with the null hypothesis

res_Kmeans_null = KmeansInference::kmeans_estimation(X_null, K, seed = 1) 
cl_null = res_Kmeans_null$cluster[[res_Kmeans_null$iter]] 
df_null = data.frame(x = X_null[, 1],
                     y = X_null[, 2],
                     col_vec = as.factor(cl_null),
                     shape_vec = as.factor(true_y_null))
plot_null = ggplot(df_null, aes(x = x, y = y)) + 
  geom_point(aes(color = col_vec), size = 5) + 
  scale_color_manual(values = c("1" = "#E27396", 
                                "2" = "#9C84B3", 
                                "3" = "#7AC1C1", 
                                "4" = "#FADB75")) + 
  labs(title = "Null", 
       color = "Clusters Formed", shape = "True Clusters", 
       x = "Bill Length (mm)", 
       y = "Bill Depth (mm)") + 
  plot_details + 
  theme(legend.position = "none") 

##------------------------------------------------------------------------------
## Visualizes the data set consistent with the alternative hypothesis

res_Kmeans_alter = KmeansInference::kmeans_estimation(X_alter, K, seed = 1) 
cl_alter = res_Kmeans_alter$cluster[[res_Kmeans_alter$iter]] 
df_alter = data.frame(x = X_alter[, 1],
                     y = X_alter[, 2],
                     col_vec = as.factor(cl_alter),
                     shape_vec = as.factor(true_y_alter))
plot_alter = ggplot(df_alter, aes(x = x, y = y)) + 
  geom_point(aes(shape = shape_vec, color = col_vec), size = 5) + 
  scale_shape_manual(values = c("Adelie" = 16, 
                                "Gentoo" = 15)) + 
  scale_color_manual(values = c("1" = "#E27396", 
                                "2" = "#9C84B3", 
                                "3" = "#7AC1C1", 
                                "4" = "#FADB75")) + 
  labs(title = "Alternative", 
       color = "Clusters Formed", shape = "True Clusters", 
       x = "Bill Length (mm)", 
       y = "Bill Depth (mm)") + 
  plot_details + 
  theme(legend.position = "none") 

##------------------------------------------------------------------------------
## Combines and saves plots

plot_comb = ((plot_null | plot_alter) &
  theme(legend.position = "right",
        legend.key.size = unit(2, "cm"), 
        legend.text.align = 0, 
        plot.title = element_text(hjust = 0.5, size = 50), 
        axis.title = element_text(size = 40, color = "black")) &
  #---
  # Citations: the line below for changing the size of points in the legend has 
  # been modified from https://stackoverflow.com/questions/20415963/how-to-
  # increase-the-size-of-points-in-legend-of-ggplot2
  guides(color = guide_legend(override.aes = list(size = 7)), 
         shape = guide_legend(override.aes = list(size = 7)))) +
  #---
  plot_layout(guides = "collect")

to_save = paste0("Figures/figure_14.png")
if (dir.exists("Figures") == 0) {
  dir.create("Figures")
  ggsave(plot_comb, filename = to_save, width = 21, height = 8)
} else {
  ggsave(plot_comb, filename = to_save, width = 21, height = 8)
}

##------------------------------------------------------------------------------