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

K = 10

##------------------------------------------------------------------------------
## Subsets the penguins dataset to bill depth and bill length of female penguins 
## of the Adelie and Gentoo species and standardizes

penguins = na.omit(penguins)
penguins = penguins[penguins$sex == "female", ]
X_orig = as.matrix(penguins[, c(3, 4)])
true_y_orig = as.vector(as.matrix(penguins[, 1]))
# subsets to the Adelie and Gentoo species
ind = which(true_y_orig %in% c("Adelie", "Gentoo"))
X_orig = X_orig[ind, ]
true_y_orig = true_y_orig[ind]
# standardizes
X_orig = apply(X_orig, 2, function(x){(x - mean(x)) / sd(x)})

##------------------------------------------------------------------------------
## Visualizes the first two dimensions of the data set

res_Kmeans_orig = KmeansInference::kmeans_estimation(X_orig, K, seed = 1) 
cl_orig = res_Kmeans_orig$cluster[[res_Kmeans_orig$iter]] 
df_orig = data.frame(x = X_orig[, 1],
                      y = X_orig[, 2],
                      col_vec = as.factor(cl_orig),
                      shape_vec = as.factor(true_y_orig))
plot_orig = ggplot(df_orig, aes(x = x, y = y)) + 
  geom_point(aes(shape = shape_vec, color = col_vec), size = 5) + 
  scale_shape_manual(values = c("Adelie" = 16, 
                                "Gentoo" = 15)) + 
  scale_color_manual(values = c("1" = "#E27396", 
                                "2" = "#9C84B3", 
                                "3" = "#7AC1C1", 
                                "4" = "#FADB75", 
                                "5" = "#6290C3", 
                                "6" = "#273469", 
                                "7" = "#30343F", 
                                "8" = "#69747C", 
                                "9" = "#89023E", 
                                "10" = "#EE8F90")) + 
  labs(title = "Visualization of Penguins Data", 
       color = "Clusters Formed", shape = "True Clusters", 
       x = "Bill Length (mm)", 
       y = "Bill Depth (mm)") + 
  plot_details + 
  theme(plot.title = element_text(hjust = 0.5, size = 50), 
        axis.title = element_text(size = 40, 
                                  color = "black")) + 
  #---
  # Citations: the line below for changing the size of points in the legend has 
  # been modified from https://stackoverflow.com/questions/20415963/how-to-
  # increase-the-size-of-points-in-legend-of-ggplot2
  guides(color = guide_legend(override.aes = list(size = 7)), 
         shape = guide_legend(override.aes = list(size = 7))) 
  #---

##------------------------------------------------------------------------------
## saves plots

to_save = paste0("Figures/figure_15.png")
if (dir.exists("Figures") == 0) {
  dir.create("Figures")
  ggsave(plot_orig, filename = to_save, width = 15, height = 10)
} else {
  ggsave(plot_orig, filename = to_save, width = 15, height = 10)
}

##------------------------------------------------------------------------------