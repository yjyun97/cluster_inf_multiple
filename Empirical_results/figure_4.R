##------------------------------------------------------------------------------
## Imports libraries 

library(ggplot2)
library(patchwork)
source("specifications_for_plots.R")

##------------------------------------------------------------------------------
## Specifies settings 

# settings for data generating process
n = 120
q = 2 
delta = 6
ss = 1
ss_diag = ss * rep(1, q)

# settings for clustering
K = 3

##------------------------------------------------------------------------------
## Sets mu

mu_distinct_h = delta * cbind(seq(0, (K - 1)), rep(0, K))
mu_distinct_K = delta * cbind(cos(2 * pi / K * seq(0, (K - 1))), 
                              sin(2 * pi / K * seq(0, (K - 1))))
tmp = matrix(rep(diag(rep(1, K)), each = round(n / K)), nrow = n)
mu_h = tmp %*% mu_distinct_h
mu_K = tmp %*% mu_distinct_K
true_cl = rep(1:K, each = round(n / K))

##------------------------------------------------------------------------------
## Generates a data set from each distribution 

seed = 2
set.seed(seed)
X_h = t(apply(mu_h, 1, function(x) { MASS::mvrnorm(1, x, diag(ss_diag)) }))
set.seed(seed)
X_K = t(apply(mu_K, 1, function(x) { MASS::mvrnorm(1, x, diag(ss_diag)) }))

##------------------------------------------------------------------------------
## Runs K-means clustering each data set

res_Kmeans_h = KmeansInference::kmeans_estimation(X_h, K, seed = seed) 
res_Kmeans_K = KmeansInference::kmeans_estimation(X_K, K, seed = seed) 
cl_h = res_Kmeans_h$cluster[[res_Kmeans_h$iter]]
cl_K = res_Kmeans_K$cluster[[res_Kmeans_K$iter]]

##------------------------------------------------------------------------------
## Visualizes the clustering outcomes on the two data sets

df_h = data.frame(x = X_h[, 1], 
                  y = X_h[, 2],
                  col_vec = as.factor(cl_h),
                  shape_vec = as.factor(true_cl))
plot_h = ggplot(df_h, aes(x = x, y = y)) + 
  geom_point(aes(shape = shape_vec, color = col_vec), size = 5) + 
  scale_shape_manual(values = c("1" = 15, "2" = 19, "3" = 17)) + 
  scale_color_manual(values = c("1" = "#E27396", 
                                "2" = "#9C84B3", 
                                "3" = "#7AC1C1")) + 
  labs(title = "Horizontal", 
       color = "Clusters Formed", shape = "True Clusters", 
       x = NULL, y = NULL) +
  plot_details + 
  theme(legend.position = "none")

df_K = data.frame(x = X_K[, 1], 
                  y = X_K[, 2],
                  col_vec = as.factor(cl_K),
                  shape_vec = as.factor(true_cl))
plot_K = ggplot(df_K, aes(x = x, y = y)) + 
  geom_point(aes(shape = shape_vec, color = col_vec), size = 5) + 
  scale_shape_manual(values = c("1" = 15, "2" = 19, "3" = 17)) + 
  scale_color_manual(values = c("1" = "#E27396", 
                                "2" = "#9C84B3", 
                                "3" = "#7AC1C1")) + 
  labs(title = "K-gon", 
       color = "Clusters Formed", shape = "True Clusters", 
       x = NULL, y = NULL) +
  plot_details + 
  theme(legend.position = "none")

##------------------------------------------------------------------------------
## Combines and saves plots

plot_comb = (plot_h | plot_K) &
  theme(plot.title = element_text(size = 40))

to_save = paste0("Figures/figure_4.png")
if (dir.exists("Figures") == 0) {
  dir.create("Figures")
  ggsave(plot_comb, filename = to_save, width = 15, height = 7)
} else {
  ggsave(plot_comb, filename = to_save, width = 15, height = 7)
}

##------------------------------------------------------------------------------