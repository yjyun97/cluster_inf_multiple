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
g = 2

##------------------------------------------------------------------------------
## Subsets the penguins dataset to bill depth and bill length of female penguins 

penguins = na.omit(penguins)
penguins = penguins[penguins$sex == "female", ]
X_orig = as.matrix(penguins[, c(3, 4)])
true_y_orig = as.vector(as.matrix(penguins[, 1]))

##------------------------------------------------------------------------------
## Subsets further to the Adelie and Gentoo species and standardizes

ind_orig = which(true_y_orig %in% c("Adelie", "Gentoo"))
X_orig = X_orig[ind_orig, ]
X_orig = apply(X_orig, 2, function(x){(x - mean(x)) / sd(x)})
true_y_orig = true_y_orig[ind_orig]

##------------------------------------------------------------------------------
## Computes p-values

# estimates variances
ss_hat_sample = fun_ss_hat_sample(X_orig)
ss_hat_med = fun_ss_hat_med(X_orig)
# runs K-means clustering
res_Kmeans = KmeansInference::kmeans_estimation(X_orig, K, seed = 1) 
cl = res_Kmeans$cluster[[res_Kmeans$iter]] 
if (length(unique(cl)) < K) {
  print("The number of clusters defined is less than K")
} else {
  # computes p-values
  set_V_f = (fun_set_V_dep(cl, "farthest", g, X_orig))[[1]]
  set_V_c = (fun_set_V_dep(cl, "closest", g, X_orig))[[1]]
  psample_f = fun_p_sigma(X_orig, res_Kmeans, ss_hat_sample, "farthest", 
                          set_V_f) 
  psample_c = fun_p_sigma(X_orig, res_Kmeans, ss_hat_sample, "closest", 
                          set_V_c) 
  pmed_f = fun_p_sigma(X_orig, res_Kmeans, ss_hat_med, "farthest", 
                          set_V_f) 
  pmed_c = fun_p_sigma(X_orig, res_Kmeans, ss_hat_med, "closest", 
                          set_V_c) 
  psampleJ_f =  fun_p_sigma_J(X_orig, res_Kmeans, ss_hat_sample, "farthest", 
                              set_V_f)
  psampleJ_c = fun_p_sigma_J(X_orig, res_Kmeans, ss_hat_sample, "closest", 
                             set_V_c)
  pmedJ_f =  fun_p_sigma_J(X_orig, res_Kmeans, ss_hat_med, "farthest", 
                              set_V_f)
  pmedJ_c = fun_p_sigma_J(X_orig, res_Kmeans, ss_hat_med, "closest", 
                             set_V_c)
}

##------------------------------------------------------------------------------
## Prints results
## Citations: the codes for creating a LATEX table are from the 
##            xtable documentation.

# prints results
results = xtable(data.frame(c1 = c(psample_f, psample_c),
                            c2 = c(pmed_f, pmed_c),
                            c3 = c(psampleJ_f, psampleJ_c),
                            c4 = c(pmedJ_f, pmedJ_c),
                            row.names = c("Setting 1", "Setting 2")), 
                 display = c("s", "g", "g", "g", "g"))
names(results) = c("psample", "pmed", "psampleJ", "pmedJ")
print(results, math.style.exponents = TRUE)

##------------------------------------------------------------------------------
## Saves xtable result

to_save = paste0("Tables/table_2.txt")
if (dir.exists("Tables") == 0) {
  dir.create("Tables")
  print(results, file = to_save)
} else {
  print(results, file = to_save)
}

##------------------------------------------------------------------------------
## Saves workspace

to_save = "RData/table_2.RData"
if (dir.exists("RData") == 0) {
  dir.create("RData")
  save.image(file = to_save)
} else {
  save.image(file = to_save)
}

##------------------------------------------------------------------------------