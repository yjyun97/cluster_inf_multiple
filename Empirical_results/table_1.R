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
num_trial = 100

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
## Computes p-values for the data consistent with the null hypothesis

# estimates variances
ss_hat_sample_null = fun_ss_hat_sample(X_null)
ss_hat_med_null = fun_ss_hat_med(X_null)
# specifies V = V1
set_V = fun_set_V_pre(K, 1)
# creates vectors for storing p-values for each test
#vec_pstar_null = rep(0, num_trial)
vec_psample_null = rep(0, num_trial)
vec_pmed_null = rep(0, num_trial)
vec_psampleBon_null = rep(0, num_trial)
vec_pmedBon_null = rep(0, num_trial)
for (i in 1:num_trial) {
  # runs K-means clustering
  res_Kmeans = KmeansInference::kmeans_estimation(X_null, K, seed = i) 
  cl = res_Kmeans$cluster[[res_Kmeans$iter]] 
  # checks if the clustering outcome gives K clusters
  if (length(unique(cl)) < K) {
    #vec_pstar_null[i] = NA
    vec_psample_null[i] = NA
    vec_pmed_null[i] = NA
    vec_psampleBon_null[i] = NA
    vec_pmedBon_null[i] = NA
  } else {
    #vec_pstar_null[i] =  fun_p_star(X_null, res_Kmeans) 
    vec_psample_null[i] = fun_p_sigma(X_null, res_Kmeans, ss_hat_sample_null) 
    vec_pmed_null[i] = fun_p_sigma(X_null, res_Kmeans, ss_hat_med_null) 
    vec_psampleBon_null[i] = fun_p_sigma_Bon(X_null, res_Kmeans, 
                                             set_V, ss_hat_sample_null)
    vec_pmedBon_null[i] = fun_p_sigma_Bon(X_null, res_Kmeans, 
                                         set_V, ss_hat_med_null)
  }
}

#avg_pstar_null = mean(vec_pstar_null[is.na(vec_pstar_null) == 0])
avg_psample_null = mean(vec_psample_null[is.na(vec_psample_null) == 0])
avg_pmed_null = mean(vec_pmed_null[is.na(vec_pmed_null) == 0])
avg_psampleBon_null = mean(vec_psampleBon_null[is.na(vec_psampleBon_null) == 0])
avg_pmedBon_null = mean(vec_pmedBon_null[is.na(vec_pmedBon_null) == 0])

#print(c(avg_pstar_null, avg_psample_null, avg_pmed_null, 
#        avg_psampleBon_null, avg_pmedBon_null))

##------------------------------------------------------------------------------
## Computes p-values for the data consistent with the alternative hypothesis

# estimates variances
ss_hat_sample_alter = fun_ss_hat_sample(X_alter)
ss_hat_med_alter = fun_ss_hat_med(X_alter)
# specifies V = V1
set_V = fun_set_V_pre(K, 1)
# creates vectors for storing p-values for each test
#vec_pstar_alter = rep(0, num_trial)
vec_psample_alter = rep(0, num_trial)
vec_pmed_alter = rep(0, num_trial)
vec_psampleBon_alter = rep(0, num_trial)
vec_pmedBon_alter = rep(0, num_trial)
for (i in 1:num_trial) {
  # runs K-means clustering
  res_Kmeans = KmeansInference::kmeans_estimation(X_alter, K, seed = i) 
  cl = res_Kmeans$cluster[[res_Kmeans$iter]] 
  # checks if the clustering outcome gives K clusters
  if (length(unique(cl)) < K) {
    print("hello")
    #vec_pstar_alter[i] = NA
    vec_psample_alter[i] = NA
    vec_pmed_alter[i] = NA
    vec_psampleBon_alter[i] = NA
    vec_pmedBon_alter[i] = NA
  } else {
    #vec_pstar_alter[i] =  fun_p_star(X_alter, res_Kmeans) 
    vec_psample_alter[i] = fun_p_sigma(X_alter, res_Kmeans, ss_hat_sample_alter) 
    vec_pmed_alter[i] = fun_p_sigma(X_alter, res_Kmeans, ss_hat_med_alter) 
    vec_psampleBon_alter[i] = fun_p_sigma_Bon(X_alter, res_Kmeans, 
                                             set_V, ss_hat_sample_alter)
    vec_pmedBon_alter[i] = fun_p_sigma_Bon(X_alter, res_Kmeans, 
                                          set_V, ss_hat_med_alter)
  }
}

#avg_pstar_alter = mean(vec_pstar_alter[is.na(vec_pstar_alter) == 0])
avg_psample_alter = mean(vec_psample_alter[is.na(vec_psample_alter) == 0])
avg_pmed_alter = mean(vec_pmed_alter[is.na(vec_pmed_alter) == 0])
avg_psampleBon_alter = 
  mean(vec_psampleBon_alter[is.na(vec_psampleBon_alter) == 0])
avg_pmedBon_alter = mean(vec_pmedBon_alter[is.na(vec_pmedBon_alter) == 0])

#print(c(avg_pstar_alter, avg_psample_alter, avg_pmed_alter, 
#        avg_psampleBon_alter, avg_pmedBon_alter))

##------------------------------------------------------------------------------
## Prints results
## Citations: the codes for creating a LATEX table are from the 
##            xtable documentation.

# prints results
# results = xtable(data.frame(c1 = c(avg_pstar_null, avg_pstar_alter),
#                             c2 = c(avg_psample_null, avg_psample_alter), 
#                             c3 = c(avg_pmed_null, avg_pmed_alter), 
#                             c4 = c(avg_psampleBon_null, avg_psampleBon_alter), 
#                             c5 = c(avg_pmedBon_null, avg_pmedBon_alter),
#                             row.names = c("Null", "Alternative")), 
#                  display = c("s", "g", "g", "g", "g", "g"))
# names(results) = c("pstar", "psample", "pmed", "psampleBon", "pmedBon")
# print(results, math.style.exponents = TRUE)

results = xtable(data.frame(c1 = c(avg_psample_null, avg_psample_alter), 
                            c2 = c(avg_pmed_null, avg_pmed_alter), 
                            c3 = c(avg_psampleBon_null, avg_psampleBon_alter), 
                            c4 = c(avg_pmedBon_null, avg_pmedBon_alter),
                            row.names = c("Null", "Alternative")), 
                 display = c("s", "g", "g", "g", "g"))
names(results) = c("psample", "pmed", "psampleBon", "pmedBon")
print(results, math.style.exponents = TRUE)

##------------------------------------------------------------------------------
## Saves xtable result

to_save = paste0("Tables/table_1.txt")
if (dir.exists("Tables") == 0) {
  dir.create("Tables")
  print(results, file = to_save)
} else {
  print(results, file = to_save)
}

##------------------------------------------------------------------------------
## Saves workspace

to_save = "RData/table_1.RData"
if (dir.exists("RData") == 0) {
  dir.create("RData")
  save.image(file = to_save)
} else {
  save.image(file = to_save)
}

##------------------------------------------------------------------------------