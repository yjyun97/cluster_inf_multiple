##------------------------------------------------------------------------------
## Reads in functions 

source("../Proposed_methods/functions_pval.R")
source("../Proposed_methods/functions_trunc.R")
source("../Proposed_methods/functions_helper.R")
source("../Proposed_methods/functions_settings.R")

##------------------------------------------------------------------------------
## Specifies settings 

# settings for the null hypothesis
choice_of_V = "all" # ("all", "farthest", "closest")

# settings for the choice of test
which_test = "psigmaBon" # ("psigma", "psigmaBon", "psigmaJ", "pstar", "pstarJ")
if (which_test == "psigmaBon") {
  which_set_V = 3 # (1, 2, 3)
}
if (choice_of_V %in% c("farthest", "closest")) {  
  g = 1 
}
which_ss_hat = "sigma" # ("sigma", "med", "sample")

# settings for data generating process
n = 120
q = 5 
which_mu = "horizontal" # c("horizontal", "Kgon")
delta = 0 # c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6)
ss = 1
ss_diag = ss * rep(1, q)

# settings for clustering
K = 10 

# settings for simulations
num_trial = 1500

##------------------------------------------------------------------------------
## Sets file name 

if (which_test %in% c("psigma", "psigmaBon", "psigmaJ") & 
    which_ss_hat == "sigma") {
  known_or_unknown = "known"
} else {
  known_or_unknown = "unknown"
}
if (choice_of_V == "all" | which_test == "psigmaBon") {
  pre_or_dep = "pre"
} else {
  pre_or_dep = "dep"
}
if (which_test == "psigma") {
  name_test = paste0("p", which_ss_hat)
} else if (which_test == "psigmaBon") {
  name_test = paste0("p", which_ss_hat, "BonV", which_set_V)
} else {
  name_test = which_test 
}
if (pre_or_dep == "pre") {
  if (name_test == "psigma" & q == 20) {
    file_name = paste(known_or_unknown, pre_or_dep, name_test, paste0("q", q), 
                      paste0("K", K), which_mu, paste0("delta", delta), 
                      sep = "_")
  } else {
    file_name = paste(known_or_unknown, pre_or_dep, name_test, paste0("K", K), 
                      which_mu, paste0("delta", delta), sep = "_")
  }
} else {
  if (which_test == "psigmaJ" & q == 20) {
    file_name = paste(known_or_unknown, pre_or_dep, name_test, choice_of_V, 
                      paste0("q", q), paste0("K", K), paste0("g", g), which_mu, 
                      paste0("delta", delta), sep = "_")
  } else {
    file_name = paste(known_or_unknown, pre_or_dep, name_test, choice_of_V, 
                      paste0("K", K), paste0("g", g), which_mu, 
                      paste0("delta", delta), sep = "_")
  }
}

##------------------------------------------------------------------------------
## Sets mu

if (which_mu == "horizontal") {
  mu_distinct = delta * cbind(seq(0, (K - 1)), rep(0, K))
} else if (which_mu == "Kgon") {
  mu_distinct = delta * cbind(cos(2 * pi / K * seq(0, (K - 1))), 
                              sin(2 * pi / K * seq(0, (K - 1))))
} else {
  stop("The specified which_mu does not exist.")
}
tmp = matrix(rep(diag(rep(1, K)), each = round(n / K)), nrow = n)
mu = tmp %*% mu_distinct
mu = cbind(mu, matrix(rep(0, n * (q - 2)), nrow = n))

##------------------------------------------------------------------------------
## Computes p-values

vec_p = rep(0, num_trial)
for (i in 1:num_trial) {
  set.seed(i)
  # generates data 
  X = t(apply(mu, 1, function(x) { MASS::mvrnorm(1, x, diag(ss_diag)) }))
  # runs K-means clustering
  res_Kmeans = KmeansInference::kmeans_estimation(X, K, seed = i) 
  cl = res_Kmeans$cluster[[res_Kmeans$iter]] 
  # checks if the clustering outcome gives K clusters
  if (length(unique(cl)) < K) {
    vec_p[i] = NA
  } else {
    # estimates sigma^2 if needed
    if (which_test %in% c("pstar", "pstarJ")) {
    } else {
      if (which_ss_hat == "sigma") {
        ss_hat = ss
      } else if (which_ss_hat == "med") {
        ss_hat = fun_ss_hat_med(X)
      } else if (which_ss_hat == "sample") {
        ss_hat = fun_ss_hat_sample(X)
      } else {
        stop("The specified which_ss does not exist.")
      }
    }
    # specifies mathcal{V} 
    if (choice_of_V == "all" & which_test != "psigmaBon") {
      set_V = fun_set_V_pre(K, 2)
    } else if (which_test == "psigmaBon") {
      set_V = fun_set_V_pre(K, which_set_V)
    } else if (choice_of_V == "farthest" | choice_of_V == "closest"){
      set_V = (fun_set_V_dep(cl, choice_of_V, g, X))[[1]]
    } else {
      stop("The specified choice_of_V does not exist.")
    }
    # checks if the clustering outcome is consistent with the alternative
    mat_V = (fun_mat_V(cl, "pre", set_V))[[1]]
    diff_means = sum(apply(mat_V, 2, function(x) { return(t(mu) %*% x) }) ** 2)
    if (diff_means < 1e-30 & delta != 0) {
      vec_p[i] = NA
      next
    } 
    # runs the tests 
    if (which_test == "psigma") {
      vec_p[i] = fun_p_sigma(X, res_Kmeans, ss_hat, choice_of_V, set_V) 
    } else if (which_test == "psigmaJ") {
      vec_p[i] = fun_p_sigma_J(X, res_Kmeans, ss_hat, choice_of_V, set_V)
    } else if (which_test == "pstar") {
      vec_p[i] = fun_p_star(X, res_Kmeans, choice_of_V, set_V) 
    } else if (which_test == "pstarJ") {
      vec_p[i] = fun_p_star_J(X, res_Kmeans, choice_of_V, set_V)
    } else if (which_test == "psigmaBon") {
      vec_p[i] = fun_p_sigma_Bon(X, res_Kmeans, set_V, ss_hat)
    } else {
      stop("The specified which_test does not exist.")
    }
  }
}

##------------------------------------------------------------------------------
## Saves the workspace 

to_save = paste0("RData_", file_name, ".RData")
if (dir.exists("RData") == 0) {
  dir.create("RData")
  save.image(file = to_save)
} else {
  save.image(file = to_save)
}

##------------------------------------------------------------------------------