##------------------------------------------------------------------------------
## Functions for specifying set_V 

fun_set_V_dep = function(cl, choice_of_V, g, X) {
  # Input(s):
  # - cl: vector of cluster assignments
  # - choice_of_V: process by which mathcal{V} is chosen 
  # - g: size of mathcal{V}
  # - X: data matrix of dimensions n by q
  # Output(s):
  # - set_V: 2 by g matrix whose columns correspond to the elements of   
  #   mathcal{V} chosen according to Setting 1 or Setting 2 of Section 4
  # - set_V_c: 2 by g_c matrix whose columns correspond to the elements of 
  #   mathcal{V}_all \ mathcal{V} chosen according to Setting 1 or Setting 2 of 
  #   Section 4
  n = nrow(X)
  if (g > choose(K, 2)) {
    stop("g can be at most K choose 2.")
  }
  g_all = choose(K, 2)
  g_c = g_all - g
  mat_V_all = matrix(0, n, g_all)
  store_pairs = matrix(0, 2, g_all)
  count = 1
  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
      mat_V_all[, count] = fun_v(cl, i, j)
      store_pairs[, count] = c(i, j)
      count = count + 1
    }
  }
  vec_dist_all = rep(0, g_all)
  for (i in 1:g_all) {
    tmp = t(X) %*% mat_V_all[, i]
    vec_dist_all[i] = sum(tmp ** 2)
  }
  if (choice_of_V == "closest") {
    vec_ind = sort(vec_dist_all, index.return = TRUE)$ix
  } else if (choice_of_V == "farthest") {
    vec_ind = sort(vec_dist_all, decreasing = TRUE, index.return = TRUE)$ix
  } else {
    stop("The specified choice_of_V is not compatible with this function.")
  } 
  set_V = matrix(store_pairs[, vec_ind[1:g]], 2, g)
  set_V_c = matrix(store_pairs[, vec_ind[(g + 1):g_all]], 2, g_c)
  to_return = list(set_V, set_V_c)
  return(to_return)
}

fun_set_V_pre = function(K, which_set_V = 1) {
  # Input(s):
  # - K: number of clusters produced by K-means clustering
  # - which_set_V: i in mathcal{V}_i defined in Section 6 
  # Output(s):
  # - set_V: matrix whose columns correspond to the elements of mathcal{V}_i 
  if (which_set_V == 1) {
    set_V = matrix(0, 2, choose(K, 2))
    count = 1
    for (i in 1:(K - 1)) {
      for (j in (i + 1):K) {
        set_V[1, count] = i
        set_V[2, count] = j
        count = count + 1
      }
    }
  } else if (which_set_V == 2) {
    set_V = rbind(1:(K - 1), 2:K)
  } else if (which_set_V == 3) {
    set_V = rbind(rep(1, K - 1), 2:K)
  } else {
    stop("The specified which_set_V is not compatible with this function.")
  }
  return(set_V)
}

##------------------------------------------------------------------------------
## Functions for computing estimates of sigma^2 

fun_ss_hat_med <- function(X){
  # Input(s): 
  # - X: data matrix of dimensions n by q 
  # Output(s): 
  # - estimate of the parameter sigma^2 proposed by Chen et al. [2023]
  X_tmp = (X - rep(1, nrow(X)) %*% matrix(apply(X, 2, median), nrow = 1)) ** 2
  to_return = (1 / qchisq(1 / 2, 1)) * median(X_tmp) 
  return(to_return)
}

fun_ss_hat_sample = function(X) {
  # Citations: 
  # - adapted from Yun and Barber[2023]'s function fun_ss_hat_all from  
  #   https://github.com/yjyun97/cluster_inf_unknown_var/blob/main/functions.R
  # Input(s): 
  # - X: data matrix of dimensions n by q 
  # Output(s): 
  # - estimate of the parameter sigma^2 proposed by Gao et al. [2022]
  n = dim(X)[1]
  q = dim(X)[2]
  center = colMeans(X)
  mat_center = matrix(rep(center, n), n, q, byrow = TRUE) 
  to_return = sum((X - mat_center) ** 2) / ((n - 1) * q)
  return(to_return)
}

##------------------------------------------------------------------------------