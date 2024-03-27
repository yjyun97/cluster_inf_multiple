##------------------------------------------------------------------------------
## Function for computing p_sigma

fun_p_sigma = function(X, res_Kmeans, ss, choice_of_V = "all", set_V = NULL) {
  # Citation(s):
  # - Uses the function TChisqRatioApprox from the package KmeansInference
  # - Uses the functions Intervals_full and interval_intersection from the 
  #   package intervals 
  # - Applies Chen and Witten [2023]'s way of addressing numerical instability 
  #   due to finite precision in the chunk of codes enclosed by a pair of #---s, 
  #   where the codes have been modified from 
  #   https://github.com/yiqunchen/KmeansInference/blob/main/R/kmeans_inf.R
  # Input(s):
  # - X: data matrix of dimensions n by q
  # - res_Kmeans: list containing various quantities outputted by K-means 
  #   clustering, specifically the output of the function kmeans_estimation 
  #   of the package KmeansInference
  # - ss: sigma^2 or its estimate
  # - choice_of_V: process by which mathcal{V} is chosen 
  # - set_V: matrix whose columns correspond to the elements of mathcal{V}
  # Output(s):
  # - p_sigma: p-value p_sigma of Section 3.2
  if (choice_of_V == "all" | choice_of_V == "pre") {
  } else if (choice_of_V %in% c("farthest", "closest")) {
    choice_of_V = "pre"
  } else {
    stop("The specified choice_of_V does not exist.")
  }
  if (!is.null(set_V)) {
    if (sum(set_V[1, ] - set_V[2, ] > 0) > 0) {
      stop("The second entry has to be larger than the first one for each 
           element of set_V.")
    }
  }
  # saves quantities for later use
  n = dim(X)[1]
  q = dim(X)[2]
  cl_final = res_Kmeans$cluster[[res_Kmeans$iter]]
  K = length(cl_final)
  if (choice_of_V == "pre" & is.null(set_V)) {
    stop("set_V has to be given for the specified choice_of_V.")
  } 
  if (choice_of_V == "pre") {
    if (ncol(set_V) == choose(K, 2)) {
      choice_of_V = "all"
    }
  }
  res_PE = fun_PE(cl_final, choice_of_V, set_V)
  PE = res_PE[[1]] 
  dim_E = res_PE[[2]]
  d = q * dim_E
  # computes T_sigma 
  T_sigma = sqrt(sum((PE %*% X) ** 2) / ss)
  # computes S_sigma 
  S_sigma = fun_S_sigma(X, PE, res_Kmeans, ss)
  #---
  # applies Chen and Witten [2023]'s way of addressing numerical instability 
  # due to finite precision. The codes have been modified from 
  # https://github.com/yiqunchen/KmeansInference/blob/main/R/kmeans_inf.R
  S_sigma = 
    intervals::interval_union(S_sigma,
                              intervals::Intervals_full(c(T_sigma - (1e-09), 
                                                          T_sigma + (1e-09))),
                              check_valid = FALSE)
  #---
  # computes p_sigma
  set_denom = S_sigma ** 2 
  set_numer_tmp = intervals::Intervals_full(c(T_sigma ** 2, Inf))
  set_numer = intervals::interval_intersection(set_numer_tmp, set_denom)
  p_sigma = KmeansInference:::TChisqRatioApprox(d, set_numer, set_denom) 
  return(p_sigma)
}

##------------------------------------------------------------------------------
## Function for computing p_{sigma, J}

fun_p_sigma_J = function(X, res_Kmeans, ss, choice_of_V = "all", set_V = NULL) {
  # Citation(s):
  # - Uses the function TChisqRatioApprox from the package KmeansInference
  # - Uses the functions Intervals_full and interval_intersection from the 
  #   package intervals 
  # - Applies Chen and Witten [2023]'s way of addressing numerical instability 
  #   due to finite precision in the chunk of codes enclosed by a pair of #---s, 
  #   where the codes have been modified from 
  #   https://github.com/yiqunchen/KmeansInference/blob/main/R/kmeans_inf.R
  # Input(s):
  # - X: data matrix of dimensions n by q
  # - res_Kmeans: list containing various quantities outputted by K-means 
  #   clustering, specifically the output of the function kmeans_estimation 
  #   of the package KmeansInference
  # - ss: sigma^2 or its estimate
  # - choice_of_V: process by which mathcal{V} is chosen 
  # - set_V: matrix whose columns correspond to the elements of mathcal{V}
  # Output(s):
  # - p_sigma_J: p-value p_sigma_J of Section 4
  if (!(choice_of_V %in% c("pre", "all", "farthest", "closest"))) {
    stop("The specified choice_of_V does not exist.")
  }
  if (choice_of_V == "pre" | choice_of_V == "all") {
    p_sigma_J = fun_p_sigma(X, res_Kmeans, ss, choice_of_V, set_V)
    return(p_sigma_J)
  } 
  if (!is.null(set_V)) {
    if (sum(set_V[1, ] - set_V[2, ] > 0) > 0) {
      stop("The second entry has to be larger than the first one for each 
           element of set_V.")
    }
  }
  # saves quantities for later use
  n = dim(X)[1]
  q = dim(X)[2]
  cl_final = res_Kmeans$cluster[[res_Kmeans$iter]]
  K = length(cl_final)
  if (ncol(set_V) == choose(K, 2)) {
    choice_of_V = "all"
    p_sigma_J = fun_p_sigma(X, res_Kmeans, ss, choice_of_V, set_V)
    return(p_sigma_J)
  }
  if (is.null(set_V)) {
    stop("set_V has to be given for the specified choice_of_V.")
  }
  res_PE = fun_PE(cl_final, choice_of_V, set_V)
  PE = res_PE[[1]] 
  dim_E = res_PE[[2]]
  d = q * dim_E
  # computes T_sigma 
  T_sigma = sqrt(sum((PE %*% X) ** 2) / ss)
  # computes S_sigma_J 
  g = ncol(set_V)
  S_sigma_J = fun_S_sigma_J(X, PE, res_Kmeans, choice_of_V, g, ss)
  #---
  # applies Chen and Witten [2023]'s way of addressing numerical instability 
  # due to finite precision. The codes have been modified from 
  # https://github.com/yiqunchen/KmeansInference/blob/main/R/kmeans_inf.R
  S_sigma_J = 
    intervals::interval_union(S_sigma_J,
                              intervals::Intervals_full(c(T_sigma - (1e-09), 
                                                          T_sigma + (1e-09))),
                              check_valid = FALSE)
  #---
  # computes p_sigma_J
  set_denom = S_sigma_J ** 2 
  set_numer_tmp = intervals::Intervals_full(c(T_sigma ** 2, Inf))
  set_numer = intervals::interval_intersection(set_numer_tmp, set_denom)
  p_sigma_J = KmeansInference:::TChisqRatioApprox(d, set_numer, set_denom) 
  return(p_sigma_J)
}

##------------------------------------------------------------------------------
## Function for computing p^star

fun_p_star = function(X, res_Kmeans, choice_of_V = "all", set_V = NULL) {
  # Citation(s):
  # - Uses the function TChisqRatioApprox from the package KmeansInference
  # - Uses the functions Intervals_full and interval_intersection from the 
  #   package intervals 
  # - Uses the function fun_F_to_chi2 from 
  #   https://github.com/yjyun97/cluster_inf_unknown_var/blob/main/functions.R
  # - Applies Chen and Witten [2023]'s way of addressing numerical instability 
  #   due to finite precision in the chunk of codes enclosed by a pair of #---s, 
  #   where the codes have been modified from 
  #   https://github.com/yiqunchen/KmeansInference/blob/main/R/kmeans_inf.R
  # Input(s):
  # - X: data matrix of dimensions n by q
  # - res_Kmeans: list containing various quantities outputted by K-means 
  #   clustering, specifically the output of the function kmeans_estimation 
  #   of the package KmeansInference
  # - choice_of_V: process by which mathcal{V} is chosen 
  # - set_V: matrix whose columns correspond to the elements of mathcal{V}
  # Output(s):
  # - p_star: p-value p_star of Section 5.2.1
  if (choice_of_V == "all" | choice_of_V == "pre") {
  } else if (choice_of_V %in% c("farthest", "closest")) {
    choice_of_V = "pre"
  } else {
    stop("The specified choice_of_V does not exist.")
  }
  if (!is.null(set_V)) {
    if (sum(set_V[1, ] - set_V[2, ] > 0) > 0) {
      stop("The second entry has to be larger than the first one for each 
           element of set_V.")
    }
  }
  # saves quantities for later use
  n = dim(X)[1]
  q = dim(X)[2]
  cl_final = res_Kmeans$cluster[[res_Kmeans$iter]]
  K = length(cl_final)
  if (choice_of_V == "pre" & is.null(set_V)) {
    stop("set_V has to be given for the specified choice_of_V.")
  } 
  if (choice_of_V == "pre") {
    if (ncol(set_V) == choose(K, 2)) {
      choice_of_V = "all"
    }
  }
  res_PE = fun_PE(cl_final, choice_of_V, set_V)
  PE = res_PE[[1]] 
  dim_E = res_PE[[2]]
  set_K = fun_set_K(cl_final, choice_of_V, set_V)
  P1 = fun_P1(cl_final, set_K)
  P2 = fun_P2(PE, P1)
  d = q * dim_E
  d_star = q * (length(which(cl_final %in% set_K)) - length(set_K))
  r_star = d_star / d
  # computes T_star
  T_star_numer = sum((PE %*% X) ** 2) / d
  T_star_denom = sum((P1 %*% X) ** 2) / d_star
  T_star = T_star_numer / T_star_denom
  # computes S_star 
  S_star = fun_S_star(X, PE, P1, P2, r_star, res_Kmeans)
  #---
  # applies Chen and Witten [2023]'s way of addressing numerical instability 
  # due to finite precision. The codes have been modified from 
  # https://github.com/yiqunchen/KmeansInference/blob/main/R/kmeans_inf.R
  S_star = 
    intervals::interval_union(S_star,
                              intervals::Intervals_full(c(T_star - (1e-09), 
                                                          T_star + (1e-09))),
                              check_valid = FALSE)
  #---
  # computes p_star
  set_denom = S_star
  set_numer_tmp = intervals::Intervals_full(c(T_star, Inf))
  set_numer = intervals::interval_intersection(set_numer_tmp, set_denom)
  set_denom_chi2 = fun_F_to_chi2(set_denom, d, d_star)
  set_numer_chi2 = fun_F_to_chi2(set_numer, d, d_star)
  p_star = KmeansInference:::TChisqRatioApprox(d, set_numer_chi2, 
                                               set_denom_chi2)
  return(p_star)
}

##------------------------------------------------------------------------------
## Function for computing p^star_J

fun_p_star_J = function(X, res_Kmeans, choice_of_V = "all", set_V = NULL) {
  # Citation(s):
  # - Uses the function TChisqRatioApprox from the package KmeansInference
  # - Uses the functions Intervals_full and interval_intersection from the 
  #   package intervals 
  # - Uses the function fun_F_to_chi2 from 
  #   https://github.com/yjyun97/cluster_inf_unknown_var/blob/main/functions.R
  # - Applies Chen and Witten [2023]'s way of addressing numerical instability 
  #   due to finite precision in the chunk of codes enclosed by a pair of #---s, 
  #   where the codes have been modified from 
  #   https://github.com/yiqunchen/KmeansInference/blob/main/R/kmeans_inf.R
  # Input(s):
  # - X: data matrix of dimensions n by q
  # - res_Kmeans: list containing various quantities outputted by K-means 
  #   clustering, specifically the output of the function kmeans_estimation 
  #   of the package KmeansInference
  # - choice_of_V: process by which mathcal{V} is chosen 
  # - set_V: matrix whose columns correspond to the elements of mathcal{V}
  # Output(s):
  # - p_star_J: p-value p_star_J of Section 5.2.2
  if (!(choice_of_V %in% c("pre", "all", "farthest", "closest"))) {
    stop("The specified choice_of_V does not exist.")
  }
  if (choice_of_V == "pre" | choice_of_V == "all") {
    p_star_J = fun_p_star(X, res_Kmeans, choice_of_V, set_V)
    return(p_star_J)
  } 
  if (!is.null(set_V)) {
    if (sum(set_V[1, ] - set_V[2, ] > 0) > 0) {
      stop("The second entry has to be larger than the first one for each 
           element of set_V.")
    }
  }
  # saves quantities for later use
  n = dim(X)[1]
  q = dim(X)[2]
  cl_final = res_Kmeans$cluster[[res_Kmeans$iter]]
  K = length(cl_final)
  if (ncol(set_V) == choose(K, 2)) {
    choice_of_V = "all"
    p_star_J = fun_p_star(X, res_Kmeans, choice_of_V, set_V)
    return(p_star_J)
  }
  if (is.null(set_V)) {
    stop("set_V has to be given for the specified choice_of_V.")
  }
  res_PE = fun_PE(cl_final, choice_of_V, set_V)
  PE = res_PE[[1]] 
  dim_E = res_PE[[2]]
  set_K = fun_set_K(cl_final, choice_of_V, set_V)
  P1 = fun_P1(cl_final, set_K)
  P2 = fun_P2(PE, P1)
  d = q * dim_E
  d_star = q * (length(which(cl_final %in% set_K)) - length(set_K))
  r_star = d_star / d
  # computes T_star
  T_star_numer = sum((PE %*% X) ** 2) / d
  T_star_denom = sum((P1 %*% X) ** 2) / d_star
  T_star = T_star_numer / T_star_denom
  # computes S_star_J 
  g = ncol(set_V)
  S_star_J = fun_S_star_J(X, PE, P1, P2, r_star, res_Kmeans, choice_of_V, g)
  #---
  # applies Chen and Witten [2023]'s way of addressing numerical instability 
  # due to finite precision. The codes have been modified from 
  # https://github.com/yiqunchen/KmeansInference/blob/main/R/kmeans_inf.R
  S_star_J = 
    intervals::interval_union(S_star_J,
                              intervals::Intervals_full(c(T_star - (1e-09), 
                                                          T_star + (1e-09))),
                              check_valid = FALSE)
  #---
  # computes p_star_J
  set_denom = S_star_J
  set_numer_tmp = intervals::Intervals_full(c(T_star, Inf))
  set_numer = intervals::interval_intersection(set_numer_tmp, set_denom)
  set_denom_chi2 = fun_F_to_chi2(set_denom, d, d_star)
  set_numer_chi2 = fun_F_to_chi2(set_numer, d, d_star)
  p_star_J = KmeansInference:::TChisqRatioApprox(d, set_numer_chi2, 
                                                 set_denom_chi2)
  return(p_star_J)
}

##------------------------------------------------------------------------------
## Function for computing p_{sigma,Bon}

fun_p_sigma_Bon = function(X, res_Kmeans, set_V, ss) {
  # Input(s):
  # - X: data matrix of dimensions n by q
  # - res_Kmeans: list containing various quantities outputted by K-means 
  #   clustering, specifically the output of the function kmeans_estimation 
  #   of the package KmeansInference
  # - set_V: matrix whose columns correspond to the elements of mathcal{V}
  # - ss: sigma^2 or its estimate
  # Output(s):
  # - p_sigma_Bon: p-value p_{sigma,Bon} of Section 3.1
  size_set_V = ncol(set_V)
  vec_p = rep(0, size_set_V)
  for (i in 1:size_set_V) {
    vec_p[i] = fun_p_sigma(X, res_Kmeans, ss, "pre", matrix(set_V[, i], 2, 1))
  }
  vec_p_adj = vec_p * size_set_V
  p_sigma_Bon = min(min(vec_p_adj), 1)
  return(p_sigma_Bon)
}

##------------------------------------------------------------------------------