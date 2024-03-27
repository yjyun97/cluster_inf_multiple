##------------------------------------------------------------------------------
## Function for computing S_sigma

fun_S_sigma = function(X, PE, res_Kmeans, ss) {
  # Citations: 
  # - adapted from the function kmeans_compute_S_iso of Chen and Witten [2023] 
  #   from
  #   https://github.com/yiqunchen/KmeansInference/blob/main/R/trunc_sets.R
  # - Uses the function solve_one_ineq_complement from the package 
  #   KmeansInference
  # - Uses the functions reduce, Intervals_full, interval_complement, and 
  #   interval_intersection from the package intervals 
  # Input(s):
  # - X: data matrix of dimensions n by q
  # - PE: projection matrix P_{mathcal{E}}
  # - res_Kmeans: list containing various quantities outputted by K-means 
  #   clustering, specifically the output of the function kmeans_estimation 
  #   of the package KmeansInference
  # - ss: sigma^2 or its estimate
  # Output(s):
  # - S_sigma: truncation set S_sigma of Section 3.2 
  # saves quantities for later use
  n = nrow(X)
  cl_center_all = res_Kmeans$centers
  cl_all = do.call(rbind, res_Kmeans$cluster)
  J_plus_1 = nrow(cl_all)
  list_intervals_c = list() 
  # computes S_i_l_0
  cl_center_init = res_Kmeans$random_init_obs
  cl_init = cl_all[1,]
  K = length(cl_center_init)
  count = 1
  for (i in 1:n){
    cl_center_init_Xi = cl_init[i]
    coef_left = fun_coef_S_sigma_init(X, PE, i, 
                                      cl_center_init[cl_center_init_Xi], ss) 
    for (l in 1:K){
      coef_right = fun_coef_S_sigma_init(X, PE, i, cl_center_init[l], ss)
      coef = fun_coef_diff(coef_left, coef_right)
      interval_c = 
        KmeansInference:::solve_one_ineq_complement(coef$lam_1, coef$lam_2, 
                                                    coef$lam_3)
      list_intervals_c[[count]] = interval_c
      count = count + 1
    }
  }
  store_size = length(list_intervals_c)
  count = 1
  # computes S_i_l_j
  if (J_plus_1 > 1){
    for (j in 1:(J_plus_1 - 1)) {
      cl_j = cl_all[j + 1, ]
      cl_j_minus_1 = cl_all[j, ]
      w_j_minus_1 = matrix(0, nrow = K, ncol = nrow(X)) 
      for (l in c(1:K)){
        denom_l = sum(cl_j_minus_1 == l) 
        numer_l = rep(0, n)
        numer_l[cl_j_minus_1 == l] = 1
        w_j_minus_1[l, ] = numer_l / denom_l 
      }
      for (i in c(1:n)){
        cl_j_Xi = cl_j[i]
        coef_left = fun_coef_S_sigma(X, PE, i, cl_j_Xi, w_j_minus_1, ss) 
        for (l in 1:K){
          coef_right = fun_coef_S_sigma(X, PE, i, l, w_j_minus_1, ss)
          coef = fun_coef_diff(coef_left, coef_right)
          interval_c = 
            KmeansInference:::solve_one_ineq_complement(coef$lam_1, coef$lam_2, 
                                                        coef$lam_3)
          list_intervals_c[[store_size + count]] = interval_c
          count = count + 1
        }
      }
    }
  }
  S_sigma_c = do.call('c', list_intervals_c)
  S_sigma_c = matrix(S_sigma_c, ncol=2, byrow=T)
  S_sigma_c = intervals::reduce(intervals::Intervals_full(S_sigma_c),
                                check_valid = FALSE)
  S_sigma = intervals::interval_complement(S_sigma_c)
  int_pos = intervals::Intervals_full(c(0, Inf)) 
  S_sigma = intervals::interval_intersection(S_sigma, int_pos)
  return(S_sigma)
}

##------------------------------------------------------------------------------
## Function for computing S_{sigma, V}

fun_S_sigma_V = function(X, PE, cl_final, choice_of_V, g, ss) {
  # Citations: 
  # - adapted from the function kmeans_compute_S_iso of Chen and Witten [2023] 
  #   from
  #   https://github.com/yiqunchen/KmeansInference/blob/main/R/trunc_sets.R
  # - Uses the function solve_one_ineq_complement from the package 
  #   KmeansInference
  # - Uses the functions reduce, Intervals_full, interval_complement, and 
  #   interval_intersection from the package intervals 
  # Input(s):
  # - X: data matrix of dimensions n by q
  # - PE: projection matrix P_{mathcal{E}}
  # - cl_final: vector of cluster assignments that are outputted by the final 
  #   iteration of K-means clustering
  # - choice_of_V: process by which mathcal{V} is chosen 
  # - g: size of mathcal{V}
  # - ss: sigma^2 or its estimate
  # Output(s):
  # - S_sigma_V: truncation set S_sigma_V of Section 4
  # saves quantities for later use
  n = length(cl_final)
  K = max(cl_final)
  res_set_V = fun_set_V_dep(cl_final, choice_of_V, g, X)
  res_mat_V = fun_mat_V(cl, choice_of_V, res_set_V[[1]], res_set_V[[2]])
  mat_V = res_mat_V[[1]]
  mat_V_c = res_mat_V[[2]]
  g_c = ncol(mat_V_c)
  list_intervals_c = list()
  # computes S_sigma_V
  count = 1 
  for (i in 1:g) {
    for (j in 1:g_c) {
      v = mat_V[, i]
      vc = mat_V_c[, j]
      coef_v = fun_coef_S_sigma_V(X, PE, v, ss)
      coef_vc = fun_coef_S_sigma_V(X, PE, vc, ss)
      if (choice_of_V == "farthest") {
        coef = fun_coef_diff(coef_vc, coef_v)
      } else if (choice_of_V == "closest") {
        coef = fun_coef_diff(coef_v, coef_vc)
      } else {
        stop("The specified choice_of_V is not compatible with this function.")
      }
      interval_c = 
        KmeansInference:::solve_one_ineq_complement(coef$lam_1, coef$lam_2,
                                                    coef$lam_3)
      list_intervals_c[[count]] = interval_c
      count = count + 1
    }
  } 
  S_sigma_V_c = do.call('c', list_intervals_c)
  S_sigma_V_c = matrix(S_sigma_V_c, ncol = 2, byrow = TRUE)
  S_sigma_V_c = 
    intervals::reduce(intervals::Intervals_full(S_sigma_V_c), 
                      check_valid = FALSE)
  S_sigma_V = intervals::interval_complement(S_sigma_V_c)
  int_pos = intervals::Intervals_full(c(0, Inf)) 
  S_sigma_V = intervals::interval_intersection(S_sigma_V, int_pos)
  return(S_sigma_V)
}

##------------------------------------------------------------------------------
## Function for computing S_{sigma, J}

fun_S_sigma_J = function(X, PE, res_Kmeans, choice_of_V, g, ss) {
  # Citation(s):
  # - Uses the function interval_intersection from the package intervals  
  # Input(s):
  # - X: data matrix of dimensions n by q
  # - PE: projection matrix P_{mathcal{E}}
  # - res_Kmeans: list containing various quantities outputted by K-means 
  #   clustering, specifically the output of the function kmeans_estimation 
  #   of the package KmeansInference
  # - choice_of_V: process by which mathcal{V} is chosen 
  # - g: size of mathcal{V}
  # - ss: sigma^2 or its estimate
  # Output(s):
  # - S_sigma_J: truncation set S_sigma_J of Section 4
  cl_final = res_Kmeans$cluster[[res_Kmeans$iter]]
  S_sigma = fun_S_sigma(X, PE, res_Kmeans, ss)
  S_sigma_V = fun_S_sigma_V(X, PE, cl_final, choice_of_V, g, ss)
  S_sigma_J = intervals::interval_intersection(S_sigma, S_sigma_V) 
  return(S_sigma_J)
}

##------------------------------------------------------------------------------
## Function for computing S^star

fun_S_star = function(X, PE, P1, P2, r_star, res_Kmeans) {
  # Citations: 
  # - adapted from the function kmeans_compute_S_iso of Chen and Witten [2023] 
  #   from
  #   https://github.com/yiqunchen/KmeansInference/blob/main/R/trunc_sets.R
  # - Uses the functions reduce, Intervals_full, interval_complement, and 
  #   interval_intersection from the package intervals 
  # Input(s):
  # - X: data matrix of dimensions n by q
  # - PE: projection matrix P_{mathcal{E}}
  # - P1: projection matrix P_1
  # - P2: projection matrix P_2
  # - r_star: r^star defined in Section 5.2
  # - res_Kmeans: list containing various quantities outputted by K-means 
  #   clustering, specifically the output of the function kmeans_estimation 
  #   of the package KmeansInference
  # Output(s):
  # - S_star: truncation set S^star of Section 5.2.1
  # saves quantities for later use
  n = nrow(X)
  cl_center_all = res_Kmeans$centers
  cl_all = do.call(rbind, res_Kmeans$cluster)
  J_plus_1 = nrow(cl_all)
  list_intervals_c = list() 
  # computes S_i_l_0
  cl_center_init = res_Kmeans$random_init_obs
  cl_init = cl_all[1,]
  K = length(cl_center_init)
  count = 1
  for (i in 1:n){
    cl_center_init_Xi = cl_init[i]
    coef_left = fun_coef_S_star_init(X, PE, P1, P2, r_star, i, 
                                     cl_center_init[cl_center_init_Xi])
    for (l in 1:K){
      coef_right = fun_coef_S_star_init(X, PE, P1, P2, r_star, i, 
                                        cl_center_init[l])
      coef = fun_coef_diff(coef_left, coef_right, num_terms = 5)
      list_root = fun_roots_sqrt_S_star(coef, r_star)
      interval_c = fun_ineq_sqrt_S_star(list_root)
      list_intervals_c[[count]] = interval_c
      count = count + 1
    }
  }
  store_size = length(list_intervals_c)
  count = 1
  # computes S_i_l_j
  if (J_plus_1 > 1){
    for (j in 1:(J_plus_1 - 1)) {
      cl_j = cl_all[j + 1, ]
      cl_j_minus_1 = cl_all[j, ]
      w_j_minus_1 = matrix(0, nrow = K, ncol = nrow(X)) 
      for (l in c(1:K)){
        denom_l = sum(cl_j_minus_1 == l) 
        numer_l = rep(0, n)
        numer_l[cl_j_minus_1 == l] = 1
        w_j_minus_1[l, ] = numer_l / denom_l 
      }
      for (i in c(1:n)){
        cl_j_Xi = cl_j[i]
        coef_left = fun_coef_S_star(X, PE, P1, P2, r_star, i, cl_j_Xi, w_j_minus_1)
        for (l in 1:K){
          coef_right = fun_coef_S_star(X, PE, P1, P2, r_star, i, l, w_j_minus_1)
          coef = fun_coef_diff(coef_left, coef_right, num_terms = 5)
          list_root = fun_roots_sqrt_S_star(coef, r_star)
          interval_c = fun_ineq_sqrt_S_star(list_root)
          list_intervals_c[[store_size + count]] = interval_c
          count = count + 1
        }
      }
    }
  }
  sqrt_S_star_c = do.call('c', list_intervals_c)
  sqrt_S_star_c = matrix(sqrt_S_star_c, ncol=2, byrow=T)
  sqrt_S_star_c = intervals::reduce(intervals::Intervals_full(sqrt_S_star_c),
                                    check_valid = FALSE)
  sqrt_S_star = intervals::interval_complement(sqrt_S_star_c)
  int_pos = intervals::Intervals_full(c(0, Inf)) 
  sqrt_S_star = intervals::interval_intersection(sqrt_S_star, int_pos)
  S_star = sqrt_S_star ** 2
  return(S_star)
}

##------------------------------------------------------------------------------
## Function for computing S^star_V

fun_S_star_V = function(X, PE, P1, P2, r_star, cl_final, choice_of_V, g) {
  # Citations: 
  # - adapted from the function kmeans_compute_S_iso of Chen and Witten [2023] 
  #   from
  #   https://github.com/yiqunchen/KmeansInference/blob/main/R/trunc_sets.R
  # - Uses the functions reduce, Intervals_full, interval_complement, and 
  #   interval_intersection from the package intervals 
  # Input(s):
  # - X: data matrix of dimensions n by q
  # - PE: projection matrix P_{mathcal{E}}
  # - P1: projection matrix P_1
  # - P2: projection matrix P_2
  # - r_star: r^star defined in Section 5.2
  # - cl_final: vector of cluster assignments that are outputted by the final 
  #   iteration of K-means clustering
  # - choice_of_V: process by which mathcal{V} is chosen 
  # - g: size of mathcal{V}
  # Output(s):
  # - S_star_V: truncation set S^star_V of Section 5.2.2
  # saves quantities for later use
  n = length(cl_final)
  K = max(cl_final)
  res_set_V = fun_set_V_dep(cl_final, choice_of_V, g, X)
  res_mat_V = fun_mat_V(cl, choice_of_V, res_set_V[[1]], res_set_V[[2]])
  mat_V = res_mat_V[[1]]
  mat_V_c = res_mat_V[[2]]
  g_c = ncol(mat_V_c)
  list_intervals_c = list()
  # computes S_star_V
  count = 1 
  for (i in 1:g) {
    for (j in 1:g_c) {
      v = mat_V[, i]
      vc = mat_V_c[, j]
      coef_v = fun_coef_S_star_V(X, PE, P1, P2, r_star, v)
      coef_vc = fun_coef_S_star_V(X, PE, P1, P2, r_star, vc)
      if (choice_of_V == "farthest") {
        coef = fun_coef_diff(coef_vc, coef_v, num_terms = 5)
      } else if (choice_of_V == "closest") {
        coef = fun_coef_diff(coef_v, coef_vc, num_terms = 5)
      } else {
        stop("The specified choice_of_V is not compatible with this function.")
      }
      list_root = fun_roots_sqrt_S_star(coef, r_star)
      interval_c = fun_ineq_sqrt_S_star(list_root)
      list_intervals_c[[count]] = interval_c
      count = count + 1
    }
  } 
  sqrt_S_star_V_c = do.call('c', list_intervals_c)
  sqrt_S_star_V_c = matrix(sqrt_S_star_V_c, ncol = 2, byrow = TRUE)
  sqrt_S_star_V_c = 
    intervals::reduce(intervals::Intervals_full(sqrt_S_star_V_c), 
                      check_valid = FALSE)
  sqrt_S_star_V = intervals::interval_complement(sqrt_S_star_V_c)
  int_pos = intervals::Intervals_full(c(0, Inf)) 
  sqrt_S_star_V = intervals::interval_intersection(sqrt_S_star_V, int_pos)
  S_star_V = sqrt_S_star_V ** 2
  return(S_star_V)
}

##------------------------------------------------------------------------------
## Function for computing S^star_J

fun_S_star_J = function(X, PE, P1, P2, r_star, res_Kmeans, choice_of_V, g) {
  # Citation(s):
  # - Uses the function interval_intersection from the package intervals  
  # Input(s):
  # - X: data matrix of dimensions n by q
  # - PE: projection matrix P_{mathcal{E}}
  # - P1: projection matrix P_1
  # - P2: projection matrix P_2
  # - r_star: r^star defined in Section 5.2
  # - res_Kmeans: list containing various quantities outputted by K-means 
  #   clustering, specifically the output of the function kmeans_estimation 
  #   of the package KmeansInference
  # - choice_of_V: process by which mathcal{V} is chosen 
  # - g: size of mathcal{V}
  # Output(s):
  # - S_star_J: truncation set S_star_J of Section 5.2.2
  cl_final = res_Kmeans$cluster[[res_Kmeans$iter]]
  S_star = fun_S_star(X, PE, P1, P2, r_star, res_Kmeans) 
  S_star_V = fun_S_star_V(X, PE, P1, P2, r_star, cl_final, choice_of_V, g)
  S_star_J = intervals::interval_intersection(S_star, S_star_V) 
  return(S_star_J)
}

##------------------------------------------------------------------------------