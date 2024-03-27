##------------------------------------------------------------------------------
## Functions for computing projection matrices 

fun_v = function(cl, k1, k2) {
  # Citation(s): 
  # - taken from Yun and Barber[2023]'s function from
  #   https://github.com/yjyun97/cluster_inf_unknown_var/blob/main/functions.R
  # Input(s):
  # - cl: vector of cluster assignments
  # - k1, k2: cluster indices 
  # Output(s):
  # - vector v_{k1,k2} of length n
  n = length(cl)
  tmp1 = (cl == k1) / sum(cl == k1)
  tmp2 = (cl == k2) / sum(cl == k2)
  to_return = tmp1 - tmp2
  return(to_return)
}

fun_mat_V = function(cl, choice_of_V = "all", set_V = NULL, set_V_c = NULL) {
  # Input(s):
  # - cl: vector of cluster assignments
  # - choice_of_V: process by which mathcal{V} is chosen 
  # - set_V: matrix whose columns correspond to the elements of mathcal{V}
  # - set_V_c: matrix whose columns correspond to the elements of 
  #   mathcal{V}_all \ mathcal{V}
  # Output(s):  
  # - list consisting of 
  #   - mat_V: matrix whose columns correspond to the vectors v_{k,k'}s for 
  #     (k,k') in mathcal{V}
  #   - mat_V_c: matrix whose columns correspond to the vectors v_{k,k's}
  #     for (k,k') in mathcal{V_all} \ mathcal{V}; NULL if choice_of_V is set to 
  #     "all" or set_V_c is set to NULL
  n = length(cl)
  K = max(cl)
  if (!(choice_of_V %in% c("pre", "all", "farthest", "closest"))) {
    stop("The specified choice_of_V does not exist.")
  }
  if (choice_of_V != "all" & is.null(set_V)) {
    stop("set_V needs to be specified.")
  }
  if (choice_of_V == "all") {
    set_V = rbind(1:(K - 1), 2:K)
    mat_V_c = NULL
  } else {
    if (!(is.null(set_V_c))) {
      g_c = ncol(set_V_c)
      mat_V_c = matrix(0, n, g_c)
      for (i in 1:g_c) {
        mat_V_c[, i] = fun_v(cl, set_V_c[1, i], set_V_c[2, i])
      }
      mat_V_c = matrix(mat_V_c, n, g_c)
    } else {
      mat_V_c = NULL
    }
  }
  g = ncol(set_V)
  mat_V = matrix(0, n, g)
  for (i in 1:g) {
    mat_V[, i] = fun_v(cl, set_V[1, i], set_V[2, i])
  }
  mat_V = matrix(mat_V, n, g)
  to_return = list(mat_V, mat_V_c)
  return(to_return)
}

fun_PE = function(cl, choice_of_V = "pre", set_V = NULL) {
  # Citation (s):
  # - Uses the function gramSchmidt from the package pracma
  # - Uses the function fast.svd from the package corpcor
  # Input(s):
  # - cl: vector of cluster assignments
  # - choice_of_V: process by which mathcal{V} is chosen 
  # - set_V: matrix whose columns correspond to the elements of mathcal{V}
  # Output(s):
  # - list consisting of 
  #   - PE: projection matrix P_{mathcal{E}}
  #   - dim_E: dimension of mathcal{E}
  res_mat_V = fun_mat_V(cl, choice_of_V, set_V = set_V)
  mat_V = res_mat_V[[1]]
  K = max(cl)
  if (ncol(mat_V) == 1) {
    mat_basis_E = mat_V / sqrt(sum(mat_V ** 2))
    dim_E = 1
  } else {
    if (choice_of_V == "all") {
      mat_basis_E = pracma::gramSchmidt(mat_V)$Q
      dim_E = K - 1
    } else {
      res_svd = corpcor::fast.svd(mat_V)
      mat_basis_E = res_svd$u
      dim_E = length(res_svd$d)
    }
  }
  PE = mat_basis_E %*% t(mat_basis_E)
  to_return = list(PE, dim_E)
  return(to_return)
}

fun_P1 = function(cl, set_K) {
  # Input(s)
  # - cl: vector of cluster assignments
  # - set_K: mathcal{K} defined in Section 5.2
  # Output(s)
  # - P1: projection matrix P_1
  n = length(cl)
  P1 = matrix(0, n, n)
  for (k in set_K) {
    size_Ck = sum(cl == k) 
    indic_k = rep(0, n)
    indic_k[which(cl == k)] = 1
    summand_k = diag(indic_k) - (indic_k %*% t(indic_k)) / size_Ck
    P1 = P1 + summand_k
  }
  return(P1)
}

fun_P2 = function(PE, P1) {
  # Input(s)
  # - PE: projection matrix P_{mathcal{E}}
  # - P1: projection matrix P_1
  # Output(s)
  # - P2: projection matrix P_2
  n = nrow(PE)
  P2 = diag(n) - PE - P1
  return(P2)
}

##------------------------------------------------------------------------------
## Functions for computing the coefficients in inequalities defining the 
## truncation sets 

fun_coef_S_sigma_init = function(X, PE, i, i_prime, ss) {
  # Input(s):
  # - X: data matrix of dimensions n by q 
  # - PE: projection matrix P_{mathcal{E}}
  # - i, i_prime: indices of observations
  # - ss: sigma^2 or its estimate
  # Output(s):
  # - list consisting of the coefficients lambda_{ii',1}, lambda_{ii',2}, and 
  #   lambda_{ii',3} of Proposition 2
  n = nrow(X)
  D = sqrt(ss) * (PE %*% X) / sqrt(sum((PE %*% X) ** 2))
  E = (diag(rep(1, n)) - PE) %*% X
  d = D[i, ] - D[i_prime, ]
  e = E[i, ] - E[i_prime, ]
  lam_1 = sum(d ** 2)
  lam_2 = 2 * t(d) %*% (e)
  lam_3 = sum(e ** 2)
  to_return = list("lam_1" = lam_1, "lam_2" = lam_2, 
                   "lam_3"= lam_3)
  return(to_return)
}

fun_coef_S_sigma = function(X, PE, i, l, mat_w, ss) {
  # Input(s):
  # - X: data matrix of dimensions n by q 
  # - PE: projection matrix P_{mathcal{E}}
  # - i: index of an observation
  # - l: index of a cluster
  # - mat_w: K by n matrix, where the sth entry of lth row corresponds to 
  #   w_{l,s}^(j-1) defined in Section 2 
  # - ss: sigma^2 or its estimate
  # Output(s):
  # - list consisting of the coefficients lambda_{ilj,1}, lambda_{ilj,2}, and 
  #   lambda_{ilj,3} of Proposition 2
  n = nrow(X)
  D = sqrt(ss) * (PE %*% X) / sqrt(sum((PE %*% X) ** 2))
  E = (diag(rep(1, n)) - PE) %*% X
  MD = colSums(diag(mat_w[l, ]) %*% D)
  ME = colSums(diag(mat_w[l, ]) %*% E)
  d = D[i, ] - MD 
  e = E[i, ] - ME
  lam_1 = sum(d ** 2)
  lam_2 = 2 * t(d) %*% (e)
  lam_3 = sum(e ** 2)
  to_return = list("lam_1" = lam_1, "lam_2" = lam_2, 
                   "lam_3"= lam_3)
  return(to_return)
}

fun_coef_S_sigma_V = function(X, PE, v, ss) {
  # Input(s):
  # - X: data matrix of dimensions n by q 
  # - PE: projection matrix P_{mathcal{E}}
  # - v: vector of length n 
  # - ss: sigma^2 or its estimate
  # Output(s):
  # - list consisting of the coefficients lambda_{v,1}, lambda_{v,2}, and 
  #   lambda_{v,3} of Proposition 4
  n = nrow(X)
  D = sqrt(ss) * (PE %*% X) / sqrt(sum((PE %*% X) ** 2))
  E = (diag(rep(1, n)) - PE) %*% X
  d = t(D) %*% v
  e = t(E) %*% v
  lam_1 = sum(d ** 2)
  lam_2 = 2 * t(d) %*% e
  lam_3 = sum(e ** 2)
  to_return = list("lam_1" = lam_1, "lam_2" = lam_2, 
                   "lam_3"= lam_3)
  return(to_return)
}

fun_coef_S_star_init = function(X, PE, P1, P2, r_star, i, i_prime) {
  # Input(s):
  # - X: data matrix of dimensions n by q 
  # - PE: projection matrix P_{mathcal{E}}
  # - P1: projection matrix P_1
  # - P2: projection matrix P_2
  # - r_star: r^star defined in Section 5.2
  # - i, i_prime: indices of observations
  # Output(s):
  # - list consisting of the coefficients lambda_{ii',1}, lambda_{ii',2}, lambda_{ii',3},
  #   lambda_{ii',4}, and lambda_{ii',5} of Proposition 6
  A = (PE %*% X) / sqrt(sum((PE %*% X) ** 2))
  B = (P1 %*% X) / sqrt(sum((P1 %*% X) ** 2))
  C = (P2 %*% X) / sqrt(sum((PE %*% X) ** 2) + sum((P1 %*% X) ** 2))
  a = A[i, ] - A[i_prime, ]
  b = B[i, ] - B[i_prime, ]
  c = C[i, ] - C[i_prime, ]
  lam_1 = sum(a ** 2) + sum(c ** 2)
  lam_2 = 2 * t(a) %*% b * sqrt(r_star)
  lam_3 = 2 * t(a) %*% c 
  lam_4 = 2 * t(b) %*% c * sqrt(r_star)
  lam_5 = (sum(b ** 2) + sum(c ** 2)) * r_star
  to_return = list("lam_1" = lam_1, "lam_2" = lam_2, 
                   "lam_3"= lam_3, "lam_4"= lam_4, "lam_5"= lam_5)
  return(to_return)
}

fun_coef_S_star = function(X, PE, P1, P2, r_star, i, l, mat_w) {
  # Input(s):
  # - X: data matrix of dimensions n by q 
  # - PE: projection matrix P_{mathcal{E}}
  # - P1: projection matrix P_1
  # - P2: projection matrix P_2
  # - r_star: r^star defined in Section 5.2
  # - i: index of an observation
  # - l: index of a cluster
  # - mat_w: K by n matrix, where the sth entry of lth row corresponds to 
  #   w_{l,s}^(j-1) defined in Section 2 
  # Output(s):
  # - list consisting of the coefficients lambda_{ijl,1}, lambda_{ijl,2}, lambda_{ijl,3},
  #   lambda_{ijl,4}, and lambda_{ijl,5} of Proposition 6
  A = (PE %*% X) / sqrt(sum((PE %*% X) ** 2))
  B = (P1 %*% X) / sqrt(sum((P1 %*% X) ** 2))
  C = (P2 %*% X) / sqrt(sum((PE %*% X) ** 2) + sum((P1 %*% X) ** 2))
  MA = colSums(diag(mat_w[l, ]) %*% A)
  MB = colSums(diag(mat_w[l, ]) %*% B)
  MC = colSums(diag(mat_w[l, ]) %*% C)
  a = A[i, ] - MA
  b = B[i, ] - MB
  c = C[i, ] - MC
  lam_1 = sum(a ** 2) + sum(c ** 2)
  lam_2 = 2 * t(a) %*% b * sqrt(r_star)
  lam_3 = 2 * t(a) %*% c 
  lam_4 = 2 * t(b) %*% c * sqrt(r_star)
  lam_5 = (sum(b ** 2) + sum(c ** 2)) * r_star
  to_return = list("lam_1" = lam_1, "lam_2" = lam_2, 
                   "lam_3"= lam_3, "lam_4"= lam_4, "lam_5"= lam_5)
  return(to_return)
}

fun_coef_diff = function(coef_left, coef_right, num_terms = 3) {
  # Citations: 
  # - adapted from the function minus_quad_ineq from 
  #   https://github.com/yiqunchen/KmeansInference/blob/main/R/trunc_sets.R
  # Input(s):
  # - coef_left: list consisting of 3 or 5 elements
  # - coef-right: list consisting of 3 or 5 elements
  # Output(s): 
  # - list consisting of 3 or 5 elements, where ith element is the difference 
  #   between ith elements of coef_left coef_right 
  if (!(num_terms %in% c(3, 5))) {
    stop("The specified num_terms is not compatible with this function.")
  }
  if (num_terms == 3) {
    to_return = list("lam_1" = coef_left$lam_1 - coef_right$lam_1, 
                     "lam_2" = coef_left$lam_2 - coef_right$lam_2, 
                     "lam_3" = coef_left$lam_3 - coef_right$lam_3)
  } else {
    to_return = list("lam_1" = coef_left$lam_1 - coef_right$lam_1, 
                     "lam_2" = coef_left$lam_2 - coef_right$lam_2, 
                     "lam_3" = coef_left$lam_3 - coef_right$lam_3, 
                     "lam_4" = coef_left$lam_4 - coef_right$lam_4, 
                     "lam_5" = coef_left$lam_5 - coef_right$lam_5)
  }
  return(to_return)
}

fun_roots_sqrt_S_star = function(coef, r_star) {
  # Citation(s):
  # - Uses the function almost.unique from the package bazar
  # Input(s):
  # - coef: list consisting of the coefficients of the function f defined in 
  #   Section A.2.2
  # - r_star: r^star defined in Section 5.2
  # Output(s):
  # - list consisting of 
  #   - root: vector consisting of the roots of f
  #   - elements of coef
  #   - r_star
  toler_roots = 1 # large value to prevent losing roots
  toler_unique = 1e-30 # small value to prevent losing roots
  lam_1 = coef$lam_1
  lam_2 = coef$lam_2
  lam_3 = coef$lam_3
  lam_4 = coef$lam_4
  lam_5 = coef$lam_5
  # gets coefficients in the last equation of Section A.2.2
  quar = lam_3 ** 2 - lam_1 ** 2
  cubic = 2 * (lam_3 * lam_4 - lam_1 * lam_2) 
  quad = lam_4 ** 2 + lam_3 ** 2 * r_star - lam_2 ** 2 - 2 * lam_1 * lam_5 
  lin = 2 * (lam_3 * lam_4 * r_star - lam_2 * lam_5)
  const = lam_4 ** 2 * r_star - lam_5 ** 2
  vec_coef = c(const, lin, quad, cubic, quar)
  # finds solutions to f1^2 = f2^2, where f1 and f2 are defined in Section A.2.2
  root = polyroot(vec_coef)
  root = Re(root[abs(Im(root)) < toler_roots]) # gets the real roots
  root = root[root > 0] # gets positive real roots 
  # finds solutions to f1 = f2 
  if (length(root) == 0) {
  } else {
    root = bazar::almost.unique(root, toler_unique)
    sub_root = rep(0, length(root))
    for (i in 1:length(root)) {
      f1 = (lam_3 * root[i] + lam_4) * sqrt(root[i] ** 2 + r_star)
      f2 = -lam_1 * root[i] ** 2 - lam_2 * root[i] - lam_5
      if (abs(f1 - f2) < toler_roots) { sub_root[i] = 1 } else { 
        sub_root[i] = 0 }
    }
    root = root[sub_root == 1]
    root = sort(root) 
  }
  to_return = list("root" = root, 
                   "lam_1" = lam_1, 
                   "lam_2" = lam_2, 
                   "lam_3"= lam_3, 
                   "lam_4"= lam_4, 
                   "lam_5"= lam_5, 
                   "r_star" = r_star)
  return(to_return)
}

fun_f_S_star = function(y, lam_1, lam_2, lam_3, lam_4, lam_5, r_star) {
  # Input(s):
  # - y: any non-negative real number
  # - lam_1, lam_2, lam_3, lam_4, and lam_5: coefficients in the 
  #   function f defined in Section A.2.2
  # Output(s):
  # - f(y), i.e., the value of f at y 
  term1 = lam_1 * y ** 2
  term2 = lam_2 * y 
  term3 = lam_3 * y * sqrt(y ** 2 + r_star)
  term4 = lam_4 * sqrt(y ** 2 + r_star)
  term5 = lam_5 
  to_return = term1 + term2 + term3 + term4 + term5 
  return(to_return)
}

fun_ineq_sqrt_S_star = function(list_root) {
  # Input(s):
  # - list_root: a list consisting of a vector whose entries are the roots of 
  #   f defined in Section A.2.2, the coefficients of f, and r_star defined in 
  #   Section 5.2
  # Output(s)
  # - vector where the ith and (i+1) entries for every odd i represents the 
  #   endpoints of an interval. The union of all such intervals form a set of 
  #   solutions to f(y)<=0
  root = list_root$root
  lam_1 = list_root$lam_1
  lam_2 = list_root$lam_2
  lam_3 = list_root$lam_3
  lam_4 = list_root$lam_4
  lam_5 = list_root$lam_5
  r_star = list_root$r_star
  if (length(root) == 0) {
    val = fun_f_S_star(1, lam_1, lam_2, lam_3, lam_4, lam_5, r_star)
    if (val <= 0) {
      to_return = c(0, 0)  
    } else {
      print(val)
      print(length(root))
      print(c(lam_1, lam_2, lam_3, lam_4, lam_5))
      to_return = c(-Inf, Inf) 
    }
  } else {
    mat_int = matrix(0, length(root) + 1, 2)
    for (i in 1:length(root)) {
      if (i == 1) {
        val = fun_f_S_star(root[1] / 2, 
                           lam_1, lam_2, lam_3, lam_4, lam_5, r_star)
        if (val <= 0) {
          mat_int[i, ] = c(0, 0)
        } else {
          mat_int[i, ] = c(0, root[1])
        }
      } else {
        val = fun_f_S_star((root[i] + root[i - 1]) / 2, 
                           lam_1, lam_2, lam_3, lam_4, lam_5, r_star)
        if (val <= 0) {
          mat_int[i, ] = c(0, 0) 
        } else {
          mat_int[i, ] = c(root[i - 1], root[i])
        }
      }
    }
    val = fun_f_S_star(root[length(root)] + 5, 
                       lam_1, lam_2, lam_3, lam_4, lam_5, r_star)
    if (val <= 0) {
      mat_int[length(root) + 1, ] = c(0, 0)
    } else {
      mat_int[length(root) + 1, ] = c(root[length(root)], Inf)
    }
    to_return = c(t(mat_int))
  }
  return(to_return)
}

fun_coef_S_star_V = function(X, PE, P1, P2, r_star, v) {
  # Input(s):
  # - X: data matrix of dimensions n by q 
  # - PE: projection matrix P_{mathcal{E}}
  # - P1: projection matrix P_1
  # - P2: projection matrix P_2
  # - r_star: r^star defined in Section 5.2
  # - v: a vector of length n
  # Output(s):
  # - list consisting of the coefficients lambda_{v,1}, lambda_{v,2}, 
  #   lambda_{v,3}, lambda_{v,4}, and lambda_{v,5} of Proposition 8
  A = (PE %*% X) / sqrt(sum((PE %*% X) ** 2))
  B = (P1 %*% X) / sqrt(sum((P1 %*% X) ** 2))
  C = (P2 %*% X) / sqrt(sum((PE %*% X) ** 2) + sum((P1 %*% X) ** 2))
  a = t(A) %*% v 
  b = t(B) %*% v
  c = t(C) %*% v
  lam_1 = sum(a ** 2) + sum(c ** 2)
  lam_2 = 2 * t(a) %*% b * sqrt(r_star)
  lam_3 = 2 * t(a) %*% c 
  lam_4 = 2 * t(b) %*% c * sqrt(r_star)
  lam_5 = (sum(b ** 2) + sum(c ** 2)) * r_star
  to_return = list("lam_1" = lam_1, "lam_2" = lam_2, 
                   "lam_3"= lam_3, "lam_4"= lam_4, "lam_5"= lam_5)
  return(to_return)
}

##------------------------------------------------------------------------------
## Miscellaneous functions for the setting of unknown variance of Section 5

fun_set_K = function(cl, choice_of_V = "pre", set_V = NULL) {
  # Input(s):
  # - cl: vector of cluster assignments
  # - choice_of_V: process by which mathcal{V} is chosen 
  # - set_V: matrix whose columns correspond to the elements of mathcal{V}
  # Output(s):
  # - set_K: mathcal{K} of Section 5.2
  K = max(cl)
  if (choice_of_V == "all") {
    set_K = 1:K
  } else if (choice_of_V %in% c("pre", "closest", "farthest")) {
    if (is.null(set_V)) {
      stop("set_V needs to be specified.")
    } 
    set_K = unique(c(set_V))
  } else {
    stop("The specified choice_of_V does not exist.")
  }
  return(set_K)
}

fun_F_to_chi2 = function(x, n1, n2) {
  # Citations: 
  # - taken from Yun and Barber[2023]'s function fun_F_to_chi2 from 
  #   https://github.com/yjyun97/cluster_inf_unknown_var/blob/main/functions.R
  # Input(s):
  # - x: support of the F distribution
  # - n1: first parameter of the F distribution 
  # - n2: second parameter of the F distribution 
  # Output(s):
  # - the corresponding support of the chi squared distribution
  # Remark(s): 
  # - Uses SFA approximation as described in  
  #   Li, Baibing, and Elaine B. Martin. "An approximation to the F
  #   distribution using the chi-square distribution." Computational 
  #   statistics & data analysis 40, no. 1 (2002): 21-26.
  numRow = nrow(x)
  numCol = ncol(x)
  to_return = matrix(0, numRow, numCol)
  for (i in 1:numRow) {
    for (j in 1:numCol) {
      if (is.infinite(x[i, j]) == 1) {
        to_return[i, j] = Inf
      } else {
        tmp1 = 2 * n2 + (n1 * x[i, j]) / 3 + (n1 - 2)
        tmp2 = 2 * n2 + (4 * n1 * x[i, j]) / 3
        lambda =  tmp1 / tmp2
        to_return[i, j] = lambda * n1 * x[i, j]
      }
    }
  }
  return(to_return)
}

##------------------------------------------------------------------------------
## Function for computing empirical power

fun_power = function(vec_p, thresh) {
  # Citations: 
  # - taken from Yun and Barber[2023]'s function fun_power from
  #   https://github.com/yjyun97/cluster_inf_unknown_var/blob/main/functions.R
  # Input(s):
  # - vec_p: vector of p-values  
  # - thresh: the significance level of the test 
  # Output(s):
  # - empirical power
  # - standard error
  vec_p_new = vec_p[is.na(vec_p) != 1]
  num_p = length(vec_p_new)
  indic = rep(0, num_p)
  indic[vec_p_new < thresh] = 1
  po = sum(indic) / num_p
  eb = sd(indic) / sqrt(num_p) 
  to_return = c(po, eb)
  return(to_return)
}

##------------------------------------------------------------------------------