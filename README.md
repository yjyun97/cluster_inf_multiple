### Summary

This repository contains the codes used for producing the empirical results in the paper "Selective inference for multiple pairs of clusters after K-means clustering" (authors: Youngjoo Yun and Yinqiu He), available at [to be filled in]. 

### Citations 

Some of the codes have been adopted from 

* the source code of the package **kmeans_estimation** written by the authors of [Chen and Witten [2023]](https://jmlr.org/papers/v24/22-0371.html), available at https://github.com/yiqunchen/KmeansInference and 
* the codes written by the authors of [Yun and Barber [2023]](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-17/issue-2/Selective-inference-for-clustering-with-unknown-variance/10.1214/23-EJS2143.full), available at https://github.com/yjyun97/cluster_inf_unknown_var. 

Additionally, the codes make use of functions from the packages **KmeansInference**, **intervals**, **MASS**, **pracma**, **corpcor**, **bazar**, **ggplot2,**  **patchwork**, **ggtext**, **xtable**, and **palmerpenguins**. Citations can be found in the relevant parts of the codes that implement the proposed methods. 

For each empirical result in the paper, the function **mvrnorm** from the package **MASS** is used for generating the data, and the function **kmeans_estimation** from the package **KmeansInference** is used for  running $K$-means clustering on the data. The output of **kmeans_estimation** is then passed in as an input to the inference procedures of the proposed methods. 

### Correspondence to the paper

The files that implement the proposed methods of the paper are located inside the folder **Proposed_methods**, and those that produce the empirical results presented in the paper can be found in the folder **Empirical_results**. We elaborate on the contents of each folder below; all of the files and folders mentioned below are located inside the respective folders. 

#### 1. Proposed methods

The functions that implement the proposed methods of the paper can be found in the files named **functions_pval.R**, **functions_trunc.R**, **functions_helper.R**, and **functions_settings**. Specifically, 

* **functions_pval.R** contains the functions that implement the proposed tests of Section 3.2, Section 4, Section 5.2.1, and Section 5.2.2, along with the baseline testing procedure of Seciton 3.1: the functions  **fun_p_sigma**, **fun_p_sigma_J**, **fun_p_star**, **fun_p_star_J**, and **fun_p_sigma_Bon** produce the p-values $p_\sigma$, $p_{\sigma,J}$, $p^* $,  $p_J^* $, and $p_{\sigma, Bon}$,​​ respectively. These functions take the following arguments as inputs.

  * **X**: the data matrix of dimentions $n$ by $q$.
  * **res_Kmeans**: the result of running $K$-means clustering on $X$, using the function **kmeans_estimation** of the package **KmeansInference** written by the authors of Chen and Witten [2023].
  * **choice_of_V**: the process by which $\mathcal{V}$ is chosen. *"all"* refers to $\mathcal{V}_{\mathrm{all}},$ *"pre"* refers to the setting where $\mathcal{V}$ is pre-specified, *"farthest"* and *"closest"* refer to the settings where $\mathcal{V}$​ is chosen in a data-dependent way according to Setting 1 and Setting 2, respectively, of Section 4. 
  * **set_V**: a matrix of  dimensions 2 by $q$ whose columns consist of elements of $\mathcal{V}.$​ 

The functions **fun_p_sigma**, **fun_p_sigma_J**, and **fun_p_sigma_Bon** further take the arguement **ss,** which can be either $\sigma^2$​​ or its plug-in estimate. 

* **functions_trunc.R** contains the functions that implement the sets $S_\sigma$, $S_{\sigma,\mathcal{V}}$, $S_J$, $S^* $, $S_{\mathcal{V}}^* ,$ and $S^*_J,$ which are used in the computation of the p-values in the functions in **functions_pval.R**. 

* **functions_helper.R** contains helper functions that are called in the functions in **functions_trunc.R** and **functions_pval.R**. 

* **functions_settings.R** contains the functions listed below. 

  * **fun_set_V_dep**:  given $g,$ which refers to the number of pairs of clusters to test for, returns $\mathcal{V}$ chosen in a data-dependent way according to either Setting 1 or Setting 2 of  Section 4. 
  * **fun_set_V_pre**: returns $\mathcal{V}_1$, $\mathcal{V}_2$, or $\mathcal{V}_3$​, whch are defined in Section 6. 
  * **fun_ss_hat_med**: returns a realization of $\hat{\sigma}^2_{\mathrm{med}},$  the estimator of $\sigma^2$ proposed by [Chen and Witten [2023]](https://jmlr.org/papers/v24/22-0371.html).
  * **fun_ss_hat_sample**: returns a realization of $\hat{\sigma}^2_{\mathrm{sample}},$  the estimator of $\sigma^2$ proposed by [Gao et al. [2022]](https://www.tandfonline.com/doi/full/10.1080/01621459.2022.2116331).

#### 2. Empirical results

The files **RData_simul_null.R** and **RData_simul_alter.R** are used for generating the simulation results of Section 6, which are saved as .RData files in the folder **RData**. The files named **figure_i.R** reads in these .RData files to generate **Figure i**, for each **i** in {1,...,15} and saves it in the folder **Figures**. Likewise, the files **table_1.R** and **table_2.R** are used for generating the results on the penguins data of Section 7. The results are saved as .RData files into the folder **RData** and as .txt files in the folder **Tables**. Finally, **specifications_for_plots** contains settings for the visualization of the plots in the figures. 
