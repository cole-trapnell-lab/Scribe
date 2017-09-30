#' Perform the cross-convergence mapping test for two genes. 
#' 
#' This function uses the CCM function from rEDM package to determine if one time series contain the necessary 
#' dynamic information to recover the influence of another causal variable. See the \package{rEDM} for details
#' of the ccm function. 
#' 
#' @param ordered_exprs_mat CellDataSet for the experiment
#' @param lib_colum The index (or name) of the column to cross map from (same as in rEDM).
#' @param target_column The index (or name) of the column to cross map to (same as in rEDM). 
#' @param RNGseed will set a seed for the random number generator, enabling reproducible runs of ccm with randomly generated libraries (same as in rEDM). 
#' @return a ggplot2 plot object
#' @import rEDM
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' gene_pairs_mat <- matrix(c('H19', 'Ccnd2', 'Ccnd2', 'Scnn1g'), ncol = 2)
#' cal_cross_map(lung, gene_pairs_mat)
#' }
#' @export
#' 
# this file keeps all other relationship estimators: 
cal_cross_map <- function(ordered_exprs_mat, lib_colum, target_column, RNGseed = 2016) { 
  # lib_xmap_target <- ccm(ordered_exprs_mat, E = 3, random_libs = TRUE, lib_column = lib_colum, #ENSG00000122180.4
  #                        target_column = target_column, lib_sizes = seq(10, 75, by = 5), num_samples = 300, RNGseed = RNGseed)
  # target_xmap_lib <- ccm(ordered_exprs_mat, E = 3, random_libs = TRUE, lib_column = target_column, #ENSG00000122180.4
  #                        target_column = lib_colum, lib_sizes = seq(10, 75, by = 5), num_samples = 300, RNGseed = RNGseed)
  # 
  lib_xmap_target <- ccm(ordered_exprs_mat, lib_column = lib_colum, E = 1, #ENSG00000122180.4
                         target_column = target_column, RNGseed = RNGseed)
  target_xmap_lib <- ccm(ordered_exprs_mat, lib_column = target_column, E = 1, #ENSG00000122180.4
                         target_column = lib_colum, RNGseed = RNGseed)
  
  lib_xmap_target_means <- ccm_means(lib_xmap_target)
  target_xmap_lib_means <- ccm_means(target_xmap_lib)
  
  return(data.frame(mean_lib_xmap_target_means = mean(lib_xmap_target_means$rho[is.finite(lib_xmap_target_means$rho)], na.rm = T), 
                    mean_target_xmap_lib_means = mean(target_xmap_lib_means$rho[is.finite(target_xmap_lib_means$rho)], na.rm = T)))
}

#' Perform the Granger causality test for two genes. 
#' 
#' This function uses the grangertest function from lmtest package to determine if one time series contain the necessary 
#' dynamic information to recover the influence of another causal variable. See the \package{lmtest} for details
#' of the grangertest function. 
#' 
#' @param ordered_exprs_mat CellDataSet for the experiment
#' @return a ggplot2 plot object
#' @import lmtest
#' @importFrom lmtest grangertest
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' gene_pairs_mat <- matrix(c('H19', 'Ccnd2', 'Ccnd2', 'Scnn1g'), ncol = 2)
#' cal_cross_map(lung, gene_pairs_mat)
#' }
#' @export
#' 
cal_grangertest <- function(ordered_exprs_mat) {
  df <- data.frame(ordered_exprs_mat)
  x0_by_x1 <- grangertest(x0 ~ x1, order = 1, data = df)
  x1_by_x0 <- grangertest(x1 ~ x0, order = 1, data = df)
  
  return(data.frame(x0_by_x1 = x0_by_x1$`Pr(>F)`[2], 
                    x1_by_x0 = x1_by_x0$`Pr(>F)`[2]))
}

#' Perform the Granger causality test for two genes. 
#' 
#' This function uses the grangertest function from vars package to determine if one time series contain the necessary 
#' dynamic information to recover the influence of another causal variable. See the \package{vars} for details
#' of the VAR or causality function. 
#' 
#' @param ordered_exprs_mat CellDataSet for the experiment
#' @return a ggplot2 plot object
#' @import vars
#' @importFrom vars VAR
#' @importFrom vars causality
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' gene_pairs_mat <- matrix(c('H19', 'Ccnd2', 'Ccnd2', 'Scnn1g'), ncol = 2)
#' cal_cross_map(lung, gene_pairs_mat)
#' }
#' @export
#' 
cal_grangertest_var <- function(ordered_exprs_mat) {
  var <- VAR(ordered_exprs_mat, p = 1, type = "const")
  x1_by_x0 <- causality(var, cause = "x0")$Granger
  x0_by_x1 <- causality(var, cause = "x1")$Granger
  
  return(data.frame(x0_by_x1 = x0_by_x1$p.value, 
                    x1_by_x0 = x1_by_x0$p.value))
}

#' Calculate the mutual information between all gene pairs   
#' 
#' This function uses the knnmi.all function from \package{parmigene} package to calculate the mutual information between all pairs of genes in the data. 
#' See the \package{parmigene} for details of the knnmi.all function. 
#'  
#' @param ordered_exprs_mat CellDataSet for the experiment
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' gene_pairs_mat <- matrix(c('H19', 'Ccnd2', 'Ccnd2', 'Scnn1g'), ncol = 2)
#' plot_scatter_pairs(lung, gene_pairs_mat)
#' }
#' @export
#' 
cal_knn_mi_parmigene <- function(ordered_exprs_mat) {
  parmigene::knnmi.all(ordered_exprs_mat)
}

# mi_res <- matrix(1, nrow = ncol(time_series), ncol = ncol(time_series))
# combn_mat <- combn(1:ncol(time_series), 2)

# combn_mat_split <- split(t(combn_mat), 1:ncol(combn_mat))

#' Calculate the mutual information between all pairs of genes. 
#' 
#' This function calculates the mutual information using the KSG estimators. 
#' 
#' @param combn_mat_split CellDataSet for the experiment
#' @param k Number of nearest neighbors used for calculating the mutual information values. 
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' gene_pairs_mat <- matrix(c('H19', 'Ccnd2', 'Ccnd2', 'Scnn1g'), ncol = 2)
#' plot_scatter_pairs(lung, gene_pairs_mat)
#' }
#' @export
#' 
cal_knn_mi <- function(combn_mat_split, k = 5) {
    mclapply(combn_mat_split, function(x, ordered_exprs_mat, k = k){
    col_names <- colnames(ordered_exprs_mat)[x]
    mi(time_series[, col_names[1]], time_series[, col_names[2]])
  }, ordered_exprs_mat = time_series, k = k, mc.cores = detectCores() - 2)
}