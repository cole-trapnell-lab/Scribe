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

cal_grangertest <- function(ordered_exprs_mat) {
  df <- data.frame(ordered_exprs_mat)
  x0_by_x1 <- grangertest(x0 ~ x1, order = 1, data = df)
  x1_by_x0 <- grangertest(x1 ~ x0, order = 1, data = df)
  
  return(data.frame(x0_by_x1 = x0_by_x1$`Pr(>F)`[2], 
                    x1_by_x0 = x1_by_x0$`Pr(>F)`[2]))
}

cal_grangertest_var <- function(ordered_exprs_mat) {
  var <- VAR(ordered_exprs_mat, p = 1, type = "const")
  x1_by_x0 <- causality(var, cause = "x0")$Granger
  x0_by_x1 <- causality(var, cause = "x1")$Granger
  
  return(data.frame(x0_by_x1 = x0_by_x1$p.value, 
                    x1_by_x0 = x1_by_x0$p.value))
}

# this function uses the parmigene package to caculate knn which is very fast (benchma)
cal_mi <- function(ordered_exprs_mat) {
  parmigene::knnmi.all(ordered_exprs_mat)
}

# mi_res <- matrix(1, nrow = ncol(time_series), ncol = ncol(time_series))
# combn_mat <- combn(1:ncol(time_series), 2)

# combn_mat_split <- split(t(combn_mat), 1:ncol(combn_mat))
mi_result <- function(combn_mat_split) {
    mclapply(combn_mat_split, function(x, ordered_exprs_mat){
    col_names <- colnames(ordered_exprs_mat)[x]
    di::mi(time_series[, col_names[1]], time_series[, col_names[2]])
  }, ordered_exprs_mat = time_series, mc.cores = detectCores() - 2)
}