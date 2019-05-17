# ver 0.1; code review and documentation at Oct 2, 2017

# think about better ways to deal with complex tree structure ()

#' Create a pseudo-time-seires
#'
#' This function takes a cds (branched or not) and then convert to a pseudotime series.  
#' @param cds A cds which has been ordered with Monocle 2
#' @param branch_points Vector for the branch points. If it is null, there is no branching in the data 
#' @return a list storing the pseudo-time-series data and a vector storing the run id for each cell (row) in the data 
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' lung_data_list <- buildPseudoTimeSeries(lung, branch_points = 1)
#' }
#' @export
buildPseudoTimeSeries <- function(cds, branch_points = NULL) {
  if(is.null(pData(cds)$Pseudotime))
    stop('Please order your cell dataset with Monocle 2 before running createPTS')
  
  if(!is.integer(branch_points) | any(branch_points > length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points)))
    stop("Branch point should be positive integer and should not be larger than the possible number of branch points. Check branch points with plot_cell_trajectory function in Monocle")
    
  if(is.null(branch_points) == FALSE) { # enable multiple branch points
    data_tmp <- c()
    run_vec_tmp <- c()
    pseudotime <- c()
    run_ind <- 0
    
    for(i in branch_points) {
      tmp <- buildBranchCellDataSet(cds, progenitor_method = 'duplicate', branch_point = i)
      sorted_tmp <- tmp[, order(pData(tmp)$Pseudotime, decreasing = T)]
      level_vec <- unique(pData(sorted_tmp)$Branch)
      
      data_tmp <- rbind(data_tmp, t(as.matrix(exprs(sorted_tmp))))
      Branch <- as.character(pData(sorted_tmp)$Branch)
      Branch[Branch == level_vec[2]] <- run_ind # make sure the first run id is 0 after reversion
      Branch[Branch == level_vec[1]] <- run_ind + 1
      
      run_vec_tmp <- c(run_vec_tmp, Branch)
      pseudotime <- c(pseudotime, pData(tmp)$Pseudotime)
      run_ind <- run_ind + 2
    }
    
    # run 0, 1, ...
    run_vec_tmp <- as.numeric(run_vec_tmp) 
    run_vec <- rev(run_vec_tmp)
    data <- data_tmp[rev(1:ncol(data_tmp)), ]
    
    for(i in unique(run_vec_tmp)) {
        cells_in_current_run <- which(run_vec_tmp == i)
        order_cells_in_current_run_by_pseudotime <- order(pseudotime[cells_in_current_run])
        
        data[cells_in_current_run, ] <- data_tmp[cells_in_current_run[order_cells_in_current_run_by_pseudotime], ]
    }
  }
  else { # no branching in the data 
    data <- t(as.matrix(exprs(cds)))
    data <- data[order(pData(cds)$Pseudotime), ] # order by pseudotime 
    run_vec <- rep(0, ncol(cds))
  }
  
  return(list(data = data, run_vec = run_vec))
}

#' Calculate turning point for a linear or branched trajectory 
#'
#' This function estimates the inflection point or the gene-wise branch point for each gene in a cell trajectory, without or with branch points respectively
#' 
#' @param cds_subset A cds_subset object after trajectory reconstruction
#' @param type A character determines whether or not we will return the order for the turning point or an actual value. By default, it is "order".  
#' @param branch_point If the cds_subset involves a branching process, we pass this argument to determine which branch point we should use for calculating the "turning point"
#' @param cores Numer of cores to run this function 
#' @return A updated cds_subset with a newly added column (turning_point) in pData indicates the inflection or branch time point.
#' @importFrom inflection bede
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' lung_update <- estimate_turning_point(lung, type = 'value')
#' }
#' @export
#' 
#' to do: write a function to integrate BGP with Monocle 
estimate_turning_point <- function(cds_subset, type = c('order', 'value'), branch_point = 1, cores = 1) {
  cds_exprs <- pData(cds_subset)
  cds_exprs$turning_point <- 0

  # branch trajectory 
  if(length(cds_subset@auxOrderingData$DDRTree$branch_points) > 0) {
    b_cds_subset <- buildBranchCellDataSet(cds_subset)
    
    ILRs <- calILRs(cds_subset, branch_point = branch_point, return_all = T, cores = cores) #you can return all important results from the data
    branchtimepoint <- detectBifurcationPoint(ILRs$str_logfc_df, return_cross_point = T)
  
    if(type == 'value')  {
      fData(cds_subset)[, 'turning_point'] <- branchtimepoint
    } else if(type == 'order') {
      # calculate the cell order index corresponding to the branch point (used in RDI calcultion)
      iplasts <- unlist(lapply(branchtimepoint, function(x) {
          branches <- unique(as.character(pData(b_cds_subset)$Branch))

          # take the average of the minimal pseudotimes from the first branch and the second branch 
          a <- which.min(abs(sort(pData(b_cds_subset)[pData(b_cds_subset)$Branch == branches[1], "Pseudotime"]) - abs(x)))
          b <- which.min(abs(sort(pData(b_cds_subset)[pData(b_cds_subset)$Branch == branches[2], "Pseudotime"]) - abs(x)))
          
          round(mean(c(a, b), na.rm = T))
        }
      )) # get the order of the cell corresponding the inflection point)
      fData(cds_subset)[, 'turning_point'] <- iplasts
    }
  } else{ # linear trajectory 
    rng <- range(pData(cds_subset)$Pseudotime)
    new_data <- data.frame(Pseudotime = seq(rng[1], rng[2], length.out = 100))
    cds_exprs <- genSmoothCurves(cds_subset, cores = cores, trend_formula = "~sm.ns(Pseudotime, df = 3)",
                    relative_expr = TRUE, new_data = new_data)
    for(i in 1:nrow(cds_exprs)) {
      data <- cds_exprs[i, ]; 
      inflection_point <- bede(1:length(data), data, 0)
      if(!is.finite(inflection_point$iplast))
        inflection_point <- bede(1:length(data), data, 1)
    
      # message('turning_point is ', inflection_point$iplast)
      iplast  <- inflection_point$iplast
      
      if(type == 'value')  { 
        fData(cds_subset)[i, 'turning_point'] <- iplast 
      } else if(type == 'order') {
      # get the order of the cell corresponding the inflection point 
       fData(cds_subset)[i, 'turning_point'] <- which.min(abs(sort(pData(cds_subset)$Pseudotime) - new_data$Pseudotime[iplast])) 
      }
    }
  }
  
  return(cds_subset)
}

#' Calculate uMI values
#'
#' This function estimates the uniformed mutual information value for all gene-pair combination from the genes_data (or that specified in the super_graph) 
#' 
#' @param cds_subset A cds_subset which has been ordered with Monocle 2
#' @param super_graph A graph including all possible interactions used for performing the causality inference. When it is NULL, 
#' all possible pairs of interactions for the genes_data are used, otherwise only the interactions specified in the graph are used.
#' Note that the super_graph only accept integer matrix for now (each integer corresponds to a particular gene in the genes_data).  
#' @param log
#' @param pseudo_cnt
#' @param k Number for nearest neighbors used in entropy calculation.
#' @param methodWhich 2D density estimator you would like to use. 1 is kde estimator and 2 is knn based estimator. Default to be 1. 
#' @param k_density The number of k nearest neighbors you would like to use when calculating the density (only applicable when method == 2 or using knn based density estimation).
#' @param bw Bindwidth used for the kernel density estimator. Currently it is not used. The bindwidth in the kde function is automatically estimated. 
#' @return a numeric matrix storing a uMI result for all gene-pair combination from the genes_data (or that specified in the super_graph). 
#'
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' AT1_lung <- lung[, pData(lung)$State %in% c(2, 3)]
#' rdi_res <- calculate_umi(AT1_lung, delays = 10, method = 1)
#' }
#' @export
#' 
#' Third column is the conditional RDI value.
#' 
calculate_umi <- function(cds_subset, super_graph = NULL, log = TRUE, pseudo_cnt = 1, k = 5L, method = 1L, k_density = 5L, bw = 0.0) {
  
  if(length(unique(pData(cds_subset)$State)) > 1)
    warning('calculate_rdi only deals with linear trajectory but Your cds seems like a branched trajectory!')
  
  sz <- sizeFactors(cds_subset)
  sz[is.na(sz)] <- 1
  
  if(log) {
    genes_data <- log(t(as.matrix(exprs(cds_subset))) / sz + pseudo_cnt) # Rows are samples while columns are features
  }
  else {
    genes_data <- t(as.matrix(exprs(cds_subset))) / sz # Rows are samples while columns are features
  }
  if(!all(is.finite(genes_data))) {
    stop('your data includes non finite values!')
    # finite_genes <- apply(genes_data, 2, function(x) all(is.finite(x)))
    # genes_data <- genes_data[, finite_genes] # remove non-infinite gene 
  }
  
  n_genes <- ncol(genes_data); 
  n_samples <- nrow(genes_data); 
  # set.seed(2017)
  # noise = matrix(rnorm(mean = 0, sd = 1e-10, n_genes * n_samples), nrow = n_samples)
  
  if(is.null(super_graph)) 
  {
    tmp <- expand.grid(1:n_genes, 1:n_genes, stringsAsFactors = F)
    super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to c++ index
  }
  else
  {
    if(class(super_graph[1, 1]) != "character"){
      if(max(super_graph) >  n_genes - 1 | min(super_graph) < 0)
      {
        stop("super_graph should only include integer less than the number of cells (or only include gene names from genes_data)")
      }
    }
    else {
      if(!(all(unlist(super_graph) %in% colnames(genes_data))))
      {
        stop("super_graph should only include integer less than the number of cells (or only include gene names from genes_data)")
      }
      super_graph[, 1] <- match(super_graph[, 1], colnames(genes_data))
      super_graph[, 2] <- match(super_graph[, 2], colnames(genes_data))
    }
  }
  
  uMI <- calculate_umi_cpp_wrap(as.matrix(genes_data), as.matrix(super_graph), k = k, method = method, k_density = k_density, bw = bw) # + noise
  
  # assign the name for each dimnames
  dimnames(uMI) <- list(colnames(genes_data), colnames(genes_data)) 
  
  return(uMI); 
}

#' Calculate RDI values
#'
#' This function estimates the RDI value for all gene-pair combination from the genes_data (or that specified in the super_graph) 
#' in the pseudotime/time series data.
#' 
#' @param cds_subset A cds_subset which has been ordered with Monocle 2
#' @param delays A vector of time delays used during information transfer estimation between genes 
#' @param super_graph A graph including all possible interactions used for performing the causality inference. When it is NULL, 
#' all possible pairs of interactions for the genes_data are used, otherwise only the interactions specified in the graph are used.
#' Note that the super_graph only accept integer matrix for now (each integer corresponds to a particular gene in the genes_data).  
#' @param turning_points pseudo-time/time series for the gene expression data
#' @param method a list returned from calculate_rdi function. 
#' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
#' @param log
#' @param pseudo_cnt
#' @return a list storing a matrix for RDI results, maximal rdi values and the delays corresponding to the maximal rdi values. 
#'
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' AT1_lung <- lung[, pData(lung)$State %in% c(2, 3)]
#' rdi_res <- calculate_rdi(AT1_lung, delays = 10, method = 1)
#' }
#' @export
#' 
#' Third column is the conditional RDI value.
#' 
calculate_rdi <- function(cds_subset, delays, super_graph = NULL, turning_points = 0, method = 1, uniformalize = FALSE, log = TRUE, pseudo_cnt = 1, smoothing = FALSE) {

  if(length(unique(pData(cds_subset)$State)) > 1)
    warning('calculate_rdi only deals with linear trajectory but Your cds seems like a branched trajectory!')
  
  sz <- sizeFactors(cds_subset)
  sz[is.na(sz)] <- 1
  
  if(log) {
    genes_data <- log(t(as.matrix(exprs(cds_subset)[, order(pData(cds_subset)$Pseudotime)])) / sz[order(pData(cds_subset)$Pseudotime)] + pseudo_cnt) # Rows are samples while columns are features
  }
  else {
    genes_data <- t(as.matrix(exprs(cds_subset)[, order(pData(cds_subset)$Pseudotime)])) / sz[order(pData(cds_subset)$Pseudotime)] # Rows are samples while columns are features
  }
  if(smoothing) {
    Pseudotime <- pData(cds_subset)$Pseudotime 
    
    smoothed_exprs <- apply(genes_data, 2, function(x) {
      df <- data.frame(Expression = x, 
                       Pseudotime = Pseudotime)
      fit <- loess(Expression ~ Pseudotime, df) #, span = 0.1
      predict(fit)
    })
    
    genes_data <- smoothed_exprs
  }
  if(!all(is.finite(genes_data))) {
    stop('your data includes non finite values!')
    # finite_genes <- apply(genes_data, 2, function(x) all(is.finite(x)))
    # genes_data <- genes_data[, finite_genes] # remove non-infinite gene 
  }
  
  n_genes <- ncol(genes_data); 
  n_samples <- nrow(genes_data); 
  # set.seed(2017)
  # noise = matrix(rnorm(mean = 0, sd = 1e-10, n_genes * n_samples), nrow = n_samples)
  
  if(is.null(super_graph)) 
  {
    tmp <- expand.grid(1:n_genes, 1:n_genes, stringsAsFactors = F)
    super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to c++ index
  }
  else
  {
    if(class(super_graph[1, 1]) != "character"){
      if(max(super_graph) >  n_genes - 1 | min(super_graph) < 0)
      {
        stop("super_graph should only include integer less than the number of cells (or only include gene names from genes_data)")
      }
    }
    else {
      if(!(all(unlist(super_graph) %in% colnames(genes_data))))
      {
        stop("super_graph should only include integer less than the number of cells (or only include gene names from genes_data)")
      }
      super_graph[, 1] <- match(super_graph[, 1], colnames(genes_data))
      super_graph[, 2] <- match(super_graph[, 2], colnames(genes_data))
    }
  }

  RDI_list <- calculate_rdi_cpp_wrap(as.matrix(genes_data), delays, as.matrix(super_graph), turning_points, method, uniformalize) # + noise
  
  # assign the name for each dimnames
  dimnames(RDI_list$RDI) <- list(colnames(genes_data), rep(colnames(genes_data), length(delays))) 
  dimnames(RDI_list$max_rdi_value) <- list(colnames(genes_data), colnames(genes_data))
  dimnames(RDI_list$max_rdi_delays) <- list(colnames(genes_data), colnames(genes_data)) 

  return(RDI_list); 
}

# think about passing a pseudotime ordered cds instead of a matrix 
#' Calculate conditionally RDI value
#'
#' This function estimates the conditional RDI value for all gene-pair combination from the genes_data (or that specified in the super_graph) 
#' in the pseudotime/time series data, conditioned on the top top_incoming_k incoming nodes. 
#' The result from the calculate_rdi function will be used to identify the proper delay for the gene-pair under test corresponds to that 
#' with highest RDI as well as the proper delays for the k incoming nodes which corresponds to that with highest RDI values.
#' 
#' @param cds_subset A cds_subset which has been ordered with Monocle 2
#' @param super_graph A graph including all possible interactions used for performing the causality inference. When it is NULL, 
#' all possible pairs of interactions for the genes_data are used, otherwise only the interactions specified in the graph are used.
#' Note that the super_graph only accept integer matrix for now (each integer corresponds to a particular gene in the genes_data).  
#' @param rdi_list a list returned from calculate_rdi function. 
#' @param top_incoming_k The number of genes to be conditioned on when calculating the conditional RDI values 
#' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
#' @return a dataframe storing conditional RDI results. First two columns are the id names for the genes.
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' rdi_res <- calculate_rdi(exprs(lung)[, order(pData(lung)$Pseudotime)], delays = 5)
#' lung_res_cRDI <- calculate_conditioned_rdi(lung, rdi_list = rdi_res)
#' }
#' @export
#' 
#' Third column is the conditional RDI value.
#' 
calculate_conditioned_rdi <- function(cds_subset, super_graph = NULL, rdi_list, top_incoming_k = 1, uniformalize = FALSE, log = TRUE, pseudo_cnt = 1) {
  if(length(unique(pData(cds_subset)$State)) > 1)
    warning('calculate_rdi only deals with linear trajectory but Your cds_subset seems like a branched trajectory!')
  
  sz <- sizeFactors(cds_subset)
  sz[is.na(sz)] <- 1
  
  if(log) {
    genes_data <- log(t(as.matrix(exprs(cds_subset)[, order(pData(cds_subset)$Pseudotime)])) / sz[order(pData(cds_subset)$Pseudotime)] + pseudo_cnt) # Rows are samples while columns are features
  } else {
    genes_data <- t(as.matrix(exprs(cds_subset)[, order(pData(cds_subset)$Pseudotime)])) / sz[order(pData(cds_subset)$Pseudotime)] # Rows are samples while columns are features
  }
  if(!all(is.finite(genes_data))) {
    stop('your data includes non finite values!')
    # finite_genes <- apply(genes_data, 2, function(x) all(is.finite(x)))
    # genes_data <- genes_data[, finite_genes] # remove non-infinite gene 
  }
  
  n_genes <- ncol(genes_data); 
  n_samples <- nrow(genes_data); 
  # set.seed(2017)
  # noise = matrix(rnorm(mean = 0, sd = 1e-10, n_genes * n_samples), nrow = n_samples)
  
  if(is.null(super_graph)) 
  {
    tmp <- expand.grid(1:n_genes, 1:n_genes, stringsAsFactors = F)
    super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to c++ index
  }
  else
  {
    if(class(super_graph[1, 1]) != "character"){
      if(max(super_graph) > nrow(genes_data) | min(super_graph) < 0)
        stop("super_graph should only include gene names from genes_data")
    }
    else {
      if((! all(unique(super_graph) %in% colnames(genes_data))))
        stop("super_graph should only include gene names from genes_data")
      
      super_graph[, 1] <- match(super_graph[, 1], colnames(genes_data))
      super_graph[, 2] <- match(super_graph[, 2], colnames(genes_data))
    }
  }
  max_rdi_value <- rdi_list$max_rdi_value
  max_rdi_delays <- rdi_list$max_rdi_delays
  cRDI_mat <- calculate_conditioned_rdi_cpp_wrap(genes_data, as.matrix(super_graph), # + noise
                                            max_rdi_value, max_rdi_delays, top_incoming_k, uniformalize)

  dimnames(cRDI_mat) <- list(colnames(genes_data), colnames(genes_data))
  
  return(cRDI_mat); 
}

####################################################################################################################################################################################
# update this to multiple run 
####################################################################################################################################################################################

#' @title
#' calculate_rdi_multiple_run
#' @description
#' This subroutine estimates the RDI value for all gene-pair combination from the genes_data (or that specified in the super_graph) 
#' in the pseudotime/time series data. If multiple directly related experiments are conducted, we can pass the information of each 
#' experiment by the run_vec argument and this function will concatenate those experiments based on run_vec. 
#' 
#' @param cds_subset A cds_subset which has been ordered with Monocle 2
#' @param run_vec A vector describes which run (lineage) does the current cell come from. It has the some length as the row length of genes_data. 
#' @param delays Time lags used to estimate the RDI values (be ready to update this as a matrix or a vector with the same length as the run_vec once we have a way to directly calculate the proper time lag need to be used)
#' @param super_graph A graph including all possible interactions used for performing the causality inference. When it is NULL, 
#' all possible pairs of interactions for the genes_data are used, otherwise only the interactions specified in the graph are used.
#' Note that the super_graph only accept integer matrix for now (each integer corresponds to a particular gene in the genes_data).  
#' @param method either 1, 2 represents either using (RDI: restricted direction information) or LMI (lagged mutual information) 
#' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
#'
#' @details
#' \code{calculate_rdi_multiple_run} takes a time / Pseudo-time series expression data as well as the time lags and then calculate the restricted direct information between each pair of genes under different delays. 
#' @return a list with four components: matrix for RDI (dimension is number of genes X length of delays times number of genes), vector of delays, max_rdi_value and max_rdi_delays 
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' lung_dup <- buildBranchCellDataSet(lung)
#' run_vec <- pData(lung_dup)$Branch
#' Branch_uniq <- unique(run_vec)
#' run_vec[run_vec == Branch_uniq[1]] <- 0
#' run_vec[run_vec == Branch_uniq[2]] <- 1
#' lung_res <- calculate_rdi_multiple_run(lung, run_vec = run_vec)
#' }
#' @export
#' 
calculate_rdi_multiple_run <- function(cds_subset, run_vec = NULL, delays, super_graph = NULL, method = 1, uniformalize = FALSE, log = TRUE, pseudo_cnt = 1) {
  sz <- sizeFactors(cds_subset)
  sz[is.na(sz)] <- 1
  
  if(log) {
    genes_data <- log(t(as.matrix(exprs(cds_subset)[, order(pData(cds_subset)$Pseudotime)])) / sz[order(pData(cds_subset)$Pseudotime)] + 1) # Rows are samples while columns are features
  } else {
    genes_data <- t(as.matrix(exprs(cds_subset)[, order(pData(cds_subset)$Pseudotime)]))  / sz[order(pData(cds_subset)$Pseudotime)] # Rows are samples while columns are features
  }
  
  if(!all(is.finite(genes_data))) {
    stop('your data includes non finite values!')
    # finite_genes <- apply(genes_data, 2, function(x) all(is.finite(x)))
    # genes_data <- genes_data[, finite_genes] # remove non-infinite gene 
  }
  
  if(is.null(run_vec)) {
    run_vec <- rep(0, nrow(genes_data)) 
  } else {
    run_vec <- run_vec[order(pData(cds_subset)$Pseudotime)]
  }

  if(min(run_vec) != 0 | any(diff(unique(run_vec)) != 1) | !is.vector(run_vec))
    stop("The run_vec starts from 0 to the total number of runs")
  
  n_genes <- ncol(genes_data); 
  n_samples <- nrow(genes_data); 
  # set.seed(2017)
  # noise = matrix(rnorm(mean = 0, sd = 1e-10, n_genes * n_samples), nrow = n_samples)
  
  if(is.null(super_graph)) 
  {
    tmp <- expand.grid(1:n_genes, 1:n_genes, stringsAsFactors = F)
    super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to c++ index
  }
  else
  {
    if(class(super_graph[1, 1]) != "character"){
      if(max(super_graph) >  n_genes - 1 | min(super_graph) < 0)
      {
        stop("super_graph should only include integer less than the number of genes (or only include gene names from genes_data)")
      }
    }
    else {
      if(!(all(unlist(super_graph) %in% colnames(genes_data))))
      {
        stop("super_graph should only include integer less than the number of cells (or only include gene names from genes_data)")
      }
      super_graph[, 1] <- match(super_graph[, 1], colnames(genes_data))
      super_graph[, 2] <- match(super_graph[, 2], colnames(genes_data))
    }
  }
  
  RDI_list <- calculate_rdi_multiple_run_cpp_wrap(as.matrix(genes_data), delays, run_vec, as.matrix(super_graph), turning_points = 0, method, uniformalize) 
  
  # assign the name for each dimnames
  dimnames(RDI_list$RDI) <- list(colnames(genes_data), rep(colnames(genes_data), length(delays))) 
  dimnames(RDI_list$max_rdi_value) <- list(colnames(genes_data), colnames(genes_data))
  dimnames(RDI_list$max_rdi_delays) <- list(colnames(genes_data), colnames(genes_data)) 
  
  return(RDI_list); 
}


#' Calculate conditionally RDI value for multiple runs (replicates) of experiments 
#'
#' This function estimates the RDI value for all gene-pair combination from the genes_data (or that specified in the super_graph) 
#' in the pseudotime/time series data, conditioned on the top top_incoming_k incoming nodes. If multiple directly related experiments 
#' are conducted, we can pass the information of each experiment by the run_vec argument and this function will concatenate those experiments
#' based on run_vec. 
#' 
#' @param cds_subset A cds_subset which has been ordered with Monocle 2
#' @param run_vec A vector describes which run (lineage) does the current cell come from. It has the some length as the row length of genes_data. 
#' @param super_graph A graph including all possible interactions used for performing the causality inference. When it is NULL, 
#' all possible pairs of interactions for the genes_data are used, otherwise only the interactions specified in the graph are used.
#' Note that the super_graph only accept integer matrix for now (each integer corresponds to a particular gene in the genes_data).  
#' @param rdi_list a list returned from calculate_rdi function. 
#' @param top_incoming_k The number of genes to be conditioned on when calculating the conditional RDI values 
#' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
#' @return a dataframe storing conditional RDI results. First two columns are the id names for the genes.
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' lung_dup <- buildBranchCellDataSet(lung)
#' run_vec <- pData(lung_dup)$Branch
#' Branch_uniq <- unique(run_vec)
#' run_vec[run_vec == Branch_uniq[1]] <- 0
#' run_vec[run_vec == Branch_uniq[2]] <- 1
#' lung_res <- calculate_rdi_multiple_run(lung, run_vec = run_vec)
#' lung_res_cRDI <- calculate_rdi_multiple_run(lung, run_vec = run_vec, rdi_list = lung_res)
#' }
#' @export
#' Third column is the conditional RDI value.
#' 
calculate_conditioned_rdi_multiple_run <- function(cds_subset, run_vec = NULL, super_graph = NULL, rdi_list, top_incoming_k = 1, uniformalize = FALSE, log = TRUE, pseudo_cnt = 1) {
  sz <- sizeFactors(cds_subset)
  sz[is.na(sz)] <- 1
  
  if(log) {
    genes_data <- log(t(as.matrix(exprs(cds_subset)[, order(pData(cds_subset)$Pseudotime)])) / sz[order(pData(cds_subset)$Pseudotime)] + pseudo_cnt) # Rows are samples while columns are features
  } else {
    genes_data <- t(as.matrix(exprs(cds_subset)[, order(pData(cds_subset)$Pseudotime)])) / sz[order(pData(cds_subset)$Pseudotime)] # Rows are samples while columns are features
  }
  
  if(!all(is.finite(genes_data))) {
    stop('your data includes non finite values!')
    # finite_genes <- apply(genes_data, 2, function(x) all(is.finite(x)))
    # genes_data <- genes_data[, finite_genes] # remove non-infinite gene 
  }
  
  if(is.null(run_vec)) {
    run_vec <- rep(0, nrow(genes_data)) 
  } else {
    run_vec <- run_vec[order(pData(cds_subset)$Pseudotime)]
  }

  if(min(run_vec) != 0 | any(diff(unique(run_vec)) != 1) | !is.vector(run_vec))
    stop("The run_vec starts from 0 to the total number of runs")
  
  n_genes <- ncol(genes_data); 
  n_samples <- nrow(genes_data); 
  # set.seed(2017)
  # noise = matrix(rnorm(mean = 0, sd = 1e-10, n_genes * n_samples), nrow = n_samples)
  
  if(is.null(super_graph)) 
  {
    tmp <- expand.grid(1:n_genes, 1:n_genes, stringsAsFactors = F)
    super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to c++ index
  }
  else
  {
    if(class(super_graph[1, 1]) != "character"){
      if(max(super_graph) > nrow(genes_data) | min(super_graph) < 0)
        stop("super_graph should only include gene names from genes_data")
    }
    else {
      if((! all(unique(super_graph) %in% colnames(genes_data))))
        stop("super_graph should only include gene names from genes_data")
      
      super_graph[, 1] <- match(super_graph[, 1], colnames(genes_data))
      super_graph[, 2] <- match(super_graph[, 2], colnames(genes_data))
    }
  }
  max_rdi_value <- rdi_list$max_rdi_value
  max_rdi_delays <- rdi_list$max_rdi_delays
  cRDI_mat <- calculate_conditioned_rdi_multiple_run_wrap(genes_data, as.matrix(super_graph), 
                                                 max_rdi_value, max_rdi_delays, run_vec, top_incoming_k, uniformalize) # + noise
  
  dimnames(cRDI_mat) <- list(colnames(genes_data), colnames(genes_data))
  
  return(cRDI_mat); 
}

#' Calculate temporal RDI or conditional RDI (cRDI) values. 
#' 
#' This function is developed to calculate the temporal RDI or cRDI which can reveal the time point where the regulator will affect the target gene.  
#' 
#' @param cds_subset A cds_subset which has been ordered with Monocle 2. Currently it only deals with a linear trajectory. 
#' @param super_graph A matrix for the valid gene-pairs we will consider for inferring the temporal RDI or cRDI. 
#' @param conditioning A logic argument to determine whether or not the conditional RDI will be calculated. 
#' @param window_size The window size we will use for calculating the temporal RDI. 
#' @param delay The time delay used in calculating the RDI or cRDI values.
#' @param smoothing A logic argument to determine whether or not we will first apply loess smoothing before we calculate the temporal RDI. Default to be TRUE. 
#' @param span The parameter α which controls the degree of smoothing, passed into the loess functioin. Default is set to be 0.1. 
#' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
#' @param verbose A logical argument to determine whether or not we will print the detailed running information. 
#' @return A list of two arrays where the first one is the temporal RDI value and second one is the conditional RDI values. In each array, the first 
#' dimension corresponds to the sliding window, the second and third dimensions correspond to the genes.  
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' plot_cell_trajectory(lung)
#' plot_cell_trajectory(lung, color_by = 'Time')
#' lung_AT1_lineage <- lung[, pData(lung)$State %in% c(1, 3)]
#' temp_causality_list <- calculate_temporal_causality(lung_AT1_lineage)
#' }
#' @export
#' 
calculate_temporal_causality <- function(cds_subset, super_graph = NULL, conditioning = FALSE, window_size = 40, delay = 1, smoothing = TRUE, 
                                         span = 0.1, uniformalize = FALSE, log = TRUE, pseudo_cnt = 1, verbose = FALSE, ...) {
  sz <- sizeFactors(cds_subset)
  sz[is.na(sz)] <- 1
  
  if(log == TRUE) {
    data <- log(t(as.matrix(exprs(cds_subset[, order(pData(cds_subset)$Pseudotime)]))) / sz[order(pData(cds_subset)$Pseudotime)] + pseudo_cnt)
  } else {
    data <- t(as.matrix(exprs(cds_subset[, order(pData(cds_subset)$Pseudotime)]))) / sz[order(pData(cds_subset)$Pseudotime)]
  }
  
  if(smoothing) {
    for(i in 1:ncol(data)) {
      df <- data.frame(Pseudotime = 1:nrow(data), Expression = data[, i])
      test <- loess(Expression ~ Pseudotime, df, span = span)
      data[, i] <-  predict(test)
    }
  }
  
  win_range <- nrow(data) - window_size - 1
  gene_num <- ncol(data)
  
  window_gene_gene_result <- array(dim = c(win_range + 1, gene_num, gene_num)) 
  crdi_window_gene_gene_result <- window_gene_gene_result 
  
  run_vec <- rep(1, nrow(data))
  
  if(is.null(super_graph)) {
    tmp <- expand.grid(1:ncol(data), 1:ncol(data), stringsAsFactors = F)
    super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index
  }

  for(i in 1:(win_range + 1)) {
    if(verbose)
     message('current window index is ', i)
    
    rdi_list <- calculate_rdi_multiple_run_cpp(expr_data = data[i:(i + window_size), ], delay = delay, run_vec = run_vec[i:(i + window_size)] - 1, 
                                               super_graph = as.matrix(super_graph), turning_points = 0, method = 1, uniformalize) #* 100 + noise
    
    if(conditioning) {
      con_rdi_res_test <- calculate_conditioned_rdi_multiple_run_wrap(data[i:(i + window_size), ], as.matrix(super_graph), as.matrix(rdi_list$max_rdi_value), 
                                                                      as.matrix(rdi_list$max_rdi_delays), run_vec[i:(i + window_size)] - 1, 1, uniformalize)
      
      crdi_window_gene_gene_result[i, , ] <- con_rdi_res_test
    } else {
      con_rdi_res_test <- rep(NA)
    }
    
    window_gene_gene_result[i, , ] <- rdi_list$max_rdi_value
    
  }
  
  return(list(rdi_res = window_gene_gene_result, crdi_res = crdi_window_gene_gene_result, gene_name_vec = row.names(cds_subset)))
}

#' Calculate temporal RDI or conditional RDI (cRDI) values. 
#' 
#' This function is developed to calculate the temporal RDI or cRDI which can reveal the time point where the regulator will affect the target gene.  
#' 
#' @param cds A cds which contains the genes from the TFs and informative_genes vectors in the gene_short_name column of the cds's pData. 
#' @param TFs A vector for the gene short name of the TFs, default to be NULL (All genes in the cds will be used). 
#' @param informative_genes A vector for the gene short name of the informative genes, default to be NULL (All genes in the cds will be used). 
#' @return A character matrix where the first column corresponds to the regulators and the second column corresponds to the targets.  
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' super_graph <- create_super_graph(lung)
#' }
#' @export
create_super_graph <- function(cds, TFs = NULL, informative_genes = NULL) {
  if(is.null(TFs) & is.null(informative_genes)) {
    tmp <- expand.grid(row.names(cds), row.names(cds), stringsAsFactors = F)
    super_graph <- tmp[tmp[, 1] != tmp[, 2], ]
  } else {
    # convert the genes into character vectors 
    TFs <- as.character(TFs)
    informative_genes <- as.character(informative_genes)
    
    if(!all(c(TFs, informative_genes) %in% fData(cds)$gene_short_name)) {
      stop('Not all your TFs and informative_genes are included in the gene_short_name column of the fData.')
    } else if (length(intersect(TFs, informative_genes)) < 1) {
      stop('There should be at least one TFs gene is included in the vector informative_gene')
      } else {
      cds <- cds[row.names(subset(fData(cds), gene_short_name %in% c(TFs, informative_genes))), ]
    }
    
    # obtain the gene ensemble IDs: 
    TFs <- row.names(subset(fData(cds), gene_short_name %in% TFs))
    informative_genes <- row.names(subset(fData(cds), gene_short_name %in% informative_genes))
    
    TF_vec_names <- intersect(TFs, informative_genes)
    target_vec_names <- setdiff(informative_genes, TFs)
    
    TF_pair <- expand.grid(TF_vec_names, TF_vec_names, stringsAsFactors = F) # between TFs 
    TF_target_pair <- expand.grid(TF_vec_names, target_vec_names, stringsAsFactors = F) # from TFs to the targets 
    
    super_graph <- rbind(TF_pair, TF_target_pair)
  }
  
  return(super_graph)
}

# this following function used to re-order the trajectory using principal curve 
reduceDimension <- function(cds,
                            max_components=2,
                            reduction_method=c("DDRTree", "ICA", 'tSNE', "SimplePPT", 'L1-graph', 'SGL-tree', "Principal.curve"),
                            norm_method = c("log", "vstExprs", "none"),
                            residualModelFormulaStr=NULL,
                            pseudo_expr=NULL,
                            relative_expr=TRUE,
                            auto_param_selection = TRUE,
                            verbose=FALSE,
                            scaling = TRUE,
                            ...){
  extra_arguments <- list(...)
  set.seed(2016) #ensure results from RNG sensitive algorithms are the same on all calls
  
  FM <- monocle:::normalize_expr_data(cds, norm_method, pseudo_expr)
  
  #FM <- FM[unlist(sparseApply(FM, 1, sd, convert_to_dense=TRUE)) > 0, ]
  xm <- Matrix::rowMeans(FM)
  xsd <- sqrt(Matrix::rowMeans((FM - xm)^2))
  FM <- FM[xsd > 0,]
  
  if (is.null(residualModelFormulaStr) == FALSE) {
    if (verbose)
      message("Removing batch effects")
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr),
                                       data = pData(cds), drop.unused.levels = TRUE)
    
    fit <- limma::lmFit(FM, X.model_mat, ...)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
  }else{
    X.model_mat <- NULL
  }
  
  if(scaling){
    FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
    FM <- FM[!is.na(row.names(FM)), ]
  } else FM <- as.matrix(FM)
  
  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }
  
  FM <- FM[apply(FM, 1, function(x) all(is.finite(x))), ] #ensure all the expression values are finite values
  if (is.function(reduction_method)) {
    reducedDim <- reduction_method(FM, ...)
    colnames(reducedDim) <- colnames(FM)
    monocle:::reducedDimW(cds) <- as.matrix(reducedDim)
    monocle:::reducedDimA(cds) <- as.matrix(reducedDim)
    monocle:::reducedDimS(cds) <- as.matrix(reducedDim)
    monocle:::reducedDimK(cds) <- as.matrix(reducedDim)
    dp <- as.matrix(dist(reducedDim))
    cellPairwiseDistances(cds) <- dp
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    minSpanningTree(cds) <- dp_mst
    cds@dim_reduce_type <- "function_passed"
  }
  else{
    reduction_method <- match.arg(reduction_method)
    if (reduction_method == "tSNE") {
      #first perform PCA
      if (verbose)
        message("Remove noise by PCA ...")
      
      # # Calculate the variance across genes without converting to a dense
      # # matrix:
      # FM_t <- Matrix::t(FM)
      # cell_means <- Matrix::rowMeans(FM_t)
      # cell_vars <- Matrix::rowMeans((FM_t - cell_means)^2)
      # Filter out genes that are constant across all cells:
      #genes_to_keep <- expression_vars > 0
      #FM <- FM[genes_to_keep,]
      #expression_means <- expression_means[genes_to_keep]
      #expression_vars <- expression_vars[genes_to_keep]
      # Here✬s how to take the top PCA loading genes, but using
      # sparseMatrix operations the whole time, using irlba.
      
      
      if("num_dim" %in% names(extra_arguments)){ #when you pass pca_dim to the function, the number of dimension used for tSNE dimension reduction is used
        num_dim <- extra_arguments$num_dim #variance_explained
      }
      else{
        num_dim <- 50
      }
      
      irlba_res <- prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1),
                                center = TRUE, scale. = TRUE)
      irlba_pca_res <- irlba_res$x
      
      # irlba_res <- irlba(FM,
      #                    nv=min(num_dim, min(dim(FM)) - 1),
      #                        nu=0,
      #                        center=cell_means,
      #                        scale=sqrt(cell_vars),
      #                        right_only=TRUE)
      # irlba_pca_res <- irlba_res$v
      # row.names(irlba_pca_res) <- genes_to_keep
      
      # pca_res <- prcomp(t(FM), center = T, scale = T)
      # std_dev <- pca_res$sdev
      # pr_var <- std_dev^2
      # prop_varex <- pr_var/sum(pr_var)
      # prop_varex <- irlba_res$sdev^2 / sum(irlba_res$sdev^2)
      
      topDim_pca <- irlba_pca_res#[, 1:num_dim]
      
      # #perform the model formula transformation right before tSNE:
      # if (is.null(residualModelFormulaStr) == FALSE) {
      #   if (verbose)
      #     message("Removing batch effects")
      #   X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr),
      #                                      data = pData(cds), drop.unused.levels = TRUE)
      
      #   fit <- limma::lmFit(topDim_pca, X.model_mat, ...)
      #   beta <- fit$coefficients[, -1, drop = FALSE]
      #   beta[is.na(beta)] <- 0
      #   topDim_pca <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
      # }else{
      #   X.model_mat <- NULL
      # }
      
      #then run tSNE
      if (verbose)
        message("Reduce dimension by tSNE ...")
      
      tsne_res <- Rtsne::Rtsne(as.matrix(topDim_pca), dims = max_components, pca = F,...)
      
      tsne_data <- tsne_res$Y[, 1:max_components]
      row.names(tsne_data) <- colnames(tsne_data)
      
      reducedDimA(cds) <- t(tsne_data) #this may move to the auxClusteringData environment
      
      #set the important information from densityClust to certain part of the cds object:
      cds@auxClusteringData[["tSNE"]]$pca_components_used <- num_dim
      cds@auxClusteringData[["tSNE"]]$reduced_dimension <- t(tsne_data)
      #cds@auxClusteringData[["tSNE"]]$variance_explained <- prop_varex
      
      cds@dim_reduce_type <- "tSNE"
    }
    
    else if (reduction_method == "ICA") {
      # FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
      # FM <- FM[!is.na(row.names(FM)), ]
      
      if (verbose)
        message("Reducing to independent components")
      init_ICA <- ica_helper(Matrix::t(FM), max_components,
                             use_irlba = TRUE, ...)
      x_pca <- Matrix::t(Matrix::t(FM) %*% init_ICA$K)
      W <- Matrix::t(init_ICA$W)
      weights <- W
      A <- Matrix::t(solve(weights) %*% Matrix::t(init_ICA$K))
      colnames(A) <- colnames(weights)
      rownames(A) <- rownames(FM)
      S <- weights %*% x_pca
      rownames(S) <- colnames(weights)
      colnames(S) <- colnames(FM)
      monocle:::reducedDimW(cds) <- as.matrix(W)
      monocle:::reducedDimA(cds) <- as.matrix(A)
      monocle:::reducedDimS(cds) <- as.matrix(S)
      monocle:::reducedDimK(cds) <- as.matrix(init_ICA$K)
      adjusted_S <- Matrix::t(reducedDimS(cds))
      dp <- as.matrix(dist(adjusted_S))
      cellPairwiseDistances(cds) <- dp
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      minSpanningTree(cds) <- dp_mst
      cds@dim_reduce_type <- "ICA"
    }
    else if (reduction_method == "DDRTree") {
      # FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
      # FM <- FM[!is.na(row.names(FM)), ]
      
      if (verbose)
        message("Learning principal graph with DDRTree")
      
      # TODO: DDRTree should really work with sparse matrices.
      if(auto_param_selection & ncol(cds) >= 100){
        if("ncenter" %in% names(extra_arguments)) #avoid overwrite the ncenter parameter
          ncenter <- extra_arguments$ncenter
        else
          ncenter <- cal_ncenter(ncol(FM))
        #add other parameters...
        ddr_args <- c(list(X=FM, dimensions=max_components, ncenter=ncenter, verbose = verbose),
                      extra_arguments[names(extra_arguments) %in% c("initial_method", "maxIter", "sigma", "lambda", "param.gamma", "tol")])
        #browser()
        ddrtree_res <- do.call(DDRTree, ddr_args)
      } else{
        ddrtree_res <- DDRTree(FM, max_components, verbose = verbose, ...)
      }
      if(ncol(ddrtree_res$Y) == ncol(cds))
        colnames(ddrtree_res$Y) <- colnames(FM) #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      else
        colnames(ddrtree_res$Y) <- paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      
      colnames(ddrtree_res$Z) <- colnames(FM)
      monocle:::reducedDimW(cds) <- ddrtree_res$W
      monocle:::reducedDimS(cds) <- ddrtree_res$Z
      monocle:::reducedDimK(cds) <- ddrtree_res$Y
      cds@auxOrderingData[["DDRTree"]]$objective_vals <- ddrtree_res$objective_vals
      
      adjusted_K <- Matrix::t(reducedDimK(cds))
      dp <- as.matrix(dist(adjusted_K))
      cellPairwiseDistances(cds) <- dp
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      minSpanningTree(cds) <- dp_mst
      cds@dim_reduce_type <- "DDRTree"
      cds <- monocle:::findNearestPointOnMST(cds)
    }
    else if(reduction_method == "Principal.curve") {
      dm <- destiny::DiffusionMap(t(FM))
      
      if (verbose)
        message("Learning principal graph with Principal.curve")
      
      diam_pc_tmp <- princurve::principal.curve(dm@eigenvectors[, 1:max_components]) 
      diam_pc <- as.data.frame(diam_pc_tmp$s)
      
      # FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
      # FM <- FM[!is.na(row.names(FM)), ]
      
      # colnames(ddrtree_res$Y) <- colnames(FM) #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      
      monocle:::reducedDimW(cds) <- dm@eigenvectors[, 1:2]
      diam_pc <- t(diam_pc)
      colnames(diam_pc) <- colnames(cds)
      
      monocle:::reducedDimS(cds) <- diam_pc
      monocle:::reducedDimK(cds) <- diam_pc
      cds@auxOrderingData[["DDRTree"]]$pc_res <- diam_pc
      
      adjusted_K <- Matrix::t(reducedDimK(cds))
      dp <- as.matrix(dist(adjusted_K))
      cellPairwiseDistances(cds) <- dp
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      minSpanningTree(cds) <- dp_mst
      cds@dim_reduce_type <- "DDRTree"
      cds <- monocle:::findNearestPointOnMST(cds)
    }   
    else if(reduction_method == "SimplePPT") {
      if("initial_method" %in% names(extra_arguments)){ #need to check whether or not the output match what we want
        tryCatch({
          reduced_dim_res <- extra_arguments$initial_method(FM) #variance_explained
          reduced_dim_res
        }, error = function(e) {
          error('Your initial method throws numerical errors!')
        })
      }
      else{
        if(verbose)
          message('running PCA (no further scaling or center) ...')
        
        dm <- destiny::DiffusionMap(t(FM))
        reduced_dim_res <- dm@eigenvectors[, 1:max_components] 
      }
      if(dim(reduced_dim_res)[1] != ncol(FM) & dim(reduced_dim_res)[2] < max_components )
        error("Your initial method don't generate result match the required dimension nrow(FM) * > max_components")
      
      if(verbose)
        message('running SimplePPT ...')
      
      simplePPT_args <- c(list(X=t(reduced_dim_res[, 1:max_components]), verbose = verbose),
                          extra_arguments[names(extra_arguments) %in% c("lambda", "bandwidth", "maxIter")])
      #browser()
      simplePPT_res <- do.call(simplePPT::principal_tree, simplePPT_args)
      
      colnames(simplePPT_res$MU) <- colnames(FM) #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      DCs <- t(reduced_dim_res[, 1:max_components])
      colnames(DCs) <- colnames(FM) #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      
      monocle:::reducedDimW(cds) <- DCs
      monocle:::reducedDimS(cds) <- simplePPT_res$MU
      monocle:::reducedDimK(cds) <- simplePPT_res$MU
      cds@auxOrderingData[["DDRTree"]]$objective_vals <- tail(simplePPT_res$history$objs, 1)
      
      adjusted_K <- Matrix::t(reducedDimK(cds))
      dp <- as.matrix(dist(adjusted_K))
      cellPairwiseDistances(cds) <- dp
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      minSpanningTree(cds) <- dp_mst
      cds@dim_reduce_type <- "DDRTree"
      cds <- monocle:::findNearestPointOnMST(cds)
    }
    else if(reduction_method == "L1-graph") {
      if("initial_method" %in% names(extra_arguments)){ #need to check whether or not the output match what we want
        tryCatch({
          reduced_dim_res <- extra_arguments$initial_method(FM) #variance_explained
          reduced_dim_res
        }, error = function(e) {
          error('Your initial method throws numerical errors!')
        })
      }
      else{
        if(verbose)
          message('running PCA (no further scaling or center) ...')
        dm <- destiny::DiffusionMap(t(FM))
        reduced_dim_res <- dm@eigenvectors[, 1:max_components] 
      }
      if(dim(reduced_dim_res)[1] != ncol(FM) & dim(reduced_dim_res)[2] < max_components )
        error("Your initial method don't generate result match the required dimension nrow(FM) * > max_components")
      
      if(verbose)
        message('running L1-graph ...')
      
      X <- t(reduced_dim_res[, 1:max_components])
      # D <- nrow(X); N <- ncol(X)
      # Z <- X
      
      if('C0' %in% names(extra_arguments)){
        C0 <- extra_arguments$C0
      }
      else
        C0 <- X
      Nz <- ncol(C0)
      
      # print(extra_arguments)
      if('nn' %in% names(extra_arguments))
        G <- get_knn(C0, K = extra_arguments$nn)
      else
        G <- get_knn(C0, K = 5)
      
      l1graph_args <- c(list(X = t(reduced_dim_res[, 1:max_components]), C0 = C0, G = G$G, gstruct = 'l1-graph', verbose = verbose),
                        extra_arguments[names(extra_arguments) %in% c('maxiter', 'eps', 'lambda', 'gamma', 'sigma', 'nn')])
      
      l1_graph_res <- do.call(principal_graph, l1graph_args)
      
      colnames(l1_graph_res$C) <- colnames(FM)[1:ncol(l1_graph_res$C)] #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      DCs <- t(reduced_dim_res[, 1:max_components])
      colnames(DCs) <- colnames(FM) #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      
      colnames(l1_graph_res$W) <- colnames(FM)[1:ncol(l1_graph_res$C)] #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      rownames(l1_graph_res$W) <- colnames(FM)[1:ncol(l1_graph_res$C)] #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      
      
      # row.names(l1_graph_res$X) <- colnames(cds)
      monocle:::reducedDimW(cds) <- l1_graph_res$W
      monocle:::reducedDimS(cds) <- DCs
      monocle:::reducedDimK(cds) <- l1_graph_res$C
      cds@auxOrderingData[["DDRTree"]]$objective_vals <- tail(l1_graph_res$objs, 1)
      cds@auxOrderingData[["DDRTree"]]$W <- l1_graph_res$W
      cds@auxOrderingData[["DDRTree"]]$P <- l1_graph_res$P
      
      adjusted_K <- Matrix::t(reducedDimK(cds))
      dp <- as.matrix(dist(adjusted_K))
      cellPairwiseDistances(cds) <- dp
      
      W <- l1_graph_res$W
      dimnames(l1_graph_res$W) <- list(paste('cell_', 1:nrow(W), sep = ''), paste('cell_', 1:nrow(W), sep = ''))
      W[W < 1e-5] <- 0
      gp <- graph.adjacency(W, mode = "undirected", weighted = TRUE)
      # dp_mst <- minimum.spanning.tree(gp)
      minSpanningTree(cds) <- gp
      cds@dim_reduce_type <- "DDRTree"
      cds <- monocle:::findNearestPointOnMST(cds)
    }
    else if(reduction_method == "SGL-tree") {
      if("initial_method" %in% names(extra_arguments)){ #need to check whether or not the output match what we want
        tryCatch({
          reduced_dim_res <- extra_arguments$initial_method(FM) #variance_explained
          reduced_dim_res
        }, error = function(e) {
          error('Your initial method throws numerical errors!')
        })
      }
      else{
        if(verbose)
          message('running PCA (no further scaling or center) ...')
        dm <- destiny::DiffusionMap(t(FM))
        reduced_dim_res <- dm@eigenvectors[, 1:max_components] 
      }
      if(dim(reduced_dim_res)[1] != ncol(FM) & dim(reduced_dim_res)[2] < max_components )
        error("Your initial method don't generate result match the required dimension nrow(FM) * > max_components")
      
      if(verbose)
        message('running SGL-tree ...')
      
      X <- t(reduced_dim_res[, 1:max_components])
      D <- nrow(X); N <- ncol(X)
      Z <- X
      C0 <- Z
      Nz <- ncol(C0)
      
      if('K' %in% names(extra_arguments))
        G <- get_knn(C0, K = extra_arguments$K)
      else
        G <- get_knn(C0, K = 5)
      
      l1graph_args <- c(list(X = t(reduced_dim_res[, 1:max_components]), C0 = C0, G = G$G, gstruct = 'span-tree', verbose = verbose),
                        extra_arguments[names(extra_arguments) %in% c('maxiter', 'eps', 'lambda', 'gamma', 'sigma', 'nn')])
      
      l1_graph_res <- do.call(principal_graph, l1graph_args)
      
      colnames(l1_graph_res$C) <- colnames(FM) #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      DCs <- t(reduced_dim_res[, 1:max_components])
      colnames(DCs) <- colnames(FM) #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      
      monocle:::reducedDimW(cds) <- DCs
      monocle:::reducedDimS(cds) <- DCs #1_graph_res$X
      monocle:::reducedDimK(cds) <- l1_graph_res$C
      cds@auxOrderingData[["DDRTree"]]$objective_vals <- tail(l1_graph_res$objs, 1)
      cds@auxOrderingData[["DDRTree"]]$W <- l1_graph_res$W
      cds@auxOrderingData[["DDRTree"]]$P <- l1_graph_res$P
      
      adjusted_K <- Matrix::t(reducedDimK(cds))
      dp <- as.matrix(dist(adjusted_K))
      cellPairwiseDistances(cds) <- dp
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      minSpanningTree(cds) <- dp_mst
      cds@dim_reduce_type <- "DDRTree"
      cds <- monocle:::findNearestPointOnMST(cds)
    }
    else if(reduction_method == "spl"){
      message('This option is not ready yet')
    }
    else {
      stop("Error: unrecognized dimensionality reduction method")
    }
  }
  cds
}

#' Sparsify the direct network with biologically inspired constraints  
#' 
#' This function implements the direct network sparsifier algorithm described in the supplementary file of Scribe manuscript. 
#' It relies on the lpSolveAPI package to solve the linear programming problem.  
#' 
#' @param W A matrix which contains the causality score between genes. This matrix is asymmetric. 
#' @param alphae The exponent of the power law distribution.
#' @param lambda The parameter of the exponential distribution.
#' @param B The budget parameter introduced to prevent trivial solution.
#' @return A new asymmeric matrix which contains a sparse direct network structure estimated based on the network sparsifier.
#' @examples
#' \dontrun{
#' sparse_W <- dn_sparsifer(W)
#' }
#' @export
dn_sparsifer <- function(W, alpha = 2, lambda = 20, B = 184) {
  N <- nrow(W)
  
  tmp <- matrix(-1, nrow = N, ncol = N) + diag(1, N)
  A <- matrix(as.vector(tmp), nrow = 1)
  b <- -B
  
  # initialize Theta 
  Theta <- matrix(0, nrow = N, ncol = N)
  w <- matrix(unlist(W), nrow = 1)
  sortw <- sort(w, decreasing = T)
  sigma <- sortw[B]
  Theta[W > sigma] <- 1
  
  idx <- which(sortw < 1e-10)
  xi <- sortw[idx[1] - 1] * matrix(1, nrow = N, ncol = 1)
  
  iter <- 1 
  
  while(1) {
    row_norm <- matrix(rowSums(Theta), ncol = 1)
    col_norm <- matrix(colSums(Theta), nrow = 1)
    
    eta <- repmat(1 / (row_norm + xi), 1, N) + repmat(1 / (col_norm + t(xi)), N, 1)
    
    # termination condition 
    obj1 <- -sum(sum(W * (Theta - lambda * tmp)))
    obj2 <- alpha * sum(log ((row_norm + xi)))
    obj <- obj1 + obj2 
    
    if(iter > 1) {
      obj_diff <- old_obj - obj 
      rel_diff <- abs(obj_diff / old_obj)
      message('iter = ', iter, ' obj = ', obj, ' obj_diff = ', obj_diff, ' rel_diff = ', rel_diff)
      if(rel_diff < 1e-6) {
        message('converge!')
        break 
      }
    }
    old_obj <- obj 
    
    ######################################################################################################################################################### 
    # solve linear reprogramming 
    ######################################################################################################################################################### 
    fmat <- alpha * eta - lambda * tmp - Theta
    f <- matrix(unlist(fmat), nrow = 1)
    lprec <- make.lp(length(b), length(f))
    set.objfn(lprec, f)
    for(i in 1:nrow(A)) {
      add.constraint(lprec, A[i, ], "<=", b)
    }
    
    set.bounds(lprec, lower = rep(0, N * N), upper = rep(1, N*N)) #				set.bounds(lprec, lower = c(rep(0, nw), -Inf*rep(1, K*D)), columns = 1:length(b))
    #solve the model, if this return 0 an optimal solution is found
    solve(lprec)
    
    obj_W <- get.objective(lprec)
    theta <- get.variables(lprec)
    
    Theta <- matrix(theta, nrow = N, ncol = N)
    Theta[Theta < 1e-6] <- 0 
    
    iter <- iter + 1
  }
  
  return(Theta)
}
