#' Calculate turning point for a linear or branched trajectory 
#'
#' This function estimates the inflection point or the gene-wise branch point for each gene in a cell trajectory, without or with branch points respectively
#' @param cds_subset a cds object after trajectory reconstruction 
#' @return a updated cds with a newly added column (turning_point) in pData indicates the inflection or branch time point.
#' to do: write a function to integrate BGP with Monocle 
estimate_turning_point <- function(cds_subset, type = c('order', 'value'), branch_point = 1, cores = 1) {
  cds_exprs <- pData(cds_subset)
  cds_exprs$turning_point <- 0
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
          a <- which.min(abs(sort(pData(b_cds_subset)[pData(b_cds_subset)$Branch == branches[1], "Pseudotime"]) - abs(x)))
          b <- which.min(abs(sort(pData(b_cds_subset)[pData(b_cds_subset)$Branch == branches[2], "Pseudotime"]) - abs(x)))
          
          round(mean(c(a, b), na.rm = T))
        }
      )) # get the order of the cell corresponding the inflection point)
      fData(cds_subset)[, 'turning_point'] <- iplasts
    }
    # print(iplasts)
  } else{
    rng <- range(pData(cds_subset)$Pseudotime)
    new_data <- data.frame(Pseudotime = seq(rng[1], rng[2], length.out = 100))
    cds_exprs <- genSmoothCurves(cds_subset, cores = cores, trend_formula = "~sm.ns(Pseudotime, df = 3)",
                    relative_expr = TRUE, new_data = new_data)
    for(i in 1:nrow(cds_exprs)) {
      data <- cds_exprs[i, ]; 
      inflection_point <- bede(1:length(data), data, 0)
      if(!is.finite(inflection_point$iplast))
        inflection_point <- bede(1:length(data), data, 1)
    }
    
    message('turning_point is ', inflection_point$iplast)
    iplast  <- inflection_point$iplast
    
    if(type == 'value')  { 
      fData(cds_subset)[i, 'turning_point'] <- iplast # get the order of the cell corresponding the inflection point 
    } else if(type == 'order') {
     fData(cds_subset)[i, 'turning_point'] <- which.min(abs(sort(pData(cds_subset)$Pseudotime) - iplast)) # get the order of the cell corresponding the inflection point 
    }
  }
  
  return(cds_subset)
}

# lung_update <- estimate_turning_point(lung)
# lung_update <- estimate_turning_point(lung, type = 'value')

#' Calculate conditionally RDI value
#'
#' This function estimates the RDI value for all gene-pair combination from the genes_data (or that specified in the super_graph) 
#' in the pseudotime/time series data, conditioned on the top top_incoming_k incoming nodes. 
#' The result from the calculate_rdi function will be used to identify the proper delay for the gene-pair under test corresponds to that 
#' with highest RDI as well as the proper delays for the k incoming nodes which corresponds to that with highest RDI values.
#' @param genes_data Time / Pseudotime series expression data (Rows are samples while columns are features).  
#' @param super_graph A graph including all possible interactions used for performing the causality inference. When it is NULL, 
#' all possible pairs of interactions for the genes_data are used, otherwise only the interactions specified in the graph are used.
#' Note that the super_graph only accept integer matrix for now (each integer corresponds to a particular gene in the genes_data).  
#' @param genes_data pseudo-time/time series for the gene expression data
#' @param rdi_list a list returned from calculate_rdi function. 
#' @param top_incoming_k The number of genes to be conditioned on when calculating the conditional RDI values 
#' @return a dataframe storing conditional RDI results. First two columns are the id names for the genes.
#' Third column is the conditional RDI value.
#' 
calculate_conditioned_rdi <- function(genes_data, super_graph = NULL, rdi_list, top_incoming_k = 1) {
  if(!all(is.finite(genes_data))) {
    warning('your data includes non finite values. The associated genes are removed!')
    finite_genes <- apply(genes_data, 2, function(x) all(is.finite(x)))
    genes_data <- genes_data[, finite_genes] # remove non-infinite gene 
  }
  
  n_genes <- ncol(genes_data); 
  n_samples <- nrow(genes_data); 
  set.seed(2017)
  noise = matrix(rnorm(mean = 0, sd = 1e-10, n_genes * n_samples), nrow = n_samples)
  
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
  cRDI_mat <- calculate_conditioned_rdi_cpp_wrap(genes_data + noise, as.matrix(super_graph),
                                            max_rdi_value, max_rdi_delays, top_incoming_k)

  dimnames(cRDI_mat) <- list(colnames(genes_data), colnames(genes_data))
  #cRDI <- calculate_conditioned_rdi_cpp_wrap(t(lung_exprs_AT1[1:30, 1:50]) + noise, as.matrix(super_graph), max_rdi_value, max_rdi_delays, 2L)
  
  return(cRDI_mat); 
}

#' @title
#' calculate_rdi
#' @description
#' This subroutine calculates RDI of Pseudo-time series expression data 
#' 
#' @param Time / Pseudo-time series expression data where each row is a sample, each column is a feature 
#' 
#' @param delays Time lags used to estimate the RDI values 
#' 
#' @param method either 1, 2 represents either using (RDI: restricted direction information) or LMI (lagged mutual information) 
#'
#' @details
#' \code{calculate_rdi} takes a time / Pseudo-time series expression data as well as the time lags and then calculate the restricted direct information between each pair of genes under different delays. 
#' @return a list with four components: matrix for RDI (dimension is number of genes X length of delays times number of genes), vector of delays, max_rdi_value and max_rdi_delays 
#' @export
#' 
calculate_rdi <- function(genes_data, delays, super_graph = NULL, turning_points = NULL, method = 1) {
  if(!all(is.finite(genes_data))) {
    warning('your data includes non finite values. The associated genes are removed!')
    finite_genes <- apply(genes_data, 2, function(x) all(is.finite(x)))
    genes_data <- genes_data[, finite_genes] # remove non-infinite gene 
  }
  
  n_genes <- ncol(genes_data); 
  n_samples <- nrow(genes_data); 
  set.seed(2017)
  noise = matrix(rnorm(mean = 0, sd = 1e-10, n_genes * n_samples), nrow = n_samples)
  
  if(is.null(super_graph)) 
  {
    tmp <- expand.grid(1:n_genes, 1:n_genes, stringsAsFactors = F)
    super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to c++ index
  }
  else
  {
    if(class(super_graph[1, 1]) != "character"){
      if(max(super_graph) >  n_samples - 1 | min(super_graph) < 0)
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

  RDI_list <- calculate_rdi_cpp_wrap(as.matrix(genes_data + noise), delays, as.matrix(super_graph), turning_points, method)
  
  # assign the name for each dimnames
  dimnames(RDI_list$RDI) <- list(colnames(genes_data), rep(colnames(genes_data), length(delays))) 
  dimnames(RDI_list$max_rdi_value) <- list(colnames(genes_data), colnames(genes_data))
  dimnames(RDI_list$max_rdi_delays) <- list(colnames(genes_data), colnames(genes_data)) 

  return(RDI_list); 
}

####################################################################################################################################################################################
# update this to multiple run 
####################################################################################################################################################################################

#' Calculate conditionally RDI value
#'
#' This function estimates the RDI value for all gene-pair combination from the genes_data (or that specified in the super_graph) 
#' in the pseudotime/time series data, conditioned on the top top_incoming_k incoming nodes. 
#' The result from the calculate_rdi function will be used to identify the proper delay for the gene-pair under test corresponds to that 
#' with highest RDI as well as the proper delays for the k incoming nodes which corresponds to that with highest RDI values.
#' @param genes_data Time / Pseudotime series expression data (Rows are samples while columns are features).  
#' @param run_vec A vector describes which run (lineage) does the current cell come from. It has the some length as the row length of genes_data. 
#' @param super_graph A graph including all possible interactions used for performing the causality inference. When it is NULL, 
#' all possible pairs of interactions for the genes_data are used, otherwise only the interactions specified in the graph are used.
#' Note that the super_graph only accept integer matrix for now (each integer corresponds to a particular gene in the genes_data).  
#' @param rdi_list a list returned from calculate_rdi function. 
#' @param top_incoming_k The number of genes to be conditioned on when calculating the conditional RDI values 
#' @return a dataframe storing conditional RDI results. First two columns are the id names for the genes.
#' Third column is the conditional RDI value.
#' 
calculate_conditioned_rdi_multiple_run <- function(genes_data, run_vec = NULL, super_graph = NULL, rdi_list, top_incoming_k = 1) {
  if(!all(is.finite(genes_data))) {
    warning('your data includes non finite values. The associated genes are removed!')
    finite_genes <- apply(genes_data, 2, function(x) all(is.finite(x)))
    genes_data <- genes_data[, finite_genes] # remove non-infinite gene 
  }
  
  if(is.null(run_vec))
    run_vec <- rep(0, nrow(genes_data)) 
  if(min(run_vec) != 0 | any(diff(unique(run_vec)) != 1) | !is.vector(run_vec))
    stop("The run_vec starts from 0 to the total number of runs")
  
  n_genes <- ncol(genes_data); 
  n_samples <- nrow(genes_data); 
  set.seed(2017)
  noise = matrix(rnorm(mean = 0, sd = 1e-10, n_genes * n_samples), nrow = n_samples)
  
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
  cRDI_mat <- calculate_multiple_run_conditioned_rdi_wrap(genes_data + noise, as.matrix(super_graph), 
                                                 max_rdi_value, max_rdi_delays, run_vec, top_incoming_k)
  
  dimnames(cRDI_mat) <- list(colnames(genes_data), colnames(genes_data))
  #cRDI <- calculate_conditioned_rdi_cpp_wrap(t(lung_exprs_AT1[1:30, 1:50]) + noise, as.matrix(super_graph), max_rdi_value, max_rdi_delays, 2L)
  
  return(cRDI_mat); 
}

#' @title
#' calculate_rdi
#' @description
#' This subroutine calculates RDI of Pseudo-time series expression data 
#' 
#' @param genes_data Time / Pseudotime series expression data (Rows are samples while columns are features).  
#' @param run_vec A vector describes which run (lineage) does the current cell come from. It has the some length as the row length of genes_data. 
#' @param delays Time lags used to estimate the RDI values (be ready to update this as a matrix or a vector with the same length as the run_vec once we have a way to directly calculate the proper time lag need to be used)
#' @param super_graph A graph including all possible interactions used for performing the causality inference. When it is NULL, 
#' all possible pairs of interactions for the genes_data are used, otherwise only the interactions specified in the graph are used.
#' Note that the super_graph only accept integer matrix for now (each integer corresponds to a particular gene in the genes_data).  
#' @param method either 1, 2 represents either using (RDI: restricted direction information) or LMI (lagged mutual information) 
#'
#' @details
#' \code{calculate_rdi} takes a time / Pseudo-time series expression data as well as the time lags and then calculate the restricted direct information between each pair of genes under different delays. 
#' @return a list with four components: matrix for RDI (dimension is number of genes X length of delays times number of genes), vector of delays, max_rdi_value and max_rdi_delays 
#' @export
#' 
calculate_rdi_multiple_run <- function(genes_data, run_vec = NULL, delays, super_graph = NULL, method = 1) {
  if(!all(is.finite(genes_data))) {
    warning('your data includes non finite values. The associated genes are removed!')
    finite_genes <- apply(genes_data, 2, function(x) all(is.finite(x)))
    genes_data <- genes_data[, finite_genes] # remove non-infinite gene 
  }
  
  if(is.null(run_vec))
    run_vec <- rep(0, nrow(genes_data)) 
  if(min(run_vec) != 0 | any(diff(unique(run_vec)) != 1) | !is.vector(run_vec))
    stop("The run_vec starts from 0 to the total number of runs")
  
  n_genes <- ncol(genes_data); 
  n_samples <- nrow(genes_data); 
  set.seed(2017)
  noise = matrix(rnorm(mean = 0, sd = 1e-10, n_genes * n_samples), nrow = n_samples)
  
  if(is.null(super_graph)) 
  {
    tmp <- expand.grid(1:n_genes, 1:n_genes, stringsAsFactors = F)
    super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to c++ index
  }
  else
  {
    if(class(super_graph[1, 1]) != "character"){
      if(max(super_graph) >  n_samples - 1 | min(super_graph) < 0)
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
  
  RDI_list <- calculate_rdi_multiple_run_cpp_wrap(as.matrix(genes_data + noise), delays, run_vec, as.matrix(super_graph), method)
  
  # assign the name for each dimnames
  dimnames(RDI_list$RDI) <- list(colnames(genes_data), rep(colnames(genes_data), length(delays))) 
  dimnames(RDI_list$max_rdi_value) <- list(colnames(genes_data), colnames(genes_data))
  dimnames(RDI_list$max_rdi_delays) <- list(colnames(genes_data), colnames(genes_data)) 
  
  return(RDI_list); 
}

# think about better ways to deal with complex tree structure ()

#' Create a pseudo-time-seires
#'
#' This function takes a cds (branched or not) and then convert to a pseudotime serious 
#' @param cds A cds which has been ordered with Monocle  
#' @param branch_points Vector for the branch points. If it is null, there is no branching in the data 
#' @return a list storing the pseudo-time-series data and a vector storing the run id for each cell (row) in the data 
#' 
# return(list(data = data, run_vec = run_vec))
createPTS <- function(cds, branch_points = NULL) {
  if(is.null(pData(cds)$Pseudotime))
    stop('Please order you cell dataset with Monocle before running createPTS')
  
  if(!is.integer(branch_points) | any(branch_points > length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points)))
    stop("Branch point should be positive integer and should not be larger than the possible number of branch points. Check branch points with plot_cell_trajectory function in Monocle")
    
  if(is.null(branch_points) == FALSE) { # enable multiple branch points
    data_tmp <- c()
    run_vec_tmp <- c()
    pseudotime <- c()
    run_ind <- 0
    
    for(i in branch_points) {
      tmp <- buildBranchCellDataSet(cds, progenitor_method = 'duplicate', branch_point = i)
      level_vec <- levels(pData(tmp)$Branch)
      
      data_tmp <- rbind(data_tmp, t(exprs(tmp)))
      Branch <- as.character(pData(tmp)$Branch)
      Branch[Branch == level_vec[1]] <- run_ind
      Branch[Branch == level_vec[2]] <- run_ind + 1
      
      run_vec_tmp <- c(run_vec_tmp, Branch)
      pseudotime <- c(pseudotime, pData(tmp)$Pseudotime)
      run_ind <- run_ind + 2
    }
    
    # run 0, 1, ...
    run_vec_tmp <- as.numeric(run_vec_tmp) 
    run_vec <- run_vec_tmp[order(run_vec_tmp)]
    data_tmp <- data_tmp[order(run_vec_tmp), ]
    data <- data_tmp[order(run_vec_tmp), ]
    
    for(i in unique(run_vec_tmp)) {
        cells_in_current_run <- which(run_vec_tmp == i)
        order_cells_in_current_run_by_pseudotime <- order(pseudotime[cells_in_current_run])
        
        data[cells_in_current_run, ] <- data_tmp[cells_in_current_run[order_cells_in_current_run_by_pseudotime], ]
    }
  }
  else { # no branching in the data 
    data <- t(exprs(cds))
    data <- data[order(pData(cds)$Pseudotime), ] # order by pseudotime 
    run_vec <- rep(0, ncol(cds))
  }
  
  return(list(data = data, run_vec = run_vec))
}

# 
rdi_crdi_pseudotime <- function(data, window_size = 40, delay = 1) {
  win_range <- nrow(data) - window_size - 1
  gene_num <- ncol(data)
  
  window_gene_gene_result <- array(dim = c(win_range + 1, gene_num, gene_num)) 
  crdi_window_gene_gene_result <- window_gene_gene_result 
  
  run_vec <- rep(1, nrow(data))
  tmp <- expand.grid(1:ncol(data), 1:ncol(data), stringsAsFactors = F)
  super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index
  super_graph <- super_graph[, c(2, 1)]
  
  # .Call('_Scribe_calculate_rdi_multiple_run_cpp', PACKAGE = 'Scribe', expr_data, delays, run_vec, super_graph, turning_points, method)
  
  for(i in 1:(win_range + 1)) {
    message('current i is ', i)
    rdi_list <- calculate_rdi_multiple_run_cpp(expr_data = data[i:(i + window_size), ], delay = c(5), run_vec = run_vec[i:(i + window_size)] - 1, super_graph = as.matrix(super_graph), turning_points = 0, method = 1) #* 100 + noise
    con_rdi_res_test <- calculate_multiple_run_conditioned_rdi_wrap(data[i:(i + window_size), ], as.matrix(super_graph), as.matrix(rdi_list$max_rdi_value), as.matrix(rdi_list$max_rdi_delays), run_vec[i:(i + window_size)] - 1, 1)
    
    window_gene_gene_result[i, , ] <- rdi_list$max_rdi_value
    crdi_window_gene_gene_result[i, , ] <- con_rdi_res_test
    
  }
  
  return(list(rdi_res = window_gene_gene_result, crdi_res = crdi_window_gene_gene_result))
}

# this following function used to re-order the trajectory using principal curve 
reduceDimension <- function(cds,
                            max_components=2,
                            reduction_method=c("DDRTree", "ICA", 'tSNE', "Principal.curve", "SimplePPT", 'L1-graph', 'SGL-tree'),
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
  
  FM <- normalize_expr_data(cds, norm_method, pseudo_expr)
  
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
    reducedDimW(cds) <- as.matrix(reducedDim)
    reducedDimA(cds) <- as.matrix(reducedDim)
    reducedDimS(cds) <- as.matrix(reducedDim)
    reducedDimK(cds) <- as.matrix(reducedDim)
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
      reducedDimW(cds) <- as.matrix(W)
      reducedDimA(cds) <- as.matrix(A)
      reducedDimS(cds) <- as.matrix(S)
      reducedDimK(cds) <- as.matrix(init_ICA$K)
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
      reducedDimW(cds) <- ddrtree_res$W
      reducedDimS(cds) <- ddrtree_res$Z
      reducedDimK(cds) <- ddrtree_res$Y
      cds@auxOrderingData[["DDRTree"]]$objective_vals <- ddrtree_res$objective_vals
      
      adjusted_K <- Matrix::t(reducedDimK(cds))
      dp <- as.matrix(dist(adjusted_K))
      cellPairwiseDistances(cds) <- dp
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      minSpanningTree(cds) <- dp_mst
      cds@dim_reduce_type <- "DDRTree"
      cds <- findNearestPointOnMST(cds)
    }
    else if(reduction_method == "Principal.curve") {
      dm <- destiny::DiffusionMap(t(FM))
      
      if (verbose)
        message("Learning principal graph with Principal.curve")
      
      diam_pc <- princurve::principal.curve(dm@eigenvectors[, 1:2]) 
      diam_pc <- as.data.frame(diam_pc$s)
      
      # FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
      # FM <- FM[!is.na(row.names(FM)), ]
      
      colnames(ddrtree_res$Y) <- colnames(FM) #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      
      reducedDimW(cds) <- dm@eigenvectors[, 1:2]
      diam_pc$s <- t(diam_pc$s)
      colnames(diam_pc$s) <- colnames(cds)
      
      reducedDimS(cds) <- diam_pc$s
      reducedDimK(cds) <- diam_pc$s
      cds@auxOrderingData[["Principal.curve"]]$pc_res <- diam_pc
      
      adjusted_K <- Matrix::t(reducedDimK(cds))
      dp <- as.matrix(dist(adjusted_K))
      cellPairwiseDistances(cds) <- dp
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      minSpanningTree(cds) <- dp_mst
      cds@dim_reduce_type <- "Principal.curve"
      cds <- findNearestPointOnMST(cds)
    }
    else {
      stop("Error: unrecognized dimensionality reduction method")
    }
  }
  cds
}
