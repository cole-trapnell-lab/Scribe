#' Causal network inference limitted to the interaction between informative transcription factors and those to other informative target genes. 
#' 
#' This function implements a procedure to build a causal network based on some prior information of the transcription factors (TFs). It calculates the information 
#' transfer from TFs to the targets as well as the information transfer between transcription factors while avoiding the calculation of information transfer between 
#' the targets and from targets to the TFs. This function accepts a vector for the gene short names for all the transcription factors and another vector for the 
#' informative genes selected through BEAM or methods.
#' 
#' @param cds CellDataSet for the experiment
#' @param TF A vector of the gene short names for all the transcription factors. 
#' @param informative_genes A vector of the informative genes used for network inference, which is identified through BEAM or other methods. It should include some TFs and normally other target genes.  
#' @param delays A vector of time delays used during information transfer estimation between genes 
#' @param include_conditioning A logic argument to determine whether or not we should also calculate the conditioning RDI values 
#' @param smoothing A logic argument to determine whether or not we should also smooth the data using loess before we calculate the RDI values 
#' @param cluster_TFs A logic value to determine whether or not we should cluster the TFs' expression (for either the group of TFs or targets) before estimating the RDI/cRDI value. 
#' @param cluster_targets A logic value to determine whether or not we should cluster the target genes' expression (for either the group of TFs or targets) before estimating the RDI/cRDI value. 
#' @param cluster_TFs_num Number of TF clusters you will use for estimating the information transfer.  
#' @param cluster_targets_num Number of target clusters you will use for estimating the information transfer.  
#' @param cluster_num_method Which method you will use to determine the proper number of cluster for either the TFs or the targets. 
#' Options include NULL (default, don't use any method), mcclust (model-based optimal number of clusters: TF) and finally pamk (number of clusters estimated by optimum average silhouette width). 
#' @param hclust_method Which cluster method in hclust you will use. Default to be ward.D2.
#' @param norm_method Determines how to transform expression values prior to estimate information transfer
#' @param scale Determine whether or not to scale data prior to estimate information transfer
#' @param scale_max The maximum value (in standard deviations) to show in the heatmap. Values larger than this are set to the max.
#' @param scale_min The minimum value (in standard deviations) to show in the heatmap. Values smaller than this are set to the min.
#' @return a matrix storing the RDI values for interactions from TFs to the targets as well as between transcription factors.
#' @importFrom Mclust mclust
#' @importFrom pamk fpc
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' TF <- c('H19', 'Ccnd2', 'Ccnd2', 'Scnn1g'); TF_id <- row.names(subset(fData(lung), gene_short_name %in% TF))
#' informative_genes_id <- row.names(lung)
#' wired(lung, TF_id, informative_genes_id)
#' }
#' @export
#' 
wired <- function(cds, TF, informative_genes, 
                    delays = c(1, 10, 15), 
                    include_conditioning = FALSE,
                    smoothing = FALSE, 
                    cluster_TFs = FALSE,
                    cluster_targets = FALSE, 
                    cluster_TFs_num = NULL,
                    cluster_targets_num = NULL, 
                    cluster_num_method = c(NULL, "mcclust", "pamk"),
                    hclust_method = "ward.D2",
                    norm_method = c('log', 'vstExprs'),
                    scale = FALSE,
                    scale_max = 3, 
                    scale_min = -3,
                    ...) {
  
  pseudocount <- 1 

  # 1. read in the TF list
  TF_vec_names <- intersect(TF, informative_genes)
  target_vec_names <- setdiff(informative_genes, TF)
  
  unique_gene <- unique(c(TF, informative_genes))
  gene_name_ids <- intersect(row.names(cds), unique_gene)

  if(length(unique_gene) != length(gene_name_ids)) {
    stop("The cds you provided doesn't include all genes from the TF and informative_genes vector!")
  }

  cds_subset <- cds[gene_name_ids, ]
  pData(cds_subset)$Pseudotime <- order(pData(cds_subset)$Pseudotime)

  exprs_data <- exprs(cds_subset)[, pData(cds_subset)$Pseudotime] # get the expression matrix 
  
  # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
  if(norm_method == 'vstExprs' && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == FALSE) {
    exprs_data = vstExprs(cds_subset, expr_matrix=exprs_data)
  }     
  else if(norm_method == 'log') {
    exprs_data = log10(exprs_data+pseudocount)
  }
  
  annotation_TF_cluster = NULL
  annotation_target_cluster = NULL
  cRDI <- NULL 

  if(cluster_TFs | cluster_targets | is.null(cluster_TFs_num) == FALSE | is.null(cluster_TFs_num) == FALSE) {
    if(smoothing) {
      for(i in 1:ncol(exprs_data)) {
        df <- data.frame(Pseudotime = 1:ncol(exprs_data), Expression = exprs_data[i, ])
        test <- loess(Expression ~ Pseudotime, df)
        exprs_data[i, ] <-  predict(test)
      }
    }

    cds_subset@assayData$exprs <- exprs_data # change the cds to smoothed exprs_data 

    # Row-center the data.
    m <- exprs_data
    m=m[!apply(m,1,sd)==0,]
    if(scale) {
      m=Matrix::t(scale(Matrix::t(m),center=TRUE))
      m=m[is.na(row.names(m)) == FALSE,]
      m[is.nan(m)] = 0
      m[m>scale_max] = scale_max
      m[m<scale_min] = scale_min
    }
    
    # TF_vec_names <- intersect(TF_vec_names, row.names(m))
    # target_vec_names <- intersect(target_vec_names, row.names(m))
    
    m_tfs <- m[TF_vec_names, ]
    m_targets <- m[target_vec_names, ]

    if(cluster_num_method == 'mcclust') {
      clust_num_check_tfs <- mclust::Mclust(t(m_tfs), G=1:min(10, length(TF_vec_names) / 2)) #(round(length(m_tfs) / 2))
      cluster_TFs_num <- dim(clust_num_check_tfs$z)[2]

      
      clust_num_check_targets <- mclust::Mclust(t(m_targets), G=1:min(10, length(target_vec_names) / 2))
      cluster_targets_num <- dim(clust_num_check_targets$z)[2]

      cat("model-based optimal number of clusters: TF: ", cluster_TFs_num, ", targets: ", cluster_targets_num, "\n")
    } else if(cluster_num_method == 'pamk') {
      dissimilarity_mat_tfs <- 1 - cor(t(m_tfs)); 
      dissimilarity_mat_targets <- 1 - cor(t(m_targets));      

      clust_num_check_tfs <- fpc::pamk(dissimilarity_mat_tfs, diss = T)
      clust_num_check_targets <- fpc::pamk(dissimilarity_mat_targets, diss = T)

      cluster_TFs_num <- clust_num_check_tfs$nc      
      cluster_targets_num <- clust_num_check_targets$nc
      
      cat("number of clusters estimated by optimum average silhouette width: TF: ", cluster_TFs_num, ", targets: ", cluster_targets_num, "\n")
    }
    # cluster TFs
    if(cluster_TFs | is.null(cluster_TFs_num) == FALSE) {
      # replace the hc clustering with k-means clustering 
      row_dist <- as.dist((1 - cor(Matrix::t(m_tfs)))/2)
      row_dist[is.na(row_dist)] <- 1
      
      m_hclust <- hclust(row_dist, method = hclust_method) # $order, ]
      
      annotation_TF_cluster <- data.frame(Cluster=factor(cutree(m_hclust, cluster_TFs_num)), row.names = row.names(m_tfs)) 

      # calculate the average of gene expression for each cluster and then run RDI between the TFs and the target.  
      m_TFs_clusters <- matrix(nrow = ncol(cds_subset), ncol = cluster_TFs_num)
      for(cluster_ind in 1:cluster_TFs_num) {
        gene_inds <- annotation_TF_cluster$Cluster == cluster_ind
        df <- data.frame(Pseudotime = 1:ncol(exprs_data), Expression = as.vector(t(exprs_data[gene_inds, ])))
        test <- loess(Expression ~ Pseudotime, df)
        
        m_TFs_clusters[, cluster_ind] <- predict(test) # colMeans(m_tfs[gene_inds, ])
      }

      colnames(m_TFs_clusters) <- paste0("TFs_", 1:cluster_TFs_num)
    } else {
      m_TFs_clusters <- m_tfs
    }
    # cluster TFs
    if(cluster_targets | is.null(cluster_targets) == FALSE) {
      row_dist <- as.dist((1 - cor(Matrix::t(m_targets)))/2)
      row_dist[is.na(row_dist)] <- 1
      
      m_hclust <- hclust(row_dist, method = hclust_method) # $order, ]
      
      annotation_target_cluster <- data.frame(Cluster=factor(cutree(m_hclust, cluster_targets_num)), row.names = row.names(m_targets)) 

      # calculate the average of gene expression for each cluster and then run RDI between the TFs and the target.  
      m_targets_clusters <- matrix(nrow = ncol(cds_subset), ncol = cluster_targets_num)
      for(cluster_ind in 1:cluster_targets_num) {
        gene_inds <- annotation_target_cluster$Cluster == cluster_ind
        df <- data.frame(Pseudotime = 1:ncol(exprs_data), Expression = as.vector(t(exprs_data[gene_inds, ])))
        test <- loess(Expression ~ Pseudotime, df)

        m_targets_clusters[, cluster_ind] <- predict(test) # colMeans(m_targets[gene_inds, ])
      }

      colnames(m_targets_clusters) <- paste0("target_", 1:cluster_targets_num)
    } else {
      m_targets_clusters <- m_targets
    }
    # 2. build the super-graph:

    TF_vec_names <- colnames(m_TFs_clusters)
    target_vec_names <- colnames(m_targets_clusters)

    exprs_data <- cbind(m_TFs_clusters, m_targets_clusters) # RDI requires row: cells; column: genes

    TF_pair <- expand.grid(TF_vec_names, TF_vec_names, stringsAsFactors = F) # between TFs 
    TF_target_pair <- expand.grid(TF_vec_names, target_vec_names, stringsAsFactors = F) # from TFs to the targets 

    tmp <- rbind(TF_pair, TF_target_pair)
    
    tmp[, 1] <- match(tmp[, 1], colnames(exprs_data)) # convert to positions 
    tmp[, 2] <- match(tmp[, 2], colnames(exprs_data))
    
    all_pairwise_gene <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index

    RDI_res <- calculate_rdi_cpp_wrap(as.matrix(exprs_data), delays = delays, super_graph = as.matrix(all_pairwise_gene), turning_points = 0, method = 1, uniformalize = F) # run RDI 
    
    if(include_conditioning) {
      cRDI_res <- calculate_conditioned_rdi_cpp_wrap(as.matrix(exprs_data), super_graph = as.matrix(all_pairwise_gene), 
                                                     max_rdi_value = RDI_res$max_rdi_value, max_rdi_delays = RDI_res$max_rdi_delays, k = 1, uniformalize = FALSE) # run cRDI 
    }
  } else {
    # 2. build the super-graph:
    cds_subset@assayData$exprs <- exprs_data # change the cds to smoothed exprs_data 
    exprs_data <- t(exprs_data) 
    
    TF_pair <- expand.grid(TF_vec_names, TF_vec_names, stringsAsFactors = F) # between TFs 
    TF_target_pair <- expand.grid(TF_vec_names, target_vec_names, stringsAsFactors = F) # from TFs to the targets 

    tmp <- rbind(TF_pair, TF_target_pair)
    
    tmp[, 1] <- match(tmp[, 1], colnames(exprs_data))
    tmp[, 2] <- match(tmp[, 2], colnames(exprs_data))
    
    all_pairwise_gene <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index

    RDI_res <- calculate_rdi(cds_subset, delays = delays, super_graph = all_pairwise_gene, log = FALSE, ...) # run RDI 
  
    if(include_conditioning) {
      cRDI_res <- calculate_conditioned_rdi(exprs_data, rdi_list = RDI_res, ...) # run cRDI 
    }
  }

  return(list(exprs_data = exprs_data, 
    RDI_res = RDI_res, 
    cRDI = cRDI,
    annotation_TF_cluster = annotation_TF_cluster, 
    annotation_target_cluster = annotation_target_cluster))
}
