# # test the comprehensive dataset and show the result :

# # read the file and then calculate the RDI values for each pair of genes 
# neuronsame_cell_two_stage_dropout_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/neuronsame_cell_two_stage_dropout_pseudotime.txt', sep = '\t')
# # calculate AUC and then make the boxplot: 
# run_ids <- unique(neuronsame_cell_two_stage_dropout_pseudotime$V2)

# for(id in run_ids) {
#   subset_exprs_mat <- as.matrix(subset(neuronsame_cell_two_stage_dropout_pseudotime, V2 == 'R145')[, -c(1:3)])
#   noise = matrix(rnorm(mean = 0, sd = 1e-10, nrow(subset_exprs_mat) * ncol(subset_exprs_mat)), nrow = nrow(subset_exprs_mat))
  
#   subset_exprs_mat_noise <- subset_exprs_mat + noise
#   subset_exprs_mat_noise <- subset_exprs_mat_noise - min(subset_exprs_mat_noise)
  
#   a <- Sys.time()
#   RDI_parallel_res <- calculate_rdi(t(subset_exprs_mat_noise), delays = c(1, 2, 3))
#   b <- Sys.time()
  
#   a <- Sys.time()
#   cRDI_parallel_res <- calculate_conditioned_rdi(t(subset_exprs_mat) + noise, rdi_list = RDI_parallel_res)
#   b <- Sys.time()
  
# }

# function for implementing Cole's code:
#' Causal network inference for between informative transcription factor and some informative target genes. 
#' 
#' This function implements a procedure to build a causal network based on some prior information of the transcription factors (TFs). It calculates the information 
#' transfer from TFs to the targets as well as the information transfer between transcription factors while avoiding the calculation of information transfer between 
#' the targets and from targets to the TFs. This function accepts a vector for the gene short names for all the transcription factors and another vector for the 
#' informative genes selected through BEAM or methods.
#' 
#' @param cds CellDataSet for the experiment
#' @param TF A vector of the gene short names for all the transcription factors. 
#' @param informative_genes A vector of the informative genes used for network inference, which is identified through BEAM or other methods
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' gene_pairs_mat <- matrix(c('H19', 'Ccnd2', 'Ccnd2', 'Scnn1g'), ncol = 2)
#' wired(lung, TF_vec_id, informative_genes)
#' }
#' @export
#' 
wired <- function(cds, TF, informative_genes) {
  # 1. read in the TF list
  TF_vec_names <- intersect(TF, informative_genes)
  target_vec_names <- setdiff(informative_genes, TF)
  
  # 2. build the super-graph:
  TF_pair <- expand.grid(TF_vec_names, TF_vec_names, stringsAsFactors = F) # between TFs 
  TF_target_pair <- expand.grid(TF_vec_names, target_vec_names, stringsAsFactors = F) # from TFs to the targets 

  # 3. run LMI (lagged mutual information)
  
  # 4. fit the data into multiple clusters and get the average expression pattern in each cluster 
  # newdata <- data.frame(Pseudotime = sort(pData(cds)$Pseudotime))
  # smooth_res <- genSmoothCurves(cds, new_data = newdata)
  # dissimilarity_mat <- 1 - cor(t(smooth_res))
  # library(fpc)
  # pamk.best <- pamk(dissimilarity_mat, diss = T)
  # cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
  # plot(pam(dissimilarity_mat, pamk.best$nc))
  # 
  # library(mclust)
  # # Run the function to see how many clusters
  # # it finds to be optimal, set it to search for
  # # at least 1 model and up 20.
  # d_clust <- Mclust(as.matrix(d), G=1:20)
  # m.best <- dim(d_clust$z)[2]
  # cat("model-based optimal number of clusters:", m.best, "\n")
  # # 4 clusters
  # plot(d_clust)
  # 
  # fitting_res <- plot_pseudotime_heatmap(cds, return_heatmap = T)
  # annotation_row <- data.frame(Cluster=factor(cutree(fitting_res$tree_row, 6)))
  # 
  # smooth_avg_clust_exp <- matrix(0, nrow = nrow(newdata), ncol = 6)
  # for(i in 1:6)
  # {
  #   smooth_avg_clust_exp[, i] <- rowMeans(scale(t(exprs(cds)[row.names(subset(annotation_row, Cluster == i)), ])), na.rm = T) # smooth_res
  # }
  # smooth_avg_clust_exp <- as.data.frame(smooth_avg_clust_exp)
  # dimnames(smooth_avg_clust_exp) <- list(colnames(cds), c(paste0("cluster_", 1:6)))
  # run RDI: 
  exprs_data <- t(exprs(cds)[, ]) # get the expression matrix 
  
  # exprs_data <- t(exprs(cds)[TF_vec_names, ])
  # TF_cluster_exprs_data <- cbind(smooth_avg_clust_exp, TF_exprs_data)
  
  # TF_target_pair <- expand.grid(TF_vec_names, paste0("cluster_", 1:6), stringsAsFactors = F)
  tmp <- rbind(TF_pair, TF_target_pair)
  
  tmp[, 1] <- match(tmp[, 1], colnames(exprs_data))
  tmp[, 2] <- match(tmp[, 2], colnames(exprs_data))
  # 
  all_pairwise_gene <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index
  noise = matrix(rnorm(mean = 0, sd = 1e-10, nrow(exprs_data) * ncol(exprs_data)), nrow = nrow(exprs_data)) # add random noise to break the tie 

  a <- Sys.time()
  RDI_parallel_res <- calculate_rdi(exprs_data, delays = 10, super_graph = NULL, method = 1) # run RDI 
  b <- Sys.time()
  
  return(RDI_parallel_res)
  # colnames(RDI_parallel_res$RDI) <- colnames(exprs_data)
  # match(c("Irf8", "Gfi1"), names(sort(apply(RDI_parallel_res$RDI[TF_vec_id, paste0("cluster_", 1:6)], 1, function(x) sum(x[x > 0], na.rm = T) ))))
  # match(TF_vec_id, names(sort(apply(RDI_parallel_res$RDI[, ], 1, function(x) sum(x[x > 0], na.rm = T) ))))
  # 
  # names(sort(apply(RDI_parallel_res$RDI[, ], 1, function(x) sum(x[x > 0], na.rm = T) )))[sort(match(TF_vec_id, 
  #       names(sort(apply(RDI_parallel_res$RDI[, ], 1, function(x) sum(x[x > 0], na.rm = T) )))))]
  # # 
  # qplot(sort(match(TF_vec_id, names(sort(apply(RDI_parallel_res$RDI[, ], 1, function(x) sum(x[x > 0], na.rm = T) ))))))
  # 
  # # a <- Sys.time()
  # # RDI_parallel_res <- calculate_rdi(TF_cluster_exprs_data, delays = 10, super_graph = all_pairwise_gene, method = 1) # + noise c(5, 10, 15)
  # # b <- Sys.time()
  # # 
  # # RDI_parallel_res <- calculate_rdi(exprs_data[1:5, 1:5] + noise[1:5, 1:5], delays = 5, method = 1)
  # # 
  # # a <- Sys.time()
  # # RDI_parallel_res <- calculate_rdi(exprs_data + noise, delays = 10, method = 1)
  # # b <- Sys.time()
  # 
  # # 5. gene ranking by sum up all outgoing weights
  # sum_outgoing_edges <- apply(RDI, 1, sum) #
  # 
  # # 6. post-processing the network (use CLR or other methods to make the network sparse)
  # clr(RDI) # clr(cRDI)
  # 
  # # other method
  # 
  # # 6. visualize the network:
  # g <- igraph::graph_from_adjacency_matrix(RDI, mode = 'directed', weighted = T)
  # plot(g, layout = layout.fruchterman.reingold(g), vertex.size=2, vertex.label.size=4)
}
# 
# 
# library(InformationEstimator)
# library(monocle)
# library(destiny)
# 
# a <- Sys.time()
# res <- wired(Olsson_monocyte_cds, TF_vec_id, informative_genes)
# b <- Sys.time()
# # # cluster of genes 
# # match(c("Irf8", "Gfi1"), names(sort(apply(RDI_parallel_res$RDI[TF_vec_id, paste0("cluster_", 1:6)], 1, function(x) sum(x[x > 0]) ))))
# # 
# # mlt_smooth_avg_clust_exp <- melt(smooth_avg_clust_exp)
# # mlt_smooth_avg_clust_exp$Time <- rep(1:199, times = 6)
# # 
# # qplot(Time, value, data = mlt_smooth_avg_clust_exp, color = variable)
# # 
# # # a table between TFs to the clusters (make a heatmap)
# # sum_rdi_res <- matrix(0, nrow = length(TF_vec_id), ncol = 6)
# # cluster_avg_sum_rdi_res <- matrix(0, nrow = length(TF_vec_id), ncol = 6)
# # 
# # for(i in 1:6)
# # {
# #   for(j in 1:length(TF_vec_id))
# #   {
# #     #sum_rdi_res[j, i] <- mean(mono_RDI_10[TF_vec_id[j], row.names(subset(annotation_row, Cluster == i))], na.rm = T)
# #     cluster_avg_sum_rdi_res[j, i] <- mean(RDI_parallel_res$RDI[TF_vec_id[j], paste0('cluster_', i)], na.rm = T)
# #   }
# # }
# # 
# # dimnames(sum_rdi_res) <- list(TF_vec_id, paste0('cluster_', 1:6))
# # pheatmap::pheatmap(sum_rdi_res, cluster_rows = F, cluster_cols = F)
# # 
# # dimnames(cluster_avg_sum_rdi_res) <- list(TF_vec_id, paste0('cluster_', 1:6))
# # pheatmap::pheatmap(cluster_avg_sum_rdi_res, cluster_rows = F, cluster_cols = F)
# # 
# # ################################################################################################################################################
# # # run RDI from all genes to the six categories: 
# # ################################################################################################################################################
# # TF_exprs_data <- scale(t(exprs(cds)[, ]))
# # TF_cluster_exprs_data <- cbind(smooth_avg_clust_exp, TF_exprs_data)
# # 
# # TF_target_pair <- expand.grid(row.names(cds), paste0("cluster_", 1:6), stringsAsFactors = F)
# # tmp <- rbind(TF_pair, TF_target_pair)
# # 
# # tmp[, 1] <- match(tmp[, 1], colnames(TF_cluster_exprs_data))
# # tmp[, 2] <- match(tmp[, 2], colnames(TF_cluster_exprs_data))
# # 
# # a <- Sys.time()
# # RDI_parallel_res <- calculate_rdi(TF_cluster_exprs_data, delays = 10, super_graph = TF_target_pair, method = 1) # + noise c(5, 10, 15)
# # b <- Sys.time()
# # 
# # # run it 

# #################################################################################################################################################
# # test the errors generated during the real tests
# #################################################################################################################################################
# library(R.matlab)
# library(monocle)
# library(InformationEstimator)
# library(destiny)
# library(diffusionMap)
# library(xacHelper)
# library(reshape2)
# library(stringr)

# source("/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Scripts/function.R", echo = T)
# oligodendrocytedifferent_cell_pseudotime <- calStatistics(oligodendrocytedifferent_cell_pseudotime, reference_network_pvals) # this one crashes at R11 

# subset_exprs_mat <- as.matrix(subset(oligodendrocytedifferent_cell_pseudotime, V2 == "R11")[, -c(1:2)]) # remove first two columns: gene_name / run 
# noise = matrix(rnorm(mean = 0, sd = 1e-10, nrow(subset_exprs_mat) * ncol(subset_exprs_mat)), nrow = nrow(subset_exprs_mat)) # add noise to break tie 

# subset_exprs_mat_noise <- subset_exprs_mat + noise 
# subset_exprs_mat_noise <- subset_exprs_mat_noise - min(subset_exprs_mat_noise)

# a <- Sys.time() # calculate RDI 
# RDI_parallel_res <- calculate_rdi(t(subset_exprs_mat_noise), delays = c(1, 2, 3))
# b <- Sys.time()

# a <- Sys.time() # calculate cRDI 
# cRDI_parallel_res <- calculate_conditioned_rdi(t(subset_exprs_mat_noise), rdi_list = RDI_parallel_res)
# b <- Sys.time()

# RDI_parallel_res_list <- list()
# cRDI_parallel_res_list <- list() 
# subset(oligodendrocytedifferent_cell_pseudotime, V2 == 'R11')

# oligodendrocytedifferent_cell_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/oligodendrocytedifferent_cell_pseudotime.txt', sep = '\t')

# source("/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Scripts/analysis_sparsity_network.R", echo = T)
# set.seed(2017)
# oligodendrocytedifferent_cell_pseudotime <- calStatistics(oligodendrocytedifferent_cell_pseudotime, reference_network_pvals)

# astrocytedifferent_cell_realtime_zero_forcing <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/astrocytedifferent_cell_realtime_zero_forcing.txt', sep = '\t')
# set.seed(2017)
# astrocytedifferent_cell_realtime_zero_forcing <- calStatistics(astrocytedifferent_cell_realtime_zero_forcing, reference_network_pvals) # crash at run 22

# ######################################################################################################################################################################################################
# # a. same cell but with Pseudotime
# neurondifferent_cell_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/neurondifferent_cell_pseudotime.txt', sep = '\t')
# astrocytedifferent_cell_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/astrocytedifferent_cell_pseudotime.txt', sep = '\t')
# oligodendrocytedifferent_cell_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/oligodendrocytedifferent_cell_pseudotime.txt', sep = '\t')

# # b. different cell real time but with zero-forcing
# neurondifferent_cell_realtime_zero_forcing <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/neurondifferent_cell_realtime_zero_forcing.txt', sep = '\t')
# astrocytedifferent_cell_realtime_zero_forcing <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/astrocytedifferent_cell_realtime_zero_forcing.txt', sep = '\t')
# oligodendrocytedifferent_cell_realtime_zero_forcing <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/oligodendrocytedifferent_cell_realtime_zero_forcing.txt', sep = '\t')

# # c. different cell but with time sampled from a Guassian distribution
# neurondifferent_cell_gaussian_time <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/neurondifferent_cell_gaussian_time.txt', sep = '\t')
# astrocytedifferent_cell_gaussian_time <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/astrocytedifferent_cell_gaussian_time.txt', sep = '\t')
# oligodendrocytedifferent_cell_gaussian_time <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/oligodendrocytedifferent_cell_gaussian_time.txt', sep = '\t')

# # d. different cell but with time sampled from a Guassian distribution
# neurondifferent_cell_gaussian_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/neurondifferent_cell_gaussian_pseudotime.txt', sep = '\t')
# astrocytedifferent_cell_gaussian_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/astrocytedifferent_cell_gaussian_pseudotime.txt', sep = '\t')
# oligodendrocytedifferent_cell_gaussian_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/oligodendrocytedifferent_cell_gaussian_pseudotime.txt', sep = '\t')

# # e. different cell but with cells with time sampled from a Guassian distribution; Pseudotime is used
# neurondifferent_cell_gaussian_zero_time_forcing <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/neurondifferent_cell_gaussian_time_zero_forcing.txt', sep = '\t')
# astrocytedifferent_cell_gaussian_zero_time_forcing <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/astrocytedifferent_cell_gaussian_time_zero_forcing.txt', sep = '\t')
# oligodendrocytedifferent_cell_gaussian_zero_time_forcing <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/oligodendrocytedifferent_cell_gaussian_time_zero_forcing.txt', sep = '\t')

# # e. different cell but with cells with time sampled from a Guassian distribution; Pseudotime is used
# neurondifferent_cell_gaussian_zero_forcing_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/neurondifferent_cell_gaussian_zero_forcing_pseudotime.txt', sep = '\t')
# astrocytedifferent_cell_gaussian_zero_forcing_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/astrocytedifferent_cell_gaussian_zero_forcing_pseudotime.txt', sep = '\t')
# oligodendrocytedifferent_cell_gaussian_zero_forcing_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/oligodendrocytedifferent_cell_gaussian_zero_forcing_pseudotime.txt', sep = '\t')
# ######################################################################################################################################################################################################
# set.seed(2017)
# neurondifferent_cell_pseudotime <- calStatistics(neurondifferent_cell_pseudotime, reference_network_pvals)
# astrocytedifferent_cell_pseudotime <- calStatistics(astrocytedifferent_cell_pseudotime, reference_network_pvals)
# oligodendrocytedifferent_cell_pseudotime <- calStatistics(oligodendrocytedifferent_cell_pseudotime, reference_network_pvals) # this one crashes at R11 

# neurondifferent_cell_realtime_zero_forcing <- calStatistics(neurondifferent_cell_realtime_zero_forcing, reference_network_pvals)
# astrocytedifferent_cell_realtime_zero_forcing <- calStatistics(astrocytedifferent_cell_realtime_zero_forcing, reference_network_pvals) # crash at run 22
# oligodendrocytedifferent_cell_realtime_zero_forcing <- calStatistics(oligodendrocytedifferent_cell_realtime_zero_forcing, reference_network_pvals)

# neurondifferent_cell_gaussian_time <- calStatistics(neurondifferent_cell_gaussian_time, reference_network_pvals)
# astrocytedifferent_cell_gaussian_time <- calStatistics(astrocytedifferent_cell_gaussian_time, reference_network_pvals)
# oligodendrocytedifferent_cell_gaussian_time <- calStatistics(oligodendrocytedifferent_cell_gaussian_time, reference_network_pvals)

# neurondifferent_cell_gaussian_pseudotime <- calStatistics(neurondifferent_cell_gaussian_pseudotime, reference_network_pvals)
# astrocytedifferent_cell_gaussian_pseudotime <- calStatistics(astrocytedifferent_cell_gaussian_pseudotime, reference_network_pvals)
# oligodendrocytedifferent_cell_gaussian_pseudotime <- calStatistics(oligodendrocytedifferent_cell_gaussian_pseudotime, reference_network_pvals)

# neurondifferent_cell_gaussian_zero_time_forcing <- calStatistics(neurondifferent_cell_gaussian_zero_time_forcing, reference_network_pvals)
# astrocytedifferent_cell_gaussian_zero_time_forcing <- calStatistics(astrocytedifferent_cell_gaussian_zero_time_forcing, reference_network_pvals)
# oligodendrocytedifferent_cell_gaussian_zero_time_forcing <- calStatistics(oligodendrocytedifferent_cell_gaussian_zero_time_forcing, reference_network_pvals)

# neurondifferent_cell_gaussian_zero_forcing_pseudotime <- calStatistics(neurondifferent_cell_gaussian_zero_forcing_pseudotime, reference_network_pvals)
# astrocytedifferent_cell_gaussian_zero_forcing_pseudotime <- calStatistics(astrocytedifferent_cell_gaussian_zero_forcing_pseudotime, reference_network_pvals)
# oligodendrocytedifferent_cell_gaussian_zero_forcing_pseudotime <- calStatistics(oligodendrocytedifferent_cell_gaussian_zero_forcing_pseudotime, reference_network_pvals)

# ######################################################################################################################################################################################################
# # a. same cell but with Pseudotime
# neuronsame_cell_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/neuronsame_cell_pseudotime.txt', sep = '\t')
# astrocytesame_cell_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/astrocytesame_cell_pseudotime.txt', sep = '\t')
# oligodendrocytesame_cell_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/oligodendrocytesame_cell_pseudotime.txt', sep = '\t')

# # b. same cell real time but with zero-forcing
# neuronsame_cell_realtime_zero_forcing <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/neuronsame_cell_realtime_zero_forcing.txt', sep = '\t')
# astrocytesame_cell_realtime_zero_forcing <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/astrocytesame_cell_realtime_zero_forcing.txt', sep = '\t')
# oligodendrocytesame_cell_realtime_zero_forcing <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/oligodendrocytesame_cell_realtime_zero_forcing.txt', sep = '\t')

# # c. same cell but with time sampled from a Guassian distribution
# neuronsame_cell_gaussian_time <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/neuronsame_cell_gaussian_time.txt', sep = '\t')
# astrocytesame_cell_gaussian_time <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/astrocytesame_cell_gaussian_time.txt', sep = '\t')
# oligodendrocytesame_cell_gaussian_time <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/oligodendrocytesame_cell_gaussian_time.txt', sep = '\t')

# # d. same cell but with time sampled from a Guassian distribution
# neuronsame_cell_gaussian_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/neuronsame_cell_gaussian_pseudotime.txt', sep = '\t')
# astrocytesame_cell_gaussian_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/astrocytesame_cell_gaussian_pseudotime.txt', sep = '\t')
# oligodendrocytesame_cell_gaussian_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/oligodendrocytesame_cell_gaussian_pseudotime.txt', sep = '\t')

# # e. same cell but with cells with time sampled from a Guassian distribution; Pseudotime is used
# neuronsame_cell_gaussian_zero_time_forcing <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/neuronsame_cell_gaussian_zero_time_forcing.txt', sep = '\t')
# astrocytesame_cell_gaussian_zero_time_forcing <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/astrocytesame_cell_gaussian_zero_time_forcing.txt', sep = '\t')
# oligodendrocytesame_cell_gaussian_zero_time_forcing <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/oligodendrocytesame_cell_gaussian_zero_time_forcing.txt', sep = '\t')

# # e. same cell but with cells with time sampled from a Guassian distribution; Pseudotime is used
# neuronsame_cell_gaussian_zero_forcing_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/neuronsame_cell_gaussian_zero_forcing_pseudotime.txt', sep = '\t')
# astrocytesame_cell_gaussian_zero_forcing_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/astrocytesame_cell_gaussian_zero_forcing_pseudotime.txt', sep = '\t')
# oligodendrocytesame_cell_gaussian_zero_forcing_pseudotime <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/oligodendrocytesame_cell_gaussian_zero_forcing_pseudotime.txt', sep = '\t')
# ######################################################################################################################################################################################################
# set.seed(2017)
# neuronsame_cell_pseudotime <- calStatistics(neuronsame_cell_pseudotime, reference_network_pvals)
# astrocytesame_cell_pseudotime <- calStatistics(astrocytesame_cell_pseudotime, reference_network_pvals)
# oligodendrocytesame_cell_pseudotime <- calStatistics(oligodendrocytesame_cell_pseudotime, reference_network_pvals)

# neuronsame_cell_realtime_zero_forcing <- calStatistics(neuronsame_cell_realtime_zero_forcing, reference_network_pvals)
# astrocytesame_cell_realtime_zero_forcing <- calStatistics(astrocytesame_cell_realtime_zero_forcing, reference_network_pvals)
# oligodendrocytesame_cell_realtime_zero_forcing <- calStatistics(oligodendrocytesame_cell_realtime_zero_forcing, reference_network_pvals)

# neuronsame_cell_gaussian_time <- calStatistics(neuronsame_cell_gaussian_time, reference_network_pvals)
# astrocytesame_cell_gaussian_time <- calStatistics(astrocytesame_cell_gaussian_time, reference_network_pvals)
# oligodendrocytesame_cell_gaussian_time <- calStatistics(oligodendrocytesame_cell_gaussian_time, reference_network_pvals)

# neuronsame_cell_gaussian_pseudotime <- calStatistics(neuronsame_cell_gaussian_pseudotime, reference_network_pvals)
# astrocytesame_cell_gaussian_pseudotime <- calStatistics(astrocytesame_cell_gaussian_pseudotime, reference_network_pvals)
# oligodendrocytesame_cell_gaussian_pseudotime <- calStatistics(oligodendrocytesame_cell_gaussian_pseudotime, reference_network_pvals)

# neuronsame_cell_gaussian_zero_time_forcing <- calStatistics(neuronsame_cell_gaussian_zero_time_forcing, reference_network_pvals)
# astrocytesame_cell_gaussian_zero_time_forcing <- calStatistics(astrocytesame_cell_gaussian_zero_time_forcing, reference_network_pvals)
# oligodendrocytesame_cell_gaussian_zero_time_forcing <- calStatistics(oligodendrocytesame_cell_gaussian_zero_time_forcing, reference_network_pvals)

# neuronsame_cell_gaussian_zero_forcing_pseudotime <- calStatistics(neuronsame_cell_gaussian_zero_forcing_pseudotime, reference_network_pvals)
# astrocytesame_cell_gaussian_zero_forcing_pseudotime <- calStatistics(astrocytesame_cell_gaussian_zero_forcing_pseudotime, reference_network_pvals)
# oligodendrocytesame_cell_gaussian_zero_forcing_pseudotime <- calStatistics(oligodendrocytesame_cell_gaussian_zero_forcing_pseudotime, reference_network_pvals)

# ######################################################################################################################################################################################################







