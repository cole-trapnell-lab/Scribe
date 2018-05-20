#test InformationEstimator

message('test vd')
library(Scribe)
vd(10)

# 0.9361576864649548

message('test get_NN_2Set')
library(Scribe)
load('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/XYZ.txt')
XYZ <- t(XYZ)
dimensioon <- 3; ND <- 100000; NQ <- 100000; k <- 5; eps <- 0; searchtypeInt <- 1; tree_type <- 1; radius <- 0; 
nn.idx <- 1L; dists <- 0.1;
res <- get_NN_2Set(XYZ[, ], XYZ[, ], 3L, 100000L, 100000L, 5L, as.double(0.0), 1L, 0L, as.double(0.0), dists, T)

data_nmslib <- read.table("/Users/xqiu/Dropbox (Personal)/Projects/Genes_Inference_in_Cell_Differentiation_Process/notebook/data_nmslib.txt", sep = '\t')
res <- get_NN_2Set(XYZ[, ], XYZ[, ], 3L, 100000L, 100000L, 5L, as.double(0.0), 1L, 0L, as.double(0.0), dists, T)
a <- Sys.time(); res <- Neighbour(data_nmslib, data_nmslib, 1000, cores = 1); b <- Sys.time()
# 

message('test entropy')
library(Scribe)
entropy(as.matrix(XYZ[, 1]), 5)

# 1.1298735949169973

message('test mi')
library(Scribe)
a <- Sys.time()
mi(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 2]), 5, normalize = 0)
mi(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 3]), 5, normalize = 0)
mi(as.matrix(XYZ[, 2]), as.matrix(XYZ[, 3]), 5, normalize = 0)
mi(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 1]), 5, normalize = 0)
mi(as.matrix(XYZ[, 3]), as.matrix(XYZ[, 3]), 5, normalize = 0)
mi(as.matrix(XYZ[, 2]), as.matrix(XYZ[, 2]), 5, normalize = 0)
b <- Sys.time()
b - a 

# 1.88421442831
# 1.4081763245
# 1.82672752049
# 10.0068027965
# 10.0068027965
# 10.0068027965

# Time difference of 16.8602 secs

message('test cmi')
library(Scribe)
a <- Sys.time()
cmi(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 2]), as.matrix(XYZ[, 3]), 5, normalize = 0)$cmi_res;
cmi(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 3]), as.matrix(XYZ[, 2]), 5, normalize = 0)$cmi_res;
cmi(as.matrix(XYZ[, 2]), as.matrix(XYZ[, 3]), as.matrix(XYZ[, 1]), 5, normalize = 0)$cmi_res;
b <- Sys.time()
b - a 

# 0.47189332607
# 0.00569163896183
# 0.411085051308

# Time difference of 35.11695 secs

run_vec <- rep(1, ncol(XYZ))
tmp <- expand.grid(1:ncol(XYZ), 1:ncol(XYZ), stringsAsFactors = F)
super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index
super_graph <- super_graph[, c(2, 1)]

rdi_list <- calculate_rdi_cpp_wrap(XYZ, delay = c(1), as.matrix(super_graph), method = 1, turning_points = 0)
crdi_list <- calculate_conditioned_rdi(XYZ, as.matrix(super_graph), rdi_list, 1)


run_vec <- rep(1, 401)
run_vec <- rep(1, ncol(neuron_sim_cds))
tmp <- expand.grid(1:ncol(neuron_sim_cds), 1:ncol(neuron_sim_cds), stringsAsFactors = F)
super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index
super_graph <- super_graph[, c(2, 1)]

res <- read.csv("/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Python_code/simulation_expr_mat_nonlinear_n_4_step_400_N_end_400", sep = '\t', row.names = 1)
data <- res
neuron_rdi_list2 <- calculate_rdi_multiple_run_cpp(t(data), delay = c(1), run_vec - 1, as.matrix(super_graph), method = 1, turning_points = 0) #* 100 + noise

neuron_rdi_list <- calculate_rdi_cpp_wrap(t(data), delay = c(1), as.matrix(super_graph), method = 1, turning_points = 0) #* 100 + noise

noise <- matrix(rnorm(ncol(data) * 2, sd = 0), ncol = 2)
mi(matrix(as.numeric(data['Pax6', ]), ncol = 1), matrix(as.numeric(data['Mash1', ]), ncol = 1), k = 5, normalize = 0)
mi(matrix(as.numeric(data['Brn2', ]) + noise[, 1], ncol = 1), matrix(as.numeric(data['Aldh1L', ]) + noise[, 1], ncol = 1), k = 5, normalize = 0)
mi(matrix(as.numeric(data['Brn2', ]), ncol = 1), matrix(as.numeric(data['Hes5', ]), ncol = 1), k = 5, normalize = 0)

res <- read.csv("/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Python_code/test_data.txt", sep = '\t')

test <- t(data[c('Brn2', 'Aldh1L'), ]) 
test <- unique(test)

mi(matrix(test[, 1], ncol = 1), matrix(test[, 2], ncol = 1), k = 5, normalize = 0)

test_res <- RANN::nn2(test, k = 6)
dists <- 0.1

# get_NN_2Set_cpp(const NumericMatrix& data, const NumericMatrix& query_data, int& D, int& ND, int& NQ, int& K, double& EPS,
#                 int& SEARCHTYPE, int& USEBDTREE, double& SQRAD, NumericVector& distances) // const NumericVector& query,  IntegerVector& nn_index, 

res <- get_NN_2Set(test, test, 2L, 402, 402, 5L, as.double(0.0), 1L, 0L, as.double(0.0), dists, T)
get_points_in_radius_cpp(test, test, 2L, 402, 402, 5L, as.double(0.0), 1L, 0L, as.double(0.0), res$distances, T)
test_res <- RANN::nn2(test, k = 6)

################################################################################################################################################################################
# test the relationship between cmi value and the pseudotime 
################################################################################################################################################################################

################################################################################################################################################################################
# real data 
################################################################################################################################################################################
lung <- load_lung()
lung_AT1 <- lung[, pData(lung)$State %in% c(2, 3)]
lung_AT1 <- lung_AT1[, order(pData(lung_AT1)$Pseudotime)]
delay <- 5
x <- exprs(lung_AT1[2, 5:(ncol(lung_AT1) - 1)])
y <- exprs(lung_AT1[45, 6:(ncol(lung_AT1))])
z <- exprs(lung_AT1[45, 5:(ncol(lung_AT1) - 1)])
cmi_res <- cmi(t(as.matrix(x)), t(as.matrix(y)), t(as.matrix(z)), k = 5)
qplot(1:length(cmi_res$information_samples), cmi_res$information_samples)

gene_num <- nrow(lung)
pseudotime_cmi_res_df <- matrix(nrow = gene_num * (gene_num - 1), ncol = length(cmi_res$information_samples))
cnt <- 1

################################################################################################################################################################################
# make a heatmap for all pairwise relationship
################################################################################################################################################################################
for(i in 1:nrow(lung)) {
  message('current i is ', i)
  for(j in 1:nrow(lung)) {
    if(i == j)
      next; 
    x <- exprs(lung_AT1[i, delay:(ncol(lung_AT1) - 1)])
    y <- exprs(lung_AT1[j, (delay + 1):(ncol(lung_AT1))])
    z <- exprs(lung_AT1[j, delay:(ncol(lung_AT1) - 1)])
    cmi_res <- cmi(t(as.matrix(x)), t(as.matrix(y)), t(as.matrix(z)), k = 5)
    pseudotime_cmi_res_df[cnt, ] <- cmi_res$information_samples
    cnt <- cnt + 1
  }
}

# make smooth curves for the temporal variation of the information samples 
# only check for the positive gene-pairs
mean_cmi <- rowMeans(pseudotime_cmi_res_df)
test <- di::smooth_genes(t(pseudotime_cmi_res_df[which(mean_cmi > 0), ]), window_size = 20)
pheatmap::pheatmap(t(test), cluster_cols = F)

################################################################################################################################################################################
# simulation data 
################################################################################################################################################################################
load('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/RData/neuron_network')
gene_name_vec <- c('Pax6', 'Mash1', 'Brn2', 'Zic1', 'Tuj1', 'Hes5', 'Scl', 'Olig2', 'Stat3', 'Myt1L', 'Aldh1L', 'Sox8', 'Mature')

cell_simulate <- readMat('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/mat_data/cell_simulate.mat')
all_cell_simulation <- cell_simulate$cell.simulate[, 1:400, ] #time 0-20 are the period two branches appear
example_data <- all_cell_simulation[, , 1]

tmp <- expand.grid(1:ncol(t(example_data)), 1:ncol(t(example_data)), stringsAsFactors = F)
super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index
super_graph <- super_graph[, c(2, 1)]

delay <- 1
x <- example_data[1, delay:(ncol(example_data) - 1)]
y <- example_data[2, (delay + 1):(ncol(example_data))]
z <- example_data[2, (delay):(ncol(example_data) - 1)]
cmi_res <- cmi(matrix(x, ncol = 1), matrix(y, ncol = 1), matrix(z, ncol = 1), k = 5)
qplot(1:length(cmi_res$information_samples), cmi_res$information_samples)
cmi_res$information_samples

# simulation data (when n = 1)
data_ori <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Python_code/simulation_expr_mat_n_equal_1.txt', sep = '\t', skip = 1)
x3d <- array(dim = c(13, 101, 20))

data_ori_gene_unique <- unique(data_ori$V1)

data_ori$V2
for(i in 1:dim(x3d)[1]) {
  for(j in 1:dim(x3d)[3]) {
    x3d[i, , j] <- as.numeric(data_ori[data_ori$V1 == data_ori_gene_unique[i] & data_ori$V2 == paste0("R", j), -c(1:2)]) # make sure data_ori$V2 has R#
  }
}

dim(x3d)  
example_data <- x3d[, , 1]
################################################################################################################################################################################
# show the gene expression 
################################################################################################################################################################################

# make a heatmap for all pairwise relationship
# real data 
gene_num <- nrow(example_data)
pseudotime_cmi_res_df <- matrix(nrow = gene_num * (gene_num - 1), ncol = length(cmi_res$information_samples))
cnt <- 1

cmi_res_df <- matrix(0, nrow = 13, ncol = 13)
for(i in 1:nrow(example_data)) {
  message("current i is ", i)
  for(j in 1:nrow(example_data)) {
    if(i == j)
      next; 
    x <- example_data[i, delay:(ncol(example_data) - 1)]
    y <- example_data[j, (delay + 1):(ncol(example_data))]
    z <- example_data[j, (delay):(ncol(example_data) - 1)]
    cmi_res <- cmi(matrix(x, ncol = 1), matrix(y, ncol = 1), matrix(z, ncol = 1), k = 5)
    pseudotime_cmi_res_df[cnt, ] <- cmi_res$information_samples
    cnt <- cnt + 1
    
    cmi_res_df[i, j] <- cmi_res$cmi_res
  }
}

# only check for the positive gene-pairs
mean_cmi <- rowMeans(pseudotime_cmi_res_df)
test <- di::smooth_genes(t(pseudotime_cmi_res_df[which(mean_cmi > 0), ]), window_size = 20)
# set the name for the gene pair below: 
gene_uniq <- unique(c(as.character(neuron_network[, 1]), as.character(neuron_network[, 2])))
all_cmbns <- expand.grid(gene_name_vec, gene_name_vec)
valid_all_cmbns <- all_cmbns[all_cmbns$Var1 != all_cmbns$Var2, ]
valid_all_cmbns_df <- data.frame(pair = paste(tolower(valid_all_cmbns$Var2), tolower(valid_all_cmbns$Var1), sep = '_'), pval = 0)
row.names(valid_all_cmbns_df) <- valid_all_cmbns_df$pair
valid_all_cmbns_df[paste(tolower(neuron_network$V1), tolower(neuron_network$V2), sep = '_'), 2] <- 1

test <- t(test)
row.names(test) <- valid_all_cmbns_df$pair[which(mean_cmi > 0)]

pheatmap::pheatmap(test[, ], cluster_cols = F, cluster_rows = F, annotation_names_row = T)

row.names(pseudotime_cmi_res_df) <- valid_all_cmbns_df$pair
pheatmap::pheatmap(pseudotime_cmi_res_df, cluster_cols = F, cluster_rows = F, annotation_names_row = T)

################################################################################################################################################################################
# show the result for cRDI 
################################################################################################################################################################################

# calculate rdi values
run_vec <- rep(1, ncol(example_data))
a <- Sys.time()
rdi_list <- calculate_rdi_multiple_run_cpp(t(example_data), delay = 1, run_vec - 1, as.matrix(super_graph), method = 1) #calculate_rdi(data_noise, delay, method = 1)
b <- Sys.time()
rdi_time <- b - a

top_k_list <- extract_top_incoming_nodes_delays(rdi_list$max_rdi_value, rdi_list$max_rdi_delays, 1)
con_pseudotime_cmi_res_df <- matrix(nrow = gene_num * (gene_num - 1), ncol = length(cmi_res$information_samples))

# di::di_single_run_conditioned
# function (x, y, z, n = 10) 
# {
#   if (is.numeric(x)) 
#     x <- as.matrix(x)
#     if (is.numeric(y)) 
#       y <- as.matrix(y)
#       if (is.numeric(z)) 
#         z <- as.matrix(z)
#         if (ncol(x) != ncol(y)) 
#           stop("The number of time samples has to be the same for X and Y")
#           if (nrow(x) != nrow(z)) 
#             stop("The number of time samples has to be the same for X and all Zs")
#             tau <- n
#             tot_len <- nrow(x) - tau
#             x_past <- x[(tau):(tau - 1 + tot_len), ]
#           yz_past <- y[(tau):(tau - 1 + tot_len), ]
#           for (i in 1:n) {
#             if (i > 1) {
#               x_past <- cbind(x[(tau - i + 1):(tau - i + tot_len), 
#               ], x_past)
#               yz_past <- cbind(y[(tau - i + 1):(tau - i + tot_len), 
#               ], yz_past)
#             }
#             for (j in 1:ncol(z)) {
#               yz_past <- cbind(z[(tau - i + 1):(tau - i + tot_len), 
#                                j], yz_past)
#             }
#           }
#           return(cmi(x_past, y[(tau + 1):(tau + tot_len), ], yz_past))
# }


cnt <- 1
cmi_res_df <- matrix(0, nrow = 13, ncol = 13)
for(i in 1:nrow(example_data)) {
  message("current i is ", i)
  for(j in 1:nrow(example_data)) {
    if(i == j)
      next; 
    
    z_top_k_ind <- top_k_list$top_incoming_nodes[j, 1] + 1
    if(z_top_k_ind == i)
      z_top_k_ind <- top_k_list$top_incoming_nodes[j, 2] + 1
    x_past <- example_data[i, delay:(ncol(example_data) - 1)]
    yz_past <- example_data[j, (delay):(ncol(example_data) - 1)]
    yz_past <- rbind(example_data[z_top_k_ind, delay:(ncol(example_data) - 1)], yz_past)
    y_past <- example_data[j, (delay + 1):(ncol(example_data))]
    
    cmi_res <- cmi(matrix(x_past, ncol = 1), matrix(y_past, ncol = 1), t(yz_past), k = 5)
    con_pseudotime_cmi_res_df[cnt, ] <- cmi_res$information_samples
    cnt <- cnt + 1
    cmi_res_df[i, j] <- cmi_res$cmi_res
  }
}

con_rdi_res_test <- calculate_multiple_run_conditioned_rdi_wrap(t(example_data), as.matrix(super_graph), as.matrix(rdi_list$max_rdi_value), as.matrix(rdi_list$max_rdi_delays), run_vec - 1, 1)

cmi_res_df - con_rdi_res_test
a <- Sys.time()
rdi_list <- calculate_rdi_multiple_run_cpp(t(example_data), delay = c(1), run_vec - 1, as.matrix(super_graph), method = 1) #* 100 + noise
# rdi_list <- calculate_rdi_multiple_run_cpp(data_rem_dup, delay = c(1), run_vec[!duplicated(data)] - 1, as.matrix(super_graph_remove_dup), method = 1) # + noise_remove_dup
b <- Sys.time()
rdi_time <- b - a

dimnames(rdi_list$max_rdi_value) <- list(uniq_gene, uniq_gene)
con_rdi_res_test <- calculate_multiple_run_conditioned_rdi_wrap(t(example_data), as.matrix(super_graph), as.matrix(rdi_list$max_rdi_value), as.matrix(rdi_list$max_rdi_delays), run_vec - 1, 1)
dimnames(con_rdi_res_test) <- list(uniq_gene, uniq_gene)

row.names(con_pseudotime_cmi_res_df) <- valid_all_cmbns_df$pair
pheatmap::pheatmap(con_pseudotime_cmi_res_df, cluster_cols = F, cluster_rows = F, annotation_names_row = T)


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

library(Scribe)
library(destiny)
library(monocle)

load('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/XYZ.txt')
XYZ <- t(XYZ)
dimensioon <- 3; ND <- 100000; NQ <- 100000; k <- 5; eps <- 0; searchtypeInt <- 1; tree_type <- 1; radius <- 0; 

library(InformationEstimator)
N <- 10000
a <- Sys.time()
tmp <- expand.grid(1:ncol(XYZ), 1:ncol(XYZ), stringsAsFactors = F)
super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index
super_graph <- super_graph[, c(2, 1)]
calculate_rdi_cpp_wrap(as.matrix(XYZ[1:N, ]), delays = 1L, super_graph = as.matrix(super_graph), turning_points = 0, method = 1)
b <- Sys.time()
# 
# # entropy_cpp(XYZ[, 1], k = 5, N = 100000)
# # vd_cpp(10)
# nn.idx = 1L; dists = 1.0;
# 
# # res <- get_NN_2Set(as.double(XYZ[, 1]), as.double(XYZ[, 1]), 1L, 100000L, 100000L, 5L, as.double(0.0), 1L, 0L, as.double(0.0), nn.idx, dists, T)
# 
# 
# res <- get_NN_2Set(XYZ[, ], XYZ[, ], 3L, 100000L, 100000L, 5L, as.double(0.0), 1L, 0L, as.double(0.0), nn.idx, dists, T)

a <- Sys.time(); res <- get_NN_2Set(XYZ[, ], XYZ[, ], 3L, 100000L, 100000L, 1000L, as.double(0.0), 1L, 0L, as.double(0.0), nn.idx, dists); b <- Sys.time()
a <- Sys.time(); res <- Neighbour(XYZ, XYZ, 1000, cores = 8); b <- Sys.time()

nr <- 1000
nc <- 10
p <- LshParameterSetter$new(nr, nc)
X <- matrix(rnorm(nr * nc), nr, nc)
tab <- LshNnTable$new(t(X), p)
tab$find_nearest_neighbor(as.vector(X[1, ]))

