################################################################################################################################################################################################
# test lmi and rdi, etc. 
################################################################################################################################################################################################
lung <- load_lung()
AT1_lung <- lung[, pData(lung)$State %in% c(2, 3)]

a <- Sys.time()
rdi_res <- calculate_rdi(AT1_lung, delays = 10, method = 1)
# rdi_res <- calculate_rdi(neuron_sim_cds, delays = 1, method = 1)
b <- Sys.time()
b - a 

a <- Sys.time()
crdi_res <- calculate_conditioned_rdi(AT1_lung, rdi_list = rdi_res)
# crdi_res <- calculate_conditioned_rdi(neuron_sim_cds, rdi_list = rdi_res)
b <- Sys.time()
b - a 

# all incoming node has the same RDI value for gene pair 0 and 157. Either of the gene may have no expression.

a <- Sys.time()
lmi_res <- calculate_rdi(AT1_lung, delays = 10, method = 2) #only one or two folder fast but we can use the very fast implementation from parinfogene package, etc. 
# lmi_res <- calculate_rdi(neuron_sim_cds, delays = 1, method = 2)
b <- Sys.time()
b - a 

a <- Sys.time()
mi_res <- parmigene::knnmi.all(exprs(AT1_lung)[, ], k = 5)
# mi_res <- parmigene::knnmi.all(exprs(neuron_sim_cds)[, ], k = 5)
b <- Sys.time()
b - a 

pdf('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Figures/supplementary_figures/RDI_LMI.pdf', height = 1.5, width = 1.5)
qplot(as.vector(lmi_res$RDI), as.vector(rdi_res$RDI), size = I(0.5), alpha = I(0.2)) + geom_abline() + xlab('LMI') + ylab('RDI') + 
  geom_smooth() + xacHelper:::monocle_theme_opts() + theme(axis.text.x = element_text(angle = 30))
dev.off()

pdf('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Figures/supplementary_figures/RDI_cRDI.pdf', height = 1.5, width = 1.5)
qplot(as.vector(rdi_res$RDI), as.vector(crdi_res), size = I(0.5), alpha = I(0.2)) + geom_abline() + xlab('RDI') + ylab('cRDI') + geom_smooth() + xacHelper:::monocle_theme_opts()
dev.off()

pdf('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Figures/supplementary_figures/MI_RDI.pdf', height = 1.5, width = 1.5)
qplot(as.vector(mi_res), as.vector(rdi_res$RDI), size = I(0.5), alpha = I(0.2)) + geom_abline() + xlab('MI') + ylab('RDI') + geom_smooth() + xacHelper:::monocle_theme_opts()
dev.off()

pdf('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Figures/supplementary_figures/MI_LMI.pdf', height = 1.5, width = 1.5)
qplot(as.vector(mi_res), as.vector(lmi_res$RDI), size = I(0.5), alpha = I(0.2)) + geom_abline() + xlab('MI') + ylab('LMI') + geom_smooth() + xacHelper:::monocle_theme_opts()
dev.off()

# write this in C ++ if this is important 
fast_lmi <- function(ordered_data, delay = 10) {
  mat_ori <- ordered_data[, -c((ncol(ordered_data) - 10):ncol(ordered_data))]
  mat_delay <- ordered_data[, 11:ncol(ordered_data)]
  # parmigene::knnmi.cross(mat_ori, mat_delay)
}

all_Olsson_granulocyte_cds_exprs <- t(exprs(Olsson_granulocyte_cds)[, order(granulocyte_dpt_res$pt)])

a <- Sys.time()
res <- parmigene::knnmi.all(exprs(Olsson_granulocyte_cds)[, order(granulocyte_dpt_res$pt)], k = 5)
b <- Sys.time()
b - a 
# if this work, we can just take the top 10 % edges and then do the cRDI or RDI 
# we can only use the extremely fast implementation of the MI code 

################################################################################################################################################################################################
# do the State Space Reconstruction on the neuron simulation 
################################################################################################################################################################################################

# lag 1, 2, 3 vs the original data 
load('/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision//RData/fig3.RData')

nao_sim_cds_subset <- nao_sim_cds[, c(1:100, 401:500, 801:900)]
qplot(exprs(nao_sim_cds_subset)['Mash1', 1:100], exprs(nao_sim_cds_subset)['Hes5', 1:100])
plot3d(exprs(nao_sim_cds_subset)['Mash1', ], exprs(nao_sim_cds_subset)['Hes5', ], exprs(nao_sim_cds_subset)['Pax6', ], col = c(rep(c('red', 'blue', 'green'), each = 100)))

plot3d(exprs(nao_sim_cds_subset)['Mash1', ], exprs(nao_sim_cds_subset)['Hes5', ], exprs(nao_sim_cds_subset)['Olig2', ], 
       col = c(rep(c('red', 'blue', 'green'), each = 100)), xlab = 'Mash1', ylab = 'Hes5', zlab = 'Olig2')
# add the state space reconstruction 

plot3d(exprs(nao_sim_cds_subset)['Mash1', c(1:98, 101:198, 201:298)], exprs(nao_sim_cds_subset)['Mash1', c(2:99, 102:199, 202:299)],
       exprs(nao_sim_cds_subset)['Mash1', c(3:100, 103:200, 203:300)], col = c(rep(c('red', 'blue', 'green'), each = 98)), 
       xlab = 'Mash1', ylab = 'Mash1 (delay 1)', zlab = 'Mash1 (delay 2)')
plot3d(exprs(nao_sim_cds_subset)['Hes5', c(1:98, 101:198, 201:298)], exprs(nao_sim_cds_subset)['Hes5', c(2:99, 102:199, 202:299)], 
       exprs(nao_sim_cds_subset)['Hes5', c(3:100, 103:200, 203:300)], col = c(rep(c('red', 'blue', 'green'), each = 98)),
       xlab = 'Hes5', ylab = 'Hes5 (delay 1)', zlab = 'Hes5 (delay 2)')
plot3d(exprs(nao_sim_cds_subset)['Scl', c(1:98, 101:198, 201:298)], exprs(nao_sim_cds_subset)['Scl', c(2:99, 102:199, 202:299)], 
       exprs(nao_sim_cds_subset)['Scl', c(3:100, 103:200, 203:300)], col = c(rep(c('red', 'blue', 'green'), each = 98)),
       xlab = 'Scl', ylab = 'Scl (delay 1)', zlab = 'Scl (delay 2)')

################################################################################################################################################################################################
# finalize the simulation analysis with drop-out, etc. 
################################################################################################################################################################################################

################################################################################################################################################################################################
# finalize the simulation analysis with drop-out, etc. 
################################################################################################################################################################################################


################################################################################################################################################################################################
# compare the C ++ code implmentation with the python implmentation 
################################################################################################################################################################################################
# transform the test dataset (gene name, run id, gene expression over time)

data_ori <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Python_code/test_data.txt')
data <- matrix(nrow = (ncol(data_ori) - 2) * 20, ncol = 3)
gene_vec <- data_ori$V1[-1]
uniq_gene <- levels(gene_vec)[-1]

for(i in 1:length(uniq_gene)) { #1:(ncol(data_ori) - 2)
  tmp <- t(as.matrix(data_ori[data_ori$V1 == uniq_gene[i], -(1:2)]))
  data[, i] <- tmp
}

run_vec <- rep(1:20, each = ncol(data_ori) - 2)

# run multiple rdi 
a <- Sys.time()
rdi_list <- calculate_rdi(data, delay = c(1, 2), method = 1)
b <- Sys.time()
rdi_time <- b - a 

tmp <- expand.grid(1:ncol(data), 1:ncol(data), stringsAsFactors = F)
super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index 

# calculate crdi values

# run multiple cRDI 
con_rdi_res_test <- calculate_multiple_run_conditioned_rdi_wrap(data, as.matrix(super_graph), as.matrix(rdi_list$max_rdi_value), as.matrix(rdi_list$max_rdi_delays), run_vec - 1, 1)

################################################################################################################################################################################################
# test the code: 2 run + 5 cells each (use this to verify our code)
################################################################################################################################################################################################

run_vec <- rep(1:2, each = 5)
tmp <- expand.grid(1:ncol(data), 1:ncol(data), stringsAsFactors = F)
super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index 
super_graph <- super_graph[, c(2, 1)]

# run multiple rdi 
a <- Sys.time()
rdi_list <- calculate_rdi_multiple_run_cpp(data[1:10, ], delay = c(1), run_vec - 1, as.matrix(super_graph), method = 1)
b <- Sys.time()
rdi_time <- b - a 

# calculate crdi values

# run multiple cRDI 
rdi_list$max_rdi_value <- matrix(c(-1, -0.016667, -0.012500, 0.075, -1, 0, -0.012500, 0.079167, -1), nrow = 3, byrow=F) # copy the result from python 
rdi_list$max_rdi_value <- matrix(c(-1, -0.012500, -0.016667, 0.075, -1, 0, -0.012500, 0.079167, -1), nrow = 3, byrow=F) # copy the result from python (meet the bug in python implementation)

con_rdi_res_test <- calculate_multiple_run_conditioned_rdi_wrap(data[1:10, ], as.matrix(super_graph), as.matrix(rdi_list$max_rdi_value), as.matrix(rdi_list$max_rdi_delays), run_vec - 1, 1)

# save the data for testing in python: 
data_for_python <- data.frame(rbind(t(data[1:5, ]), t(data[6:10, ])))
data_for_python <- cbind(data.frame(GENE_ID = rep(c("X", "Y", "Z"), 2), RUN_ID = rep(c(1), each = 3)), data_for_python); 
write.table(file = '/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Python_code/data_for_python.txt', data_for_python, sep = '\t', row.names = F, quote = F)













