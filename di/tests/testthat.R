library(testthat)
library(di)
library(netbiov)
library(reshape2)
library(pscl)
library(boot)

test_check("di")

vd(100)

hbp_exprs = read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Test_dataset/hbp_mat.txt',
                        sep = '\t', header = T, row.names = 1)

x <- hbp_exprs[1:20, 1:20]
entropy(x, k = 2)

y <- hbp_exprs[21:40, 21:40]
mi(x, y)

z <- hbp_exprs[41:60, 41:60]
cmi(x, y, z)

di_single_run(x, y)
di_single_run(x[1:15, 1:15], y[1:15, 1:15], n = 3)

rdi_single_run(x, y) #test this
rdi_single_run(t(scale(t(x + 1))), t(scale(t(y + 1)))) #test this

rdi_single_run(x[, 1], y[, 1])
rdi_single_run(scale(x[, 1]), scale(y[, 1]))

rdi_single_run(scale(x, scale = F)[, 1], scale(y, scale = F)[, 1])
rdi_single_run(scale(x, center = F)[, 1], scale(y, center = F)[, 1])
rdi_single_run(scale_0_1_range(x)[, 1], scale_0_1_range(y)[, 1])
rdi_single_run(exp(x), exp(y))
rdi_single_run((x)/c(1:10), (y)/c(1:10))

z <- cbind(x[, 1], y[, 1])
z[, 1] <- c(-0.10586709, 0.02310539, 0.97479104, 0.84491918, 0.11360705, -1.22814749, -1.9416257 , -0.94071219, -0.31573972, -0.78625982, 1.44150499, 0.59099037, 0.9799445 , -0.88648595, 0.88125935, 0.16234335, -1.25395129, -1.88845447, -0.44738441, -0.37723717)
z[, 2] <- c(0.59876415, -0.21560821, -0.06456355, 0.56697927, 0.94749493, -0.17139318, -1.4249324 , -1.94804793, -0.20550459, -0.02151109, -1.30797509, 1.96236729, 0.57850622, 1.30707304, -1.1497937 , 1.0829491 , 0.30927887, -1.84287659, -1.75129567, 0.19247962)

rdi_single_run_conditioned(x, y, z, c(2, 1), d=1)

di::di_single_run_conditioned(x, y, z, n=1)

#debug the following:
rdi_many_runs(x, y)

order_lung_mat <- exprs(lung)[, order(pData(lung)$Pseudotime)] + 1
order_lung_mat <- order_lung_mat[apply(order_lung_mat, 1, function(x) sum(x > 1)> 150), ]
order_lung_mat <- order_lung_mat[!(row.names(order_lung_mat) %in%  c('ENSMUSG00000032489.11', 'ENSMUSG00000071415.6')), ] #order_lung_mat[, ]

lung_subset <- lung[row.names(order_lung_mat)[1:5], ]
pd <- pData(lung_subset[, colnames(order_lung_mat)]) #
fd <- fData(lung_subset[rownames(order_lung_mat)[1:5] , ]) #fData(lung[, ])

write.table(file = '/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Test_dataset//lung_order_mat.txt', order_lung_mat[1:5, ], col.names = T, sep = '\t', row.names = T, quote = F)
write.table(file = '/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Test_dataset/lung_gene_feature.txt', fd, col.names = T, sep = '\t', row.names = T, quote = F)
write.table(file = '/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Test_dataset/lung_cell_feature.txt',  pd, col.names = T, sep = '\t', row.names = T, quote = F)

##############################################################################################################################################################################
# smooth the data:
##############################################################################################################################################################################
test_smooth <- smooth_genes(t(exprs(lung_subset)[, order(pData(lung_subset)$Pseudotime)]), window_size = 20)
t(exprs(lung))[order(pData(lung)$Pseudotime), toupper('ensmusg00000054766.7')]
test_smooth[, toupper('ensmusg00000054766.7')]

##############################################################################################################################################################################
# test on the CMI, RDI on the dataset:
##############################################################################################################################################################################
# test on rdi:
# test on just single run:
calculate_rdi_corr(toupper('ensmusg00000015452.8'), toupper('ensmusg00000025283.9'),
                   test_smooth[1:165, toupper('ensmusg00000015452.8')] + 1, test_smooth[1:165, toupper('ensmusg00000025283.9')] + 1,
                   c(1, 21, 41), nrow(test_smooth) * (nrow(test_smooth) - 1))
#test on multiple run for the UN-conditioned RDI:

#test on conditioned rdi:
calculate_conditioned_rdi(toupper('ensmusg00000015452.8'), toupper('ensmusg00000025283.9'),
                   test_smooth[1:165, ] + 1, 20, c('ENSMUSG00000015452.8' = 1),
                   nrow(test_smooth[1:165, ]) * (nrow(test_smooth[1:165, ]) - 1))

#test on multiple run for the conditioned rdi:

#test on parallel verision for UN-conditioned RDI:
rdi_res <- calculate_and_write_pairwise_dmi(test_smooth[1:165, 1:5], delays = c(1, 21, 41), cores = 1, verbose = T)

#test on parallel verision for conditioned RDI:

calculate_conditioned_rdi(gene_pair[[1]], gene_pair[[2]], test_smooth, delays_max[index_name], top_k_incoming_id_list[[gene_pair[[2]]]], N_operations)
con_rdi_res <- calculate_and_write_pairwise_dmi_conditioned(test_smooth[1:165, unique(rdi_res$id_1)], rdi_res, cores = 1, k = 1, verbose = T)

##############################################################################################################################################################################
# benchmark with python implementation:
##############################################################################################################################################################################
python_rdi <- read.table("/Users/xqiu/Dropbox (Personal)/Projects/Genes_Inference_in_Cell_Differentiation_Process/notebook/output\\output_RDI.txt", header = T)
python_rdi_conditioned <- read.table("/Users/xqiu/Dropbox (Personal)/Projects/Genes_Inference_in_Cell_Differentiation_Process/notebook/output\\output_RDI_conditioned.txt", header = T)

python_rdi[order(python_rdi$CORRELATION), ]
python_rdi_conditioned[order(python_rdi_conditioned$RDI), ]

row.names(python_rdi) <- paste(toupper(python_rdi$Gene_1_ID), toupper(python_rdi$Gene_2_ID), sep = '_')

qplot(rdi_res[sort(row.names(python_rdi)), 3], python_rdi[sort(row.names(python_rdi)), 5]) + xlab('python') + ylab('R')
qplot(rdi_res[sort(row.names(python_rdi)), 4], python_rdi[sort(row.names(python_rdi)), 6])
qplot(rdi_res[sort(row.names(python_rdi)), 5], python_rdi[sort(row.names(python_rdi)), 7])
qplot(rdi_res[sort(row.names(python_rdi)), 6], python_rdi[sort(row.names(python_rdi)), 9])

row.names(python_rdi_conditioned) <- paste(toupper(python_rdi_conditioned$Gene_1_ID), toupper(python_rdi_conditioned$Gene_2_ID), sep = '_')

qplot(as.numeric(con_rdi_res[sort(row.names(python_rdi_conditioned)), 3]), python_rdi_conditioned[sort(row.names(python_rdi_conditioned)), 6]) + xlab('python') + ylab('R')
qplot(as.numeric(con_rdi_res[sort(row.names(python_rdi_conditioned)), 3]), python_rdi_conditioned[sort(row.names(python_rdi_conditioned)), 6]) + geom_text(label = sort(row.names(python_rdi_conditioned)))
qplot(as.numeric(con_rdi_res[sort(row.names(python_rdi_conditioned)), 4]), python_rdi_conditioned[sort(row.names(python_rdi_conditioned)), 5])

# different delay max pairs: (numerical issues)
different_delay_pair <- which(as.numeric(con_rdi_res[sort(row.names(python_rdi_conditioned)), 4]) != python_rdi_conditioned[sort(row.names(python_rdi_conditioned)), 5])
pair_id <- con_rdi_res[sort(row.names(python_rdi_conditioned)), ][different_delay_pair, ]
rdi_res[pair_id, ]
python_rdi[pair_id, ]

##############################################################################################################################################################################
# New things we can do: (test on the lung dataset)
##############################################################################################################################################################################
# 1. AUC
#test on parallel verision for UN-conditioned RDI:
test_smooth <- smooth_genes(t(exprs(lung)[, order(pData(lung)$Pseudotime)]), window_size = 20)
test_smooth <- test_smooth[1:165, 1:10]

#s
scale_0_1_range <- function(x) {
 (x - min(x)) / (max(x) - min(x))
}

test_smooth_0_1 <- apply(test_smooth, 2, function(x) scale_0_1_range(x))

test_smooth_center <- apply(test_smooth, 2, function(x) scale(x, scale = F))

rdi_res <- calculate_and_write_pairwise_dmi(test_smooth, delays = c(1, 21, 41), cores = 2, verbose = T)
rdi_res_0_1 <- calculate_and_write_pairwise_dmi(test_smooth_0_1, delays = c(1, 21, 41), cores = 2, verbose = T)

#test on parallel verision for conditioned RDI:
con_rdi_res <- calculate_and_write_pairwise_dmi_conditioned(test_smooth[, unique(rdi_res$id_1)], rdi_res, cores = 2, k = 1, verbose = T)
scale_con_rdi_res <- calculate_and_write_pairwise_dmi_conditioned(test_smooth[, unique(rdi_res$id_1)], rdi_res, cores = 2, k = 1, verbose = T)

# 2. Make a network plot
#a. heatmap
con_rdi_res$conditioned_rdi <- as.numeric(con_rdi_res$conditioned_rdi)
con_rdi_res$delay_max <- as.numeric(con_rdi_res$delay_max)

rdi_res_mat <- dcast(con_rdi_res[, 1:3], id_1 ~ id_2)

row.names(rdi_res_mat) <- rdi_res_mat$id_1
rdi_res_mat <- rdi_res_mat[, -1]

rdi_res_mat[!is.finite(as.matrix(rdi_res_mat))] <- 0
diag(rdi_res_mat) <- 0

#pdf('lung_rdi.pdf', width = 30, height = 30)
pheatmap(rdi_res_mat[, ], useRaster = T, cluster_cols = T, cluster_rows = T) #, annotation_col = F, annotation_row = F
#dev.off()

#b. network
# create a graph:
g1 <- graph_from_adjacency_matrix(as.matrix(rdi_res_mat), mode = "directed", weighted = T)
plot(g1, layout = layout_nicely(g1))

res <- level.plot(g1)

# 3. rank the genes by the cRDI values
#calculate the sum of incoming cRDI values and the sum of outgoing cRDI values:
#row: targets of the row.name gene; column: source of the colname gene
outgoing_rdi_sum <- apply(rdi_res_mat, 1, sum)
incoming_rdi_sum <- apply(rdi_res_mat, 2, sum)

sort(outgoing_rdi_sum)
sort(incoming_rdi_sum)

# 4. calculate the fixed size for the negative binomial distribution
#see the code:

# 5. use pscl package for fitting the moving everage data
hurdle_res <- hurdle(ENSMUSG00000004044.9 ~ 1, data = round(data.frame(test_smooth)))
zeroinfl_res <- zeroinfl(ENSMUSG00000004044.9 ~ pData(lung)$Pseudotime[1:166], data = round(data.frame(test_smooth)))

qplot(pData(lung)$Pseudotime[1:166], zeroinfl_res$fitted.values) #the resulting data is very smooth

#use the moving average with the hurdle/zeroinfl model for data fitting:
pscl_smooth_genes <- function(exprs, window_size = 40) {
  win_range <- nrow(exprs) - window_size
  exprs_smooth <- exprs[-c(1:window_size), ]

  res <- apply(exprs, 2, function(x) {
    tmp <- rep(NA, win_range)
    for(i in 0:(win_range)){
      df <- data.frame(expression = round(x[i + c(1:window_size)]))
      # print(x)
      tryCatch({
        tmp[i + 1] <- zeroinfl(expression ~ 1, data = df)$fitted.values
        tmp
      },
      #warning = function(w) { FM_fit },
      error = function(e) {
        # tmp[i + 1] <- glm.nb(expression ~ 1, data = df)$fitted.values
        # tmp
      })
      if(is.na(tmp[i + 1]))
        tmp[i + 1] <- glm.nb(expression ~ 1, data = df)$fitted.values
    }
    return(tmp)
  })

  return(res)
}

raw_data <- t(exprs(lung_subset)[, order(pData(lung_subset)$Pseudotime)])
pscl_res <- pscl_smooth_genes(raw_data, window_size = 20)

moving_avg_res <- smooth_genes(raw_data, window_size = 20)

ind <- 1
p0 <- qplot(1:nrow(raw_data), raw_data[, ind]) + ggtitle('raw data')
p1 <- qplot(1:nrow(pscl_res), pscl_res[, ind]) + ggtitle('pscl')
p2 <- qplot(1:nrow(moving_avg_res), moving_avg_res[, ind]) + ggtitle('moving avg')
p3 <- qplot(pscl_res[, ind], moving_avg_res[, ind]) + ggtitle('pscl vs. moving avg')

xacHelper::multiplot(plotlist = list(p0, p1, p2, p3))

ind <- 2
p0 <- qplot(1:nrow(raw_data), raw_data[, ind]) + ggtitle('raw data')
p1 <- qplot(1:nrow(pscl_res), pscl_res[, ind]) + ggtitle('pscl')
p2 <- qplot(1:nrow(moving_avg_res), moving_avg_res[, ind]) + ggtitle('moving avg')
p3 <- qplot(pscl_res[, ind], moving_avg_res[, ind]) + ggtitle('pscl vs. moving avg')

ind <- 3
p0 <- qplot(1:nrow(raw_data), raw_data[, ind]) + ggtitle('raw data')
p1 <- qplot(1:nrow(pscl_res), pscl_res[, ind]) + ggtitle('pscl')
p2 <- qplot(1:nrow(moving_avg_res), moving_avg_res[, ind]) + ggtitle('moving avg')
p3 <- qplot(pscl_res[, ind], moving_avg_res[, ind]) + ggtitle('pscl vs. moving avg')

xacHelper::multiplot(plotlist = list(p0, p1, p2, p3))

ind <- 4
p0 <- qplot(1:nrow(raw_data), raw_data[, ind]) + ggtitle('raw data')
p1 <- qplot(1:nrow(pscl_res), pscl_res[, ind]) + ggtitle('pscl')
p2 <- qplot(1:nrow(moving_avg_res), moving_avg_res[, ind]) + ggtitle('moving avg')
p3 <- qplot(pscl_res[, ind], moving_avg_res[, ind]) + ggtitle('pscl vs. moving avg')

xacHelper::multiplot(plotlist = list(p0, p1, p2, p3))

ind <- 5
p0 <- qplot(1:nrow(raw_data), raw_data[, ind]) + ggtitle('raw data')
p1 <- qplot(1:nrow(pscl_res), pscl_res[, ind]) + ggtitle('pscl')
p2 <- qplot(1:nrow(moving_avg_res), moving_avg_res[, ind]) + ggtitle('moving avg')
p3 <- qplot(pscl_res[, ind], moving_avg_res[, ind]) + ggtitle('pscl vs. moving avg')

xacHelper::multiplot(plotlist = list(p0, p1, p2, p3))

gene_expression_cell_num <- esApply(lung, 1, function(x) sum(x > 0))
lung_raw_data <- t(exprs(lung)[, order(pData(lung_subset)$Pseudotime)])
lung_pscl_res <- pscl_smooth_genes(lung_raw_data, window_size = 20)
lung_moving_avg_res <- smooth_genes(lung_raw_data, window_size = 20)

ind <- 'ENSMUSG00000078202.2'
ind <- 10

p0 <- qplot(1:nrow(lung_raw_data), lung_raw_data[, ind]) + ggtitle('raw data')
p1 <- qplot(1:nrow(lung_pscl_res), lung_pscl_res[, ind]) + ggtitle('pscl')
p2 <- qplot(1:nrow(lung_moving_avg_res), lung_moving_avg_res[, ind]) + ggtitle('moving avg')
p3 <- qplot(lung_pscl_res[, ind], lung_moving_avg_res[, ind]) + ggtitle('pscl vs. moving avg')

xacHelper::multiplot(plotlist = list(p0, p1, p2, p3))

# 6. implment the bootstrapping algorithm for obtaining the more robust results
library(boot)

# boot(data, statistic, R, sim = "ordinary", stype = c("i", "f", "w"),
#      strata = rep(1,n), L = NULL, m = 0, weights = NULL,
#      ran.gen = function(d, p) d, mle = NULL, simple = FALSE, ...,
#      parallel = c("no", "multicore", "snow"),
#      ncpus = getOption("boot.ncpus", 1L), cl = NULL)

#random sample the data but still preserve the ordering of the column

calculate_and_write_pairwise_dmi_wraper <- function(d, i, delays = c(1,2,5,10,15,20,25), supergraph = NULL, cores = 1, verbose = F) {
  print(i)
  data <- d[i, ]
  # print('before RDI')
  res <- calculate_and_write_pairwise_dmi(data, delays = delays, supergraph = supergraph, cores = cores, verbose = verbose)
  # print('after RDI')

  res_colnames <- colnames(res)
  rdi_res <- c(); for(index in 1:length(res_colnames)){
    names <- res_colnames[index]

    if(index > 2)
      rdi_res <- c(rdi_res, as.numeric(res[, names]))
    else
      rdi_res <- c(rdi_res, as.numeric(res[, names]))
  }
  return(rdi_res)
}

boot_res <- boot(test_smooth, calculate_and_write_pairwise_dmi_wraper, R = 199, stype = 'i')
boot_all_runs <- rbind(boot_res$t0, boot_res$t)

qplot(boot_all_runs[, 41]) # show the result

#bootstrapping result:
rdi_res <- calculate_and_write_pairwise_dmi(test_smooth, delays = c(1, 21, 41), cores = 2, verbose = T)


