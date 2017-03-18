#' RDI and the spearman correlation of expression
#'
#' This function estimates the RDI and the spearman correlation of expression (expr1, expr2) for two genes,
#' id1, di2 along the pseudo-time-series
#' @param id1 ensemble id for gene 1
#' @param id2 ensemble id for gene 2
#' @param expr1 number for nearest neighbors used in entropy calculation
#' @param expr2 number for nearest neighbors used in entropy calculation
#' @param delays number for nearest neighbors used in entropy calculation
#' @param N_operations number for nearest neighbors used in entropy calculation
#' @return a vector of entropy values
#' @export
calculate_rdi_corr <- function(id1, id2, expr1, expr2, delays, N_operations, verbose = F) {
	dmi <- c()

	for(d in 1:length(delays)){
		if(class(expr1) != 'list') dmi[d] <- rdi_single_run(expr1, expr2, d = delays[d])
		else{ #multiple branches
			expr1_t <- c(); expr2_t <- c(); expr2_past <- c()
			for(run_id in 1:nrow(expr1)){
                expr1_t <- cbind(expr1_t, expr1[run_id, 1:(ncol(expr1)-d)])
                expr2_t <- cbind(expr2_t, expr2[run_id, -c(1:d)])
                expr2_past <- cbind(expr2_past, expr2[run_id, d:(ncol(expr2) - 1)])
			}
			dmi[d] <- cmi(expr1_t, expr2_t, expr2_past)
		}
	}

	if(class(expr1) != 'list') correlation <- cor(matrix(expr1, ncol = 1), matrix(expr2, ncol = 1), method = 'spearman')
	else{
        expr1_concatenated <- c(); expr2_concatenated <- c()
        for(run_id in 1:nrow(expr1)) {
            expr1_concatenated <- cbind(expr1_concatenated, expr1[run_id, ])
            expr2_concatenated <- cbind(expr2_concatenated, expr2[run_id, ])
        }
        correlation = cor(matrix(expr1_concatenated, ncol = 1), matrix(expr2_concatenated, ncol = 1), method = 'spearman')
	}

	if(verbose)
	  print(c(paste(id1, id2, sep = '\t'), dmi, correlation))

	res <- as.data.frame(t(c(id1, id2, dmi, correlation)), stringsAsFactors = F)
	row.names(res) = paste(id1, id2, sep = '_')
	colnames(res) <- c('id_1', 'id_2', paste('delays ', delays, sep = ''), 'spearman_cor')
	res[-c(1:2)] <- as.numeric(res[-c(1:2)])
	return(res)
}

#' Calculate all pairwise RDI and spearman rank correlation
#'
#' This function estimates the RDI and the spearman correlation of expression (expr1, expr2) for any two genes, id1, di2 in the pseudo-time-series data
#' @param genes_data pseudo-time/time series for the gene expression data
#' @param delays a vector of delays to scanning the data
#' @param supergraph a dataframe storing the list of possible edges
#' @param cores number of cores used for performing the calculation
#' @param verbose a logic argument to determine whether or not the details of calculation should be reported
#' @return a dataframe storing the RDI results. The first two columns correspond to the ensemble id while the last column corresponds to the spearman
#' rank correlation. The rest columns corresponds to the RDI values under different delays
#' @export
calculate_and_write_pairwise_dmi <- function(genes_data, delays = c(1,2,5,10,15,20,25), supergraph = NULL, cores = 1, verbose = F){
	if(verbose) print("Calculating the Directed Mutual Information for each pair of genes...")

	cnt <- 0
	N_operations <- nrow(genes_data) * (nrow(genes_data) - 1)

	pairwise_dmi_results <- c()
	pairwise_corr_results <- c()

	tmp <- expand.grid(colnames(genes_data), colnames(genes_data), stringsAsFactors = F)
	all_pairwise_gene <- tmp[as.character(tmp[, 1]) != as.character(tmp[, 2]), ] #don't calculate RDI on self-loop

    if(!is.null(supergraph)) { #consider the supergraph
        all_pairwise_gene_paste <- paste(all_pairwise_gene[, 1], all_pairwise_gene[, 2])
        supergraph_paste <- paste(supergraph[, 1], supergraph[, 2])
        all_pairwise_gene <- all_pairwise_gene[all_pairwise_gene_paste %in% supergraph_paste, ]
    }

	all_pairwise_gene_list <- split(all_pairwise_gene, row.names(all_pairwise_gene))

	#we may need to convert genes_data into sparseMatrix before passing to mclapply
	res <- mclapply(all_pairwise_gene_list, function(x, genes_data, delays, N_operations) {
		calculate_rdi_corr(x[[1]], x[[2]], genes_data[, x[[1]]], genes_data[, x[[2]]], delays, N_operations)
		}, genes_data = genes_data, delays = delays, N_operations = N_operations, mc.cores = cores)

	res <- do.call(rbind.data.frame, res)
	row.names(res) <- paste(res$id_1, res$id_2, sep = '_')
	return(res)
}

#' Calculate conditionally RDI value
#'
#' This function estimates the RDI value for gene-pair combination (id1, id2) in the pseudotime/time series data, conditioned on the top K incoming edges.
#' Delay for the gene-pair under test corresponds to the delay with highest RDI value while delay for the conditioned k
#' incoming edges are based on their corresponding delay value.
#' @param id1 ensemble id for gene 1
#' @param id2 ensemble id for gene 2
#' @param genes_data pseudo-time/time series for the gene expression data
#' @param d delay for current gene interaction under test corresponding to the largest RDI value
#' @param top_k_genes_delay a vector of top k incoming node's delay corresponding to the maximal RDI value
#' @param N_operations number of total edges for calculating RDI values (not used for now)
#' @return a dataframe storing conditional RDI results. First two columns are the id names for the genes.
#' Third column is the conditional RDI value.
#' @export
calculate_conditioned_rdi <- function(id1, id2, genes_data, d, top_k_genes_delay, N_operations, verbose = F){
	if(class(genes_data) != 'list'){
	  top_k_genes_data <- data.frame()
		for(id in names(top_k_genes_delay)) {
			top_k_genes_data <- rbind(top_k_genes_data, id =  genes_data[, id])
		}
		conditioned_rdi <- rdi_single_run_conditioned(genes_data[, id1], genes_data[, id2], d = d, z = matrix(as.matrix(top_k_genes_data), nrow = nrow(genes_data)), z_delays = top_k_genes_delay
		                                                )
	}
	else{
        expr1_t = c(); expr2_t = c(); past = c()
        for(run_id in 1:length(genes_data)) {
            tau <- max(c(top_k_genes_delay, d))
            tot_len <- nrow(genes_data[run_id]) - tau
            expr1_t <- genes_data[run_id][id1, (tau - d + 1):(tau - d + tot_len)]
            expr2_t <- genes_data[run_id][id2, (tau + 1):(tau+tot_len)]
            past_tmp <- genes_data[run_id][id2, (tau):(tau-1+tot_len)]

            for(id in top_k_genes_delay){
                past_tmp <- cbind(past_tmp, genes_data[run_id][id, (tau - top_k_genes_delay[id] + 1):(tau - top_k_genes_delay[id]+tot_len)])
            }
            past <- cbind(past, past_tmp)
        }

        conditioned_rdi <- cmi(matrix(expr1_t, nrow = 1), matrix(expr2_t, nrow = 1), matrix(past, nrow = 1))
	}

  if(verbose == T)
	  print(c(paste(id, id2, sep = '\t'), conditioned_rdi))

	res <- as.data.frame(t(c(id1, id2, conditioned_rdi)))
	row.names(res) = paste(id1, id2, sep = '_')
	colnames(res) <- c('id_1', 'id_2', 'conditioned_rdi')
	return(res)
}

#' Calculate all pairwise conditional RDI values
#'
#' This function estimates the RDI value for every gene-pair combination conditioned on the top K incoming edges.
#' Delay for the gene-pair under test corresponds to the delay with highest RDI value while delay for the conditioned K
#' incoming edges are based on their corresponding delay value.
#' @param genes_data expression matrix for all the genes included in the RDI calculation
#' @param rdi_res calculated RDI results, including ids under test, RDI values under different delays and spearman rank corelation
#' @param k number of top k-incoming nodes for calculating RDI values. Default is 1. Higher k requres to have more samples or longer time / "pseudotime" series
#' @param supergraph a data frame stores the list of feasible edges for RDI analysis. Possible edges outside of this graph will be ignored
#' @param cores number of cores used in the calculation
#' @param verbose a logical argument to determine wheter or not to print the detailed running information.
#' @return a dataframe storing conditional RDI results. First two columns are the id names for the genes.
#' Third column is the conditional RDI value and last column is the time delay corresponding to the largest RDI.
#' @export
calculate_and_write_pairwise_dmi_conditioned <- function(genes_data, rdi_res, k = 1, supergraph = NULL, cores = 1, verbose = F){
	k <- min(k, nrow(genes_data) - 2) # not conditioned on itself and the already included incoming node(s)

    # delay for each gene pair gives max RDI value:
    rdi_res$max_rdi <- apply(rdi_res[, 3:(ncol(rdi_res) - 1)], 1, max)
    max_rdi <- apply(rdi_res[, 3:(ncol(rdi_res) - 2)], 1, which.max)
    delays_max <- as.numeric(do.call(rbind.data.frame, strsplit(colnames(rdi_res)[3:(ncol(rdi_res) - 2)], ' '))[[2]])[max_rdi]
    rdi_res$max_val_delay <- delays_max
    names(delays_max) <- row.names(rdi_res)

    # find the top k incoming genes with highest RDI values:
    target_gene_list <- unique(rdi_res[, 'id_2'])
    top_k_plus_1_incoming_id_list <- lapply(target_gene_list, function(x) {
        rdi_res_subset <- rdi_res[rdi_res$id_2 == x, ]
        nodes <- rdi_res_subset[order(rdi_res_subset$max_rdi, decreasing = T)[1:(k + 1)], 'id_1'] #incoming edge
        delay <- rdi_res_subset[order(rdi_res_subset$max_rdi, decreasing = T)[1:(k + 1)], 'max_val_delay'] #incoming edge
        names(delay) <- nodes
        return(delay)
        })

    names(top_k_plus_1_incoming_id_list) <- target_gene_list

    tmp <- expand.grid(colnames(genes_data), colnames(genes_data), stringsAsFactors = F)
    all_pairwise_gene <- tmp[as.character(tmp[, 1]) != as.character(tmp[, 2]), ] #don't calculate RDI on self-loop

    if(!is.null(supergraph)) { #consider the supergraph
        all_pairwise_gene_paste <- paste(all_pairwise_gene[, 1], all_pairwise_gene[, 2])
        supergraph_paste <- paste(supergraph[, 1], supergraph[, 2])
        all_pairwise_gene <- all_pairwise_gene[all_pairwise_gene_paste %in% supergraph_paste, ]
    }

    all_pairwise_gene_list <- split(all_pairwise_gene, row.names(all_pairwise_gene))

    # print(top_k_plus_1_incoming_id_list)

        #we may need to convert genes_data into sparseMatrix before passing to mclapply
	res <- mclapply(all_pairwise_gene_list, function(gene_pair, genes_data, delays, N_operations) {
        index_name <- paste(gene_pair[1], gene_pair[2], sep = '_')
        # if(index_name == "ENSMUSG00000020044.7_ENSMUSG00000013089.9"){
        #   browser()
        # }
        # if(index_name == 'ENSMUSG00000015452.8_ENSMUSG00000020044.7') browser()

    top_k_plus_1 <- top_k_plus_1_incoming_id_list[[gene_pair[[2]]]] #avoid taking the incoming node for current test as the conditional incoming node
    valid_top_k <- top_k_plus_1[!(names(top_k_plus_1) %in% gene_pair[[1]])][1:(length(top_k_plus_1) - 1)]

		calculate_conditioned_rdi(gene_pair[[1]], gene_pair[[2]], genes_data, delays_max[index_name], valid_top_k, N_operations)
		}, genes_data = genes_data, mc.cores = cores)

    res <- do.call(rbind.data.frame, res)
    row.names(res) <- paste(res$id_1, res$id_2, sep = '_')
    res$delay_max <- delays_max[row.names(res)]

    return(res)
}

#' This function smoothes the data using moving average
#'
#' @param exprs a dataframe or matrix to be smoothed, row corresponds variables while column correspond to sample
#' @param window_size size of sliding window
#' @return moving average smoothed dataset
#' @export
smooth_genes <- function(exprs, window_size = 40) {
  win_range <- nrow(exprs) - window_size
  exprs_smooth <- exprs[-c(1:window_size), ]

  res <- apply(exprs, 2, function(x) {
    tmp <- rep(0, win_range)
    for(i in 0:(win_range)){
      tmp[i + 1] <- mean(x[i + c(1:window_size)])
    }
    return(tmp)
  })

  return(res)
}

