# use loess to make the plot 

#' Plot the expression dynamics for a pair of genes in pseduotime with separate kinetic curves.
#'
#' This plotting function is used to make kinetic plots for a pair of genes along the pseudotime. It accepts a matrix where each 
#' row is the gene pair and the first column is the hypothetical source or regulator while the second column represents the hypothetical 
#' target. The name in this matrix should match the name in the gene_short_name column of the cds object. Two vertical lines will 
#' be added which indicate the estimated inflection points of the two kinetic curves. The distance between the two curves implies 
#' the potential time delay betweeen two genes. Since Scribe doesn't consider the continous value of the pseudotime but merely the 
#' ordering of the cells along the pseudotime, the x-axis of this plot also only includes the ordering information.  
#'
#' @param cds CellDataSet for the experiment
#' @param gene_pairs_mat A matrix where each row is the gene pair and the first column is the hypothetical source or regulator while the second column represents the hypothetical target. The name in this matrix should match the name in the gene_short_name column of the cds object. 
#' @param min_expr the minimum (untransformed) expression level to use in plotted the genes.
#' @param cell_size the size (in points) of each cell used in the plot
#' @param n_row number of columns used to layout the faceted cluster panels
#' @param n_col number of columns used to layout the faceted cluster panels
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to be used to color each cell
#' @param trend_formula the model formula to be used for fitting the expression trend over pseudotime
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param relative_expr A logic flag to determine whether or not we should use the relative expression values
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' gene_pairs_mat <- matrix(c('H19', 'Ccnd2', 'Ccnd2', 'Scnn1g'), ncol = 2)
#' plot_gene_pairs_in_pseudotime(lung, gene_pairs_mat)
#' }
plot_gene_pairs_in_pseudotime <- function(cds_subset, 
                              gene_pairs_mat, 
                              min_expr=NULL, 
                              cell_size=1.5, 
                              n_row=NULL, 
                              n_col=1, 
                              # panel_order=NULL, 
                              color_by="State",
                              trend_formula="~ sm.ns(Pseudotime, df=3)",
                              label_by_short_name=TRUE,
                              relative_expr=TRUE,
                              vertical_jitter=NULL,
                              horizontal_jitter=NULL) {

  all_genes_in_pair <- as.vector(gene_pairs_mat)
  if(! all(all_genes_in_pair %in% fData(cds_subset)$gene_short_name)) {
    stop("cds_subset doesn't include all genes in gene_pairs_mat Make sure all genes are included in gene_short_name column of the cds")
  }

  cds_subset <- cds_subset[row.names(subset(fData(cds_subset), gene_short_name %in% as.vector(gene_pairs_mat))), ]

  # for(gene_name in unique(row.names(cds_subset))) {
  #   data <- subset(pData(cds_subset), gene_id == gene_name)
  #   
  #   data$Pseudotime <- order(data$Pseudotime)
    pData(cds_subset)[, 'Pseudotime'] <- order(pData(cds_subset)[, 'Pseudotime'])
  # }
  
  f_id <- NA
  Cell <- NA
    if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
        integer_expression <- TRUE
    }
    else {
        integer_expression <- FALSE
        relative_expr <- TRUE
    }
    if (integer_expression) {
        cds_exprs <- exprs(cds_subset)
        if (relative_expr) {
            if (is.null(sizeFactors(cds_subset))) {
                stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
            }
            cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
        }
        cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    }
    else {
        cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
    }
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_pData <- pData(cds_subset)
    cds_fData <- fData(cds_subset)
    cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    #cds_exprs$f_id <- as.character(cds_exprs$f_id)
    #cds_exprs$Cell <- as.character(cds_exprs$Cell)
    
    if (integer_expression) {
        cds_exprs$adjusted_expression <- cds_exprs$expression
    }
    else {
        cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
    }
    # trend_formula <- paste("adjusted_expression", trend_formula,
    #     sep = "")
    if (label_by_short_name == TRUE) {
        if (is.null(cds_exprs$gene_short_name) == FALSE) {
            cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
            cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
        }
        else {
            cds_exprs$feature_label <- cds_exprs$f_id
        }
    }
    else {
        cds_exprs$feature_label <- cds_exprs$f_id
    }
    cds_exprs$f_id <- as.character(cds_exprs$f_id)
    cds_exprs$feature_label <- factor(cds_exprs$feature_label)

    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
    model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula,
                        relative_expr = T, new_data = new_data)
    colnames(model_expectation) <- colnames(cds_subset)
    expectation <- ddply(cds_exprs, .(f_id, Cell), function(x) data.frame("expectation"=model_expectation[x$f_id, x$Cell]))
    cds_exprs <- merge(cds_exprs, expectation)
    #cds_exprs$expectation <- expectation#apply(cds_exprs,1, function(x) model_expectation[x$f_id, x$Cell])

    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
    # if (is.null(panel_order) == FALSE) {
    #     cds_subset$feature_label <- factor(cds_subset$feature_label,
    #         levels = panel_order)
    # }

    # add inflection point here: 
    cds_exprs$inflection_point <- 0
    for(gene_name in unique(cds_exprs$f_id)) {
      data <- subset(cds_exprs, f_id == gene_name)
      
      message('gene_name is ', gene_name)
      inflection_point <- bede(data$Pseudotime, data$expectation, 0)
      if(!is.finite(inflection_point$iplast))
        inflection_point <- bede(data$Pseudotime, data$expectation, 1)
      # inflection_point <- tryCatch(bede(data$Pseudotime, data$expectation, 0), error = function() 
      #   message('Running inflection point detection in the mode of convex/concave failed, trying concave/convex mode instead'), 
      #   bede(data$Pseudotime, data$expectation, 1))
      message('inflection_point is ', inflection_point$iplast)
      cds_exprs[row.names(data), 'inflection_point'] <- inflection_point$iplast
    }
    
    # enable us to plot two genes together for visualizing RDI 
    
    cds_exprs$pairs <- "None"
    cds_exprs$source <- F
    for(i in 1:nrow(gene_pairs_mat)) {
      x <- gene_pairs_mat[i, ]
      current_gene_feature_ind <- which(cds_exprs$feature_label %in% x) #index: is belong to the current gene pair? 
      current_gene_pairs_bool <- cds_exprs$pairs[current_gene_feature_ind] == "None" # did the pairs column already set? 
      
      cds_exprs$pairs[current_gene_feature_ind[current_gene_pairs_bool]] <- paste0(x[1], ' -> ', x[2]) # now set the element in pairs column (currently NA values) 
      cds_exprs$source[current_gene_feature_ind[current_gene_pairs_bool]] <- cds_exprs$feature_label[current_gene_feature_ind[current_gene_pairs_bool]] == x[1] #   
      
      dup_cds_subset <- cds_exprs[current_gene_feature_ind[!current_gene_pairs_bool], ] # 
      
      if(nrow(dup_cds_subset) > 1)
      {  
        dup_cds_subset$pairs <- paste0(x[1], ' -> ', x[2]) # 
        dup_cds_subset$source[dup_cds_subset$feature_label == x[1]] <- T # 
        dup_cds_subset$source[dup_cds_subset$feature_label != x[1]] <- F # 
        
        cds_exprs <- rbind(cds_exprs, dup_cds_subset) # 
      }
    }
    
    q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
    if (is.null(color_by) == FALSE) {
      q <- q + geom_point(aes_string(color = color_by, shape = "source"), size = I(cell_size), position=position_jitter(horizontal_jitter, vertical_jitter))
    } else {
      q <- q + geom_point(size = I(cell_size), aes(shape = source), position=position_jitter(horizontal_jitter, vertical_jitter))
    }
    
    q <- q + geom_line(aes(x = Pseudotime, y = expectation, linetype = source), data = cds_exprs)
    
    q <- q + facet_wrap(~pairs, nrow = n_row,
                                          ncol = n_col, scales = "free_y") + scale_y_log10() 
    if (min_expr < 1) {
      q <- q + expand_limits(y = c(min_expr, 1))
    }
    if (relative_expr) {
      q <- q + ylab("Relative Expression")
    } else {
      q <- q + ylab("Absolute Expression")
    }
    q <- q + xlab("Pseudo-time")
    q <- q + monocle:::monocle_theme_opts() + geom_vline(aes(xintercept = inflection_point, linetype = source))
    q 
}

#' Plot the expression dynamics for a pair of genes in pseduotime with separate kinetic curves across two different lineages.
#'
#' This plotting function is used to make kinetic plots for a pair of genes along the pseudotime in each of branches defined with
#' the branch point. So each pair of genes has two panels each represents a different branch. This function accepts a matrix where each 
#' row is the gene pair and the first column is the hypothetical source or regulator while the second column represents the hypothetical 
#' target. The name in this matrix should match the name in the gene_short_name column of the cds object. Two vertical lines will 
#' be added which indicate the estimated inflection points of the two kinetic curves. The distance between the two curves implies 
#' the potential time delay betweeen two genes. Since Scribe doesn't consider the continous value of the pseudotime but merely the 
#' ordering of the cells along the pseudotime, the x-axis of this plot also only includes the ordering information.  
#' 
#' @param cds CellDataSet for the experiment
#' @param gene_pairs_mat A matrix where each row is the gene pair and the first column is the hypothetical source or regulator while the second column represents the hypothetical target. The name in this matrix should match the name in the gene_short_name column of the cds object. 
#' @param min_expr the minimum (untransformed) expression level to use in plotted the genes.
#' @param cell_size the size (in points) of each cell used in the plot
#' @param n_row number of columns used to layout the faceted cluster panels
#' @param n_col number of columns used to layout the faceted cluster panels
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to be used to color each cell
#' @param trend_formula the model formula to be used for fitting the expression trend over pseudotime
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param relative_expr A logic flag to determine whether or not we should use the relative expression values
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' gene_pairs_mat <- matrix(c('H19', 'Ccnd2', 'Ccnd2', 'Scnn1g'), ncol = 2)
#' plot_gene_pairs_branched_pseudotime(lung, gene_pairs_mat)
#' }
plot_gene_pairs_branched_pseudotime <- function(cds, 
                                    gene_pairs_mat = NULL,
                                    branch_states = NULL, 
                                    branch_point=1,
                                    branch_labels = NULL,
                                    method = "fitting", 
                                    min_expr = NULL, 
                                    cell_size = 0.75,
                                    n_row = NULL, 
                                    n_col = 2, 
                                    panel_order = NULL, 
                                    color_by = "State",
                                    expression_curve_linetype_by = "Branch", 
                                    trend_formula = "~ sm.ns(Pseudotime, df=3) * Branch", 
                                    reducedModelFormulaStr = NULL, 
                                    label_by_short_name = TRUE,
                                    relative_expr = TRUE,
                                    ...)
{

  all_genes_in_pair <- as.vector(gene_pairs_mat)
  if(! all(all_genes_in_pair %in% fData(cds)$gene_short_name)) {
    stop("cds doesn't include all genes in gene_pairs_mat Make sure all genes are included in gene_short_name column of the cds")
  }

  cds <- cds[row.names(subset(fData(cds), gene_short_name %in% as.vector(gene_pairs_mat))), ]

  Branch <- NA  
  if (is.null(reducedModelFormulaStr) == FALSE) {
    pval_df <- branchTest(cds, 
                          branch_states=branch_states,
                          branch_point=branch_point,
                          fullModelFormulaStr = trend_formula,
                          reducedModelFormulaStr = "~ sm.ns(Pseudotime, df=3)", 
                          ...)
    fData(cds)[, "pval"] <- pval_df[row.names(cds), 'pval']
  }
  if("Branch" %in% all.vars(terms(as.formula(trend_formula)))) { #only when Branch is in the model formula we will duplicate the "progenitor" cells
    cds_subset <- monocle::buildBranchCellDataSet(cds = cds, 
                                         branch_states = branch_states, 
                                         branch_point=branch_point,
                                         branch_labels = branch_labels, 
                                         progenitor_method = 'duplicate',
                                         ...)
  }
  else {
    cds_subset <- cds
    pData(cds_subset)$Branch <- pData(cds_subset)$State
  }
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    integer_expression <- TRUE
  }
  else {
    integer_expression <- FALSE
  }
  if (integer_expression) {
    CM <- exprs(cds_subset)
    if (relative_expr){
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      CM <- Matrix::t(Matrix::t(CM)/sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(as.matrix(CM)))
  }
  else {
    cds_exprs <- reshape2::melt(exprs(cds_subset))
  }
  if (is.null(min_expr)) {
    min_expr <- cds_subset@lowerDetectionLimit
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_pData <- pData(cds_subset)
  
  cds_fData <- fData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  if (integer_expression) {
    cds_exprs$adjusted_expression <- round(cds_exprs$expression)
  }
  else {
    cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  }
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$feature_label <- as.factor(cds_exprs$feature_label)
  # trend_formula <- paste("adjusted_expression", trend_formula,
  #     sep = "")
  cds_exprs$Branch <- as.factor(cds_exprs$Branch) 
  
  new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, Branch = pData(cds_subset)$Branch)
  
  full_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula, 
                                            relative_expr = T, new_data = new_data)
  colnames(full_model_expectation) <- colnames(cds_subset)
  
  cds_exprs$full_model_expectation <- apply(cds_exprs,1, function(x) full_model_expectation[x[2], x[1]])
  if(!is.null(reducedModelFormulaStr)){
    reduced_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = reducedModelFormulaStr,
                                                 relative_expr = T, new_data = new_data)
    colnames(reduced_model_expectation) <- colnames(cds_subset)
    cds_exprs$reduced_model_expectation <- apply(cds_exprs,1, function(x) reduced_model_expectation[x[2], x[1]])
  }
  
  # FIXME: If you want to show the bifurcation time for each gene, this function
  # should just compute it. Passing it in as a dataframe is just too complicated
  # and will be hard on the user. 
  # if(!is.null(bifurcation_time)){
  #     cds_exprs$bifurcation_time <- bifurcation_time[as.vector(cds_exprs$gene_short_name)]
  # }
  if (method == "loess")
    cds_exprs$expression <- cds_exprs$expression + cds@lowerDetectionLimit
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  if (is.null(panel_order) == FALSE) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label,
                                      levels = panel_order)
  }
  cds_exprs$expression[is.na(cds_exprs$expression)] <- min_expr
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$full_model_expectation[is.na(cds_exprs$full_model_expectation)] <- min_expr
  cds_exprs$full_model_expectation[cds_exprs$full_model_expectation < min_expr] <- min_expr
  
  if(!is.null(reducedModelFormulaStr)){
    cds_exprs$reduced_model_expectation[is.na(cds_exprs$reduced_model_expectation)] <- min_expr
    cds_exprs$reduced_model_expectation[cds_exprs$reduced_model_expectation < min_expr] <- min_expr
  }
  
  cds_exprs$State <- as.factor(cds_exprs$State)
  cds_exprs$Branch <- as.factor(cds_exprs$Branch)
  
  cds_exprs$pairs <- "None"
  cds_exprs$source <- F
  for(i in 1:nrow(gene_pairs_mat)) {
    x <- gene_pairs_mat[i, ]
    current_gene_feature_ind <- which(cds_exprs$feature_label %in% x) #index: is belong to the current gene pair? 
    current_gene_pairs_bool <- cds_exprs$pairs[current_gene_feature_ind] == "None" # did the pairs column already set? 
    
    cds_exprs$pairs[current_gene_feature_ind[current_gene_pairs_bool]] <- paste0(x[1], ' -> ', x[2]) # now set the element in pairs column (currently NA values) 
    cds_exprs$source[current_gene_feature_ind[current_gene_pairs_bool]] <- cds_exprs$feature_label[current_gene_feature_ind[current_gene_pairs_bool]] == x[1] #   
    
    dup_cds_subset <- cds_exprs[current_gene_feature_ind[!current_gene_pairs_bool], ] # 
    
    if(nrow(dup_cds_subset) > 1)
    {  
      dup_cds_subset$pairs <- paste0(x[1], ' -> ', x[2]) # 
      dup_cds_subset$source[dup_cds_subset$feature_label == x[1]] <- T # 
      dup_cds_subset$source[dup_cds_subset$feature_label != x[1]] <- F # 
      
      cds_exprs <- rbind(cds_exprs, dup_cds_subset) # 
    }
  }
  
  q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
  # if (!is.null(bifurcation_time)) {
  #   q <- q + geom_vline(aes(xintercept = bifurcation_time),
  #                       color = "black", linetype = "longdash")
  # }
  if (is.null(color_by) == FALSE) {
    q <- q + geom_point(aes_string(color = color_by, shape = "feature_label"), size = I(cell_size))
  }
  if (is.null(reducedModelFormulaStr) == FALSE) {
    q <- q + scale_y_log10() + facet_wrap(~pairs + Branch + 
                                            pval, nrow = n_row, ncol = n_col, scales = "free_y")
  } else q <- q + scale_y_log10() + facet_wrap(~pairs + Branch,
                                               nrow = n_row, ncol = n_col, scales = "free_y")
  if (method == "loess") {
    q <- q + stat_smooth(aes(fill = Branch, color = Branch),
                         method = "loess")
  } else if (method == "fitting") {
    q <- q + geom_line(aes_string(x = "Pseudotime", y = "full_model_expectation",
                                  linetype = "source"), data = cds_exprs) #+ scale_color_manual(name = "Type", values = c(colour_cell, colour), labels = c("Pre-branch", "AT1", "AT2", "AT1", "AT2")
  }
  
  if(!is.null(reducedModelFormulaStr)) {
    q <- q + geom_line(aes_string(x = "Pseudotime", y = "reduced_model_expectation"),
                       color = 'black', linetype = 2, data =  cds_exprs)   
  }
  
  q <- q + ylab("Expression") + xlab("Pseudotime (stretched)")
  
  q <- q + monocle:::monocle_theme_opts()
  q + expand_limits(y = min_expr)
  
}

# library(RColorBrewer)
#' 
#' Plot the scatterplot for two genes across different cells. 
#' 
#' This plot function is the commomly used approach to visualize the relationship between a gene pair.
#' Please check other plots for plot_lag_drevi, plot_cross_map, plot_rdi_pairs_heatmap, plot_comb_logic_heatmap,  
#' This function accepts a matrix where each row is the gene pair and the first column is the hypothetical source 
#' or regulator while the second column represents the hypothetical target. The name in this matrix should match 
#' the name in the gene_short_name column of the cds object.
#' 
#' @param cds CellDataSet for the experiment
#' @param gene_pairs_mat A matrix where each row is the gene pair and the first column is the hypothetical source or regulator while the second column represents the hypothetical target. The name in this matrix should match the name in the gene_short_name column of the cds object. 
#' @param n_row number of columns used to layout the faceted cluster panels
#' @param n_col number of columns used to layout the faceted cluster panels
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' gene_pairs_mat <- matrix(c('H19', 'Ccnd2', 'Ccnd2', 'Scnn1g'), ncol = 2)
#' plot_scatter_gene_pairs(lung, gene_pairs_mat)
#' }
#' @export
#' 
plot_scatter_gene_pairs <- function(cds,
                            gene_pairs_mat,
                            n_row = NULL,
                            n_col = 1) { #, h = 1
  
  gene_id <- row.names(subset(fData(cds), gene_short_name %in% gene_pairs_mat[i, ]))

  buylrd = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
             "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
  
  if(ncol(gene_pairs_mat) != 2) {
    stop('Please provide a matrix for pairs of genes you want to infer causality.')
  }
  if(length(gene_id) == 0) {
    stop('Please make sure all genes you provided in gene_pairs_mat exist in your data.')
  }

  if(nrow(gene_pairs_mat) == 1) {
    exprs_res <- t(exprs(cds)[gene_id, ])
  }
  else {
    exprs_res <- NULL
    for(i in 1:nrow(gene_pairs_mat)) {
      exprs_res <- cbind(exprs_res, exprs(cds)[gene_id, ])
    }
  }
  exprs_res_df <- as.data.frame(t(exprs_res))
  colnames(exprs_res_df) <- c('V1', 'V2')
  exprs_res_df$pair <- rep(apply(gene_pairs_mat, 1, function(x) paste(x[1], x[2], sep = "_")), each = ncol(cds))

  p <- qplot(V1, V2, data = exprs_res_df, geom = 'point', log = 'xy', size = I(1), size = 2) +
    stat_density2d(geom = "tile", aes(fill = ..density..), contour = FALSE) + xlab('') + ylab('') +
    geom_point(color = I("darkgray"), size = I(0.85)) + scale_fill_gradientn(colours = buylrd) +
    monocle:::monocle_theme_opts() + facet_wrap(~pair, nrow = n_row, ncol = n_col)

  return(p)
}

# plot_gene_pairs(lung, gene_pairs_mat = t(matrix(row.names(lung)[2:3])))
#
# mi(exprs(lung)[2, ], exprs(lung)[3, ], k = 5)
# test <- rbind(t(matrix(row.names(lung)[2:3])), t(matrix(row.names(lung)[3:4])))
# plot_gene_pairs(lung, gene_pairs_mat = test)
# mi(exprs(lung)[2, ], exprs(lung)[3, ], k = 5)

# #add the principial curves:
# # ini_path <- exprs_data[sort(pData(data)$Pseudotime, index.return = T)$ix, ] #from the MST/pc tree
# # ini_path <- loess(exprs_data[, 1] ~ exprs_data[, 2])#from the lowess curve
# ini_path <- lowess(exprs_data[, 1], exprs_data[, 2]); ini_path <- as.matrix(data.frame(x = ini_path$x, y = ini_path$y))
#
# gene_pair_pc <- principal.curve(as.matrix(exprs_data), start = ini_path)

# geom_path(aes(gene_pair_pc$s[gene_pair_pc$tag, 1], gene_pair_pc$s[gene_pair_pc$tag, 2]), color = 'black') +
# theme(strip.background = element_rect(colour = "white",
#                                       fill = "white")) + theme(panel.border = element_blank(),
#                                                                axis.line = element_line()) + theme(legend.position = "none") +
#   theme(axis.title.y = element_text(size = 10)) + theme(axis.title.x = element_text(size = 10)) +
#   theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
#   theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())
# 
# # make MI, CMI plot: (use the shiny app?)
# 

# implement the support for multiple runs (run_vec) 
#' Plot the scatterplot formed by the space of x_{t - \mu}, y_t as well as y_{t - 1}. 
#' 
#' This plot function is used to visualize the space formed by x_{t - \mu}, y_t as well as y_{t - 1} which is 
#' used to calculate the RDI score I(x_t -> Y_t) = CMI(x_{t - \mu}, y_t | y_{t - 1}). Please check other ploting
#' functions (plot_lag_drevi, plot_cross_map, plot_rdi_pairs_heatmap, plot_comb_logic_heatmap) for more intuitive 
#' visualization. 
#' 
#' @param cds CellDataSet for the experiment
#' @param gene_pair A character vector consists of two gene short names while the first name is the hypothetical source and the second is the target. 
#' @param d The time delay between the source and the target gene 
#' @param run_vec A numeric vector of length equal to number of columns (cells) in the cds, encoding the runs (or replicate) id for each experiment. 
#' @return a plotly object
#' @import plotly
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' gene_pairs_mat <- matrix(c('H19', 'Ccnd2', 'Ccnd2', 'Scnn1g'), ncol = 2)
#' plot_rdi_scatter_gene_pairs(lung, gene_pairs_mat[1, ], d = 5)
#' }
#' @export
plot_rdi_scatter_gene_pairs <- function(cds, gene_pair, d = 1, run_vec = NULL){
  gene_id_1 <- row.names(subset(fData(cds), gene_short_name %in% gene_pair[1]))
  gene_id_2 <- row.names(subset(fData(cds), gene_short_name %in% gene_pair[2]))
  
  if(length(gene_id_1) == 0 | length(gene_id_2) == 0)
    stop('Please make sure all genes you provided in gene_pairs_mat exist in your data.')
  
  x <- exprs(cds)[gene_id_1, ]
  y <- exprs(cds)[gene_id_2, ]
  
  if (is.numeric(x))
    x <- as.matrix(x)
  if (is.numeric(y))
    y <- as.matrix(y)
  if (nrow(x) != nrow(y))
    stop("The number of time samples has to be the same for X and Y")
  
  plot_ly(type = 'scatter3d', x = x[1:(nrow(x) - d), ], y = y[-(1:d), ], z = y[d:(nrow(y) - 1), ], mode = 'markers')
}

# This function need to rewrite which should starts from calculatin the temporal RDI before creating the heatmap plot 
#' Plot the scatterplot for two genes across different cells. 
#' 
#' This plot function is the commomly used approach to visualize the relationship between a gene pair.
#' Please check other plots for plot_lag_drevi, plot_cross_map, plot_rdi_pairs_heatmap, plot_comb_logic_heatmap,  
#' This function accepts a matrix where each row is the gene pair and the first column is the hypothetical source 
#' or regulator while the second column represents the hypothetical target. The name in this matrix should match 
#' the name in the gene_short_name column of the cds object.
#' 
#' @param cds CellDataSet for the experiment
#' @param gene_pairs_mat A matrix where each row is the gene pair and the first column is the hypothetical source or regulator while the second column represents the hypothetical target. The name in this matrix should match the name in the gene_short_name column of the cds object. 
#' @param n_row number of columns used to layout the faceted cluster panels
#' @param n_col number of columns used to layout the faceted cluster panels
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
plot_regulation_over_time <- function(res, gene_name_vec) {
  dim(res) <- c(dim(res)[1], dim(res)[2] * dim(res)[2])
  
  all_cmbns <- expand.grid(gene_name_vec, gene_name_vec)
  valid_all_cmbns_df <- data.frame(pair = paste(tolower(all_cmbns$Var1), tolower(all_cmbns$Var2), sep = '_'), pval = 0)
  row.names(valid_all_cmbns_df) <- valid_all_cmbns_df$pair
  colnames(res) <- valid_all_cmbns_df$pair
  res <- apply(res, 1, function(x) (x - min(x)) / (max(x) - min(x)) ) # normalize by the maximal value
  pheatmap::pheatmap(res, cluster_rows = F, cluster_cols = T, annotation_names_col = T)
}

# sparse network visualization (Hive / cico plot )
# http://www.vesnam.com/Rblog/viznets3/
# https://en.wikipedia.org/wiki/Gephi
# end the paper with a hive plot (showing different network for all blood different cell type)  
# hierarchical edge bundling
# the layout from ISB 
# https://bost.ocks.org/mike/miserables/
# BioFabric
# http://www.diva-portal.org/smash/get/diva2:1034015/FULLTEXT01.pdf

#' Plot the lagged DREVI plot for pairs of genes across pseudotime. 
#' 
#' This plotting function builds on the original idea of DREVI plot but is extended in the context for causal network.  
#' It considers the time delay between the hypothetical regulators to the target genes which is parametered by `d`. 
#' Lagged DREVI plot first estimate the joint density (P(x_{t - d}, y_t)) for x_{t - d} and y_t and divide that by the 
#' marginal density P(x_{t - d}) to get the conditional density estimate (P(x_{t - d}, y_t | x_{x - d})). 
#' The 2-d density is estimated through the kde2d function from MASS package. We then calculate 
#' the z-score for each column of conditional density. Note that this plot tries to demonstrate the potential influence
#' between two variables instead of the factual influence. 
#' 
#' @param cds CellDataSet for the experiment
#' @param gene_pairs_mat A matrix where each row is the gene pair and the first column is the hypothetical source or regulator while the second column represents the hypothetical target. The name in this matrix should match the name in the gene_short_name column of the cds object. 
#' @param log A logic argument used to determine whether or not you should perform log transformation (using log(expression + 1)) before calculating density estimates, default to be TRUE. 
#' @param d The time delay between the source and target gene. 
#' @param grids The number of grid when creating the lagged DREVI plot.
#' @param n_row number of columns used to layout the faceted cluster panels
#' @param n_col number of columns used to layout the faceted cluster panels
#' @param scales The character string passed to facet function, determines whether or not the scale is fixed or free in different dimensions.
#' @param verbose A logic argument to determine whether or not we should print the detailed running information. 
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' gene_pairs_mat <- matrix(c('H19', 'Ccnd2', 'Ccnd2', 'Scnn1g'), ncol = 2)
#' plot_lag_drevi(lung, gene_pairs_mat)
#' }
#' @references Krishnaswamy, Smita, et al. "Conditional density-based analysis of T cell signaling in single-cell data." Science 346.6213 (2014): 1250689.
#' @export
#' 
plot_lag_drevi <- function(cds_subset, gene_pairs_mat, 
                           log = TRUE,
                           # nConBins = 4, 
                           d = 1,
                           grids = NULL,
                           n_row = NULL,
                           n_col = 1,
                           scales = "free",
                           verbose = FALSE) {

  gene_pairs_mat <- as.matrix(gene_pairs_mat)
  all_genes_in_pair <- as.vector(gene_pairs_mat)
  if(! all(all_genes_in_pair %in% fData(cds_subset)$gene_short_name)) {
    stop("cds_subset doesn't include all genes in gene_pairs_mat Make sure all genes are included in gene_short_name column of the cds")
  }

  cds_subset <- cds_subset[row.names(subset(fData(cds_subset), gene_short_name %in% all_genes_in_pair)), ]
  
  if(is.null(grids))
    dim_val <- (round((ncol(cds_subset) - d) / 4))
  else
    dim_val <- grids
  
  flat_res <- as.data.frame(matrix(0, nrow = dim_val^2  * nrow(gene_pairs_mat), ncol = 4))
  ridge_curve <- as.data.frame(matrix(0, nrow = dim_val * nrow(gene_pairs_mat), ncol = 3))
  colnames(flat_res) <- c("x", "y", "den", "type")
  colnames(ridge_curve) <- c("x", "y", "type")
  xy <- data.frame()

  id <- 0
  for(gene_pairs_ind in 1:nrow(gene_pairs_mat))  {
    # gene_pairs_ind <- 1
    if(verbose)
      message("current gene pair is ", gene_pairs_mat[gene_pairs_ind, 1], " -> ",  gene_pairs_mat[gene_pairs_ind, 2])
    
    gene_pairs <- gene_pairs_mat[gene_pairs_ind, ]
    f_ini_ind <- dim_val^2 * id #flat_res (normalized density results)
    r_ini_ind <- dim_val * id #ridge_curve 
    
    gene_pair_name <- paste(gene_pairs[1], gene_pairs[2], sep = ' -> ')

    x <- matrix(exprs(cds_subset)[row.names(subset(fData(cds_subset), gene_short_name %in% gene_pairs[1])), ], ncol = 1)
    y_ori <- matrix(exprs(cds_subset)[row.names(subset(fData(cds_subset), gene_short_name %in% gene_pairs[2])), ], ncol = 1)
    
    if(log) {
      x <- log(x + 1)
      y_ori <- log(y_ori + 1)
    }
    
    if(d != 0) {
      x  <- matrix(x[1:(nrow(x) - d), ], ncol = 1)
      y <- matrix(y_ori[-(1:d), ], ncol = 1)
      z <- y_ori[d:(nrow(y_ori) - 1), ]
    }
    else {
      x  <- matrix(x, ncol = 1)
      y <- matrix(y_ori, ncol = 1)
      z <- y_ori
    }
    
    exprs_res <- expression(paste("Target (", "y", phantom()[{
      paste("t")
    }], ")", ""))
    
    if(length(unique(x)) < dim_val | length(unique(y)) < dim_val) {
      # stop(paste0("Genes ", gene_pairs, "are required to express in at least ", dim_val, " cells"))
      # x <- x + rnorm(length(x), sd = 1e-10)
    }
    # # avoid the Error in MASS::kde2d(x, y, n = c(dim_val, dim_val), lims = c(min(x), max(x),:   bandwidths must be strictly positive
    # if(length(unique(x)) < dim_val | length(unique(y)) < dim_val) {
    #   stop(paste0("Genes ", gene_pairs, "are required to express in at least ", dim_val, " cells"))
    #   x <- x + rnorm(length(x), sd = 1e-10)
    # }
    # if(length(unique(y)) < dim_val)
    #   y <- y + rnorm(length(y), sd = 1e-10)
    # 
    ########################################################################################################################################################################
    # if(conditioning == T) {
    #   # # do a linear line fitting
    #   # df <- data.frame(y = y, z = z)
    #   # full_model_fit <- VGAM::vglm(as.formula("y~z"), data = df, family=gaussianff())
    #   # 
    #   # y <- resid(full_model_fit)
    #   # # 
    #   #       data <- Reduce(cbind, list(x, y, z))
    #   # 
    #   #       xyz_kde <- kde(data, gridsize = c(25,25, 25))
    #   # 
    #   #       x_meshgrid <- xyz_kde$eval.points[[1]]
    #   #       y_meshgrid <- xyz_kde$eval.points[[2]]
    #   #       z_meshgrid <- xyz_kde$eval.points[[3]]
    #   # 
    #   #       den_res <- matrix(0, nrow = length(x_meshgrid), ncol = length(y_meshgrid))
    #   # 
    #   #       nConBins <- length(z_meshgrid)
    #   #       for(conBins in 1:nConBins) {
    #   #         den_res_tmp <- xyz_kde$estimate[, , conBins]
    #   #         den_res_tmp[!is.finite(den_res_tmp)] <- 0
    #   #         den_res <- den_res + den_res_tmp / (nConBins * sum(den_res_tmp))
    #   #       }
    #   # 
    #   #       den_x <- colSums(den_res) # just calculate the sum for each column
    # }
    # ########################################################################################################################################################################
    # if(conditioning) {
    #   # exprs_res <- expression(paste("Target (", "y", phantom()[{
    #   #   paste("t")
    #   # }], "", "|", "y", phantom()[{
    #   #   paste("t", phantom() - phantom(), "1")
    #   # }], ")", ""))
    #   # 
    #   # rng_vec <- sort(z)[quantile(1:length(z), seq(0, 1, length.out = nConBins + 1))] #ensure each bin gets the same number of cells 
    #   # den_res_array <- array(0, dim = c(dim_val, dim_val, nConBins))
    #   # den_res <- matrix(0, nrow = dim_val, ncol = dim_val)
    #   # 
    #   # x_meshgrid <- c()
    #   # y_meshgrid <- c()
    #   # 
    #   # for(conBins in 2:(nConBins + 1)) {
    #   #   if(conBins != conBins) {
    #   #     x_rng <- x[z >= rng_vec[conBins - 1] & z < rng_vec[conBins]]
    #   #     y_rng <- y[z >= rng_vec[conBins - 1] & z < rng_vec[conBins]]
    #   #   }
    #   #   else {
    #   #     x_rng <- x[z >= rng_vec[conBins - 1] & z <= rng_vec[conBins]]
    #   #     y_rng <- y[z >= rng_vec[conBins - 1] & z <= rng_vec[conBins]]
    #   #   }
    #   #   
    #   #   bandwidth <- c(MASS::bandwidth.nrd(x_rng), MASS::bandwidth.nrd(y_rng))
    #   #   if(any(bandwidth == 0)) {
    #   #     max_vec <- c(max(x_rng), max(y_rng))
    #   #     bandwidth[bandwidth == 0] <- max_vec[bandwidth == 0] / dim_val
    #   #   }
    #   #   den_xy <- MASS::kde2d(x_rng, y_rng, n = c(dim_val, dim_val), lims = c(min(x), max(x), min(y), max(y)), h = bandwidth)
    #   #   den_res_array[, , conBins - 1] <- den_xy$z# as.data.frame()
    #   #   # dimnames(den_res) <- list(paste0("x_", as.character(den_xy$x)), paste0("y_", as.character(den_xy$y)))
    #   #   # den_x_res <- density(x, n = round(length(x))/ 4, from = min(x), to = max(x))
    #   #   # den_x <- den_x_res$y
    #   #   
    #   #   x_meshgrid <- den_xy$x
    #   #   y_meshgrid <- den_xy$y
    #   #   
    #   #   # message('x_meshgrid is ', x_meshgrid)
    #   #   # message('y_meshgrid is ', y_meshgrid)
    #   # }
    #   # max_ind <- 0
    #   # tmp <- 0
    #   # for(conBins in 1:(nConBins - 1)) {
    #   #   if(tmp < sum(den_res_array[, , conBins]))
    #   #     max_ind <- conBins
    #   #   
    #   #   tmp <- sum(den_res_array[, , conBins])
    #   #   # den_res_tmp <- den_res_array[, , conBins]
    #   #   # den_res_tmp[!is.finite(den_res_tmp)] <- 0
    #   #   # den_res <- den_res + den_res_tmp / nConBins
    #   #   # den_res <- den_res + den_res_tmp / (nConBins * sum(den_res_tmp)) 
    #   # }
    #   # den_res <- den_res_array[, , max_ind]
    #   # 
    #   # den_x <- rowSums(den_res) # just calculate the sum for each column 
    # } else {
      bandwidth <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))
      if(any(bandwidth == 0)) {
        max_vec <- c(max(x), max(y))
        bandwidth[bandwidth == 0] <- max_vec[bandwidth == 0] / dim_val
      }
      den_xy <- MASS::kde2d(x, y, n = c(dim_val, dim_val), lims = c(min(x), max(x), min(y), max(y)), h = bandwidth)
      den_res <- as.data.frame(den_xy$z)
      # dimnames(den_res) <- list(paste0("x_", as.character(den_xy$x)), paste0("y_", as.character(den_xy$y)))
      # den_x_res <- density(x, n = round(length(x))/ 4, from = min(x), to = max(x))
      # den_x <- den_x_res$y
      den_x <- rowSums(den_res) # just calculate the sum for each column

      x_meshgrid <- den_xy$x
      y_meshgrid <- den_xy$y
    # }
    
    max_ind <- 1
    
    for(i in 1:length(x_meshgrid)) {
      max_val <- max(den_res[i, ] / den_x[i]); min_val <- 0 #min(den_res[i, ] / den_x[i]) # 

      # print(den_x)
      if(den_x[i] == 0) {
      }
      else {
        max_ind <- which.max(den_res[i, ] / den_x[i]) 
      }
       # message("i is ", i, "j is ", j, "vector is ", c(x_meshgrid[i], y_meshgrid[max_ind], gene_pair_name))
      ridge_curve[i + r_ini_ind, ] <- c(x_meshgrid[i], y_meshgrid[max_ind], gene_pair_name)
      
      for(j in 1:length(y_meshgrid)) {
        rescaled_val <- (den_res[i, j] / den_x[i] - min_val) / (max_val - min_val)
        flat_res[(i - 1) * length(x_meshgrid) + j + f_ini_ind, ] <- c(x_meshgrid[i], y_meshgrid[j], rescaled_val, gene_pair_name)
      }
    }
    xy_tmp <- data.frame(x = x, y = y, 'type' = gene_pair_name)
    xy <- rbind(xy, xy_tmp)    
    
    id <- id + 1
  }

  # ridge_curve
  colnames(xy) <- c('x', 'y', 'type')
  flat_res[, 1:3] <- as.matrix(flat_res[, 1:3])
  ggplot(aes(as.numeric(x), as.numeric(y)), data = flat_res) +  geom_raster(aes(fill = as.numeric(den))) + 
    scale_fill_gradientn("Density", colours = terrain.colors(10)) + 
    geom_rug(aes(as.numeric(x), as.numeric(y)), data = xy, col="darkred",alpha=.1) + 
    geom_path(aes(as.numeric(x), as.numeric(y)), data  = ridge_curve, color = 'red') + facet_wrap(~type, scales = scales, nrow = n_row, ncol = n_col) + 
    xlab( expression(paste("Source (", "x", phantom()[{
      paste("t", phantom() - phantom(), d)
    }], ")", ""))) + 
    ylab(exprs_res) + monocle:::monocle_theme_opts()
}

# 
# lung <- load_lung()
# plot_rdi_drevi(log2(exprs(lung)[3, ] + 1), log2(exprs(lung)[27, ] + 1), x_label = as.character(fData(lung)$gene_short_name[3]), y_label = as.character(fData(lung)$gene_short_name[27]))
# 
# lung_gene_vec <- c("Cdk1", "Pdpn", "Ccnd2",  "Ccnd1")
# lung_gene_vec_ensemble <- row.names(subset(fData(lung), gene_short_name %in% lung_gene_vec))
# 
# gene_pairs_mat <- matrix(lung_gene_vec, nrow = 2)
# plot_rdi_pairs(lung, t(gene_pairs_mat), d = 20)

# LPS_gene_vec <- c("Stat1", "Hhex", "Stat2", "Rbl1")
# gene_pairs_mat <- matrix(LPS_gene_vec, nrow = 2)
# plot_rdi_pairs(Shalek_abs_subset_ko_LPS_new, t(gene_pairs_mat), d = 20)
# # 
# load('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/RData/analysis_scRNA_seq_Olsson.RData')
# Olsson_gene_vec <- c('Irf8', "Gfi1", 'Zeb2', 'Ets1', 'Irf5', "Irf8")
# 
# Olsson_gene_vec <- c('Zeb2', "Ets1")
# gene_pairs_mat <- matrix(Olsson_gene_vec, nrow = 2)
# plot_rdi_pairs(Olsson_granulocyte_cds, t(gene_pairs_mat), d = 1)

#' Plot the shadow mainfold from the lagged coordinates from a single variable. 
#' 
#' This plot function uses time lagged coordinates of a single variable (gene) x, (x_{t}, x_{t - d}, x_{t - 2d}) to reconstruct a shadow manifold. This manifold can 
#' be used for detecting causality.  
#' 
#' @param cds CellDataSet for the experiment
#' @param gene_name The name for gene used to reconstruct the shadow manifold. The gene name should appear in the gene_short_name column of the cds object. 
#' @param d The time delay between the source and target gene. 
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' x_1 <- all_cell_simulation[1, , 1]
#' x_2 <- all_cell_simulation[2, , 1]
#' x_3 <- all_cell_simulation[3, , 1]
#' 
#' plot_ly(type = 'scatter3d', x = x_1, y = x_2, z = x_3, mode = 'markers', color = 1:length(x_1))
#' lung <- load_lung() 
#' gene_pairs_mat <- matrix(c('H19', 'Ccnd2', 'Ccnd2', 'Scnn1g'), ncol = 2)
#' plot_shadow_manifold(lung, gene_pairs_mat[1])
#' }
#' @references Sugihara, George, et al. "Detecting causality in complex ecosystems." science 338.6106 (2012): 496-500.
#' @export
plot_shadow_manifold <- function(cds, gene_name, d = 1){
  gene_id <- row.names(subset(fData(cds), gene_short_name %in% gene_name))

  if(length(gene_id) == 0)
    stop('Please make sure the gene name you provided in gene_pairs_mat exist in your data.')
  
  if (is.numeric(x))
    x <- as.matrix(x)
  
  x_1 <- x[seq(1, length(x) - 2 * d, by = 1)]
  x_2 <- x[seq(1 + d, length(x) - d, by = 1)]
  x_3 <- x[seq(1 + 2 * d, length(x), by = 1)]
  # return(cmi(x[1:(nrow(x) - d), ], y[-(1:d), ], y[d:(nrow(y) - 1), ]))

  plot_ly(type = 'scatter3d', x = x_1, y = x_2, z = x_3, color = 1:length(x_1), mode = 'markers')
  
}

plot_gene_pair_delay <- function() {
  
}

# 
#' Plot the taylor diagram for the RDI, mutual information values. 
#' 
#' This plot function is the commomly used approach to visualize the relationship between a gene pair.
#' Please check other plots for plot_lag_drevi, plot_cross_map, plot_rdi_pairs_heatmap, plot_comb_logic_heatmap,  
#' This function accepts a matrix where each row is the gene pair and the first column is the hypothetical source 
#' or regulator while the second column represents the hypothetical target. The name in this matrix should match 
#' the name in the gene_short_name column of the cds object.
#' 
#' @param cds CellDataSet for the experiment
#' @param gene_pairs_mat A matrix where each row is the gene pair and the first column is the hypothetical source or regulator while the second column represents the hypothetical target. The name in this matrix should match the name in the gene_short_name column of the cds object. 
#' @param n_row number of columns used to layout the faceted cluster panels
#' @param n_col number of columns used to layout the faceted cluster panels
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' gene_pairs_mat <- matrix(c('H19', 'Ccnd2', 'Ccnd2', 'Scnn1g'), ncol = 2)
#' plot_taylor_diagram(cds_exprs[1, ], seq(0, 1, length.out = 100), MI = seq(2, 3, length.out = 100), Entropy = seq(4, 5, length.out = 1))
#' }
#' @references 
#' @export
#' 
# this function modifies the taylor.diagram function from plotrix package 
plot_taylor_diagram <- function (ref, RDI, MI, Reference, add = FALSE, col = "red", pch = 19, pos.cor = TRUE, 
                            xlab = "Entropy", ylab = "MI", main = "RDI Taylor Diagram", show.gamma = TRUE, 
                            ngamma = 3, gamma.col = 8, sd.arcs = 0, ref.sd = FALSE, sd.method = "sample", 
                            grad.corr.lines = c(0.2, 0.4, 0.6, 0.8, 0.9), pcex = 1, cex.axis = 1, 
                            normalize = FALSE, mar = c(5, 4, 6, 6), text = NULL, ...) 
{
  grad.corr.full <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 
                      1)
  
  R <- RDI
  sd.f <- MI
  sd.r <- Reference 
  
  # R <- cor(ref, model, use = "pairwise")
  # if (is.list(ref)) 
  #   ref <- unlist(ref)
  # if (is.list(model)) 
  #   ref <- unlist(model)
  # SD <- function(x, subn) {
  #   meanx <- mean(x, na.rm = TRUE)
  #   devx <- x - meanx
  #   ssd <- sqrt(sum(devx * devx, na.rm = TRUE)/(length(x[!is.na(x)]) - 
  #                                                 subn))
  #   return(ssd)
  # }
  # subn <- sd.method != "sample"
  # sd.r <- SD(ref, subn)
  # sd.f <- SD(model, subn)
  # if (normalize) {
  #   sd.f <- sd.f/sd.r
  #   sd.r <- 1
  # }
  maxsd <- 1.5 * max(sd.f, sd.r)
  oldpar <- par("mar", "xpd", "xaxs", "yaxs")
  if (!add) {
    if (pos.cor) {
      if (nchar(ylab) == 0) 
        ylab = "Standard deviation"
      par(mar = mar)
      plot(0, xlim = c(0, maxsd), ylim = c(0, maxsd), xaxs = "i", 
           yaxs = "i", axes = FALSE, main = main, xlab = xlab, 
           ylab = ylab, type = "n", cex = cex.axis, ...)
      if (grad.corr.lines[1]) {
        for (gcl in grad.corr.lines) lines(c(0, maxsd * 
                                               gcl), c(0, maxsd * sqrt(1 - gcl^2)), lty = 3)
      }
      segments(c(0, 0), c(0, 0), c(0, maxsd), c(maxsd, 
                                                0))
      axis.ticks <- pretty(c(0, maxsd))
      axis.ticks <- axis.ticks[axis.ticks <= maxsd]
      axis(1, at = axis.ticks, cex.axis = cex.axis)
      axis(2, at = axis.ticks, cex.axis = cex.axis)
      if (sd.arcs[1]) {
        if (length(sd.arcs) == 1) 
          sd.arcs <- axis.ticks
        for (sdarc in sd.arcs) {
          xcurve <- cos(seq(0, pi/2, by = 0.03)) * sdarc
          ycurve <- sin(seq(0, pi/2, by = 0.03)) * sdarc
          lines(xcurve, ycurve, col = "blue", lty = 3)
        }
      }
      if (show.gamma[1]) {
        if (length(show.gamma) > 1) 
          gamma <- show.gamma
        else gamma <- pretty(c(0, maxsd), n = ngamma)[-1]
        if (gamma[length(gamma)] > maxsd) 
          gamma <- gamma[-length(gamma)]
        labelpos <- seq(45, 70, length.out = length(gamma))
        for (gindex in 1:length(gamma)) {
          xcurve <- cos(seq(0, pi, by = 0.03)) * gamma[gindex] + 
            sd.r
          endcurve <- which(xcurve < 0)
          endcurve <- ifelse(length(endcurve), min(endcurve) - 
                               1, 105)
          ycurve <- sin(seq(0, pi, by = 0.03)) * gamma[gindex]
          maxcurve <- xcurve * xcurve + ycurve * ycurve
          startcurve <- which(maxcurve > maxsd * maxsd)
          startcurve <- ifelse(length(startcurve), max(startcurve) + 
                                 1, 0)
          lines(xcurve[startcurve:endcurve], ycurve[startcurve:endcurve], 
                col = gamma.col)
          if (xcurve[labelpos[gindex]] > 0) 
            boxed.labels(xcurve[labelpos[gindex]], ycurve[labelpos[gindex]], 
                         gamma[gindex], border = FALSE)
        }
      }
      xcurve <- cos(seq(0, pi/2, by = 0.01)) * maxsd
      ycurve <- sin(seq(0, pi/2, by = 0.01)) * maxsd
      lines(xcurve, ycurve)
      bigtickangles <- acos(seq(0.1, 0.9, by = 0.1))
      medtickangles <- acos(seq(0.05, 0.95, by = 0.1))
      smltickangles <- acos(seq(0.91, 0.99, by = 0.01))
      segments(cos(bigtickangles) * maxsd, sin(bigtickangles) * 
                 maxsd, cos(bigtickangles) * 0.97 * maxsd, sin(bigtickangles) * 
                 0.97 * maxsd)
      par(xpd = TRUE)
      if (ref.sd) {
        xcurve <- cos(seq(0, pi/2, by = 0.01)) * sd.r
        ycurve <- sin(seq(0, pi/2, by = 0.01)) * sd.r
        lines(xcurve, ycurve)
      }
      points(sd.r, 0, cex = pcex)
      text(cos(c(bigtickangles, acos(c(0.95, 0.99)))) * 
             1.05 * maxsd, sin(c(bigtickangles, acos(c(0.95, 
                                                       0.99)))) * 1.05 * maxsd, c(seq(0.1, 0.9, by = 0.1), 
                                                                                  0.95, 0.99))
      text(maxsd * 0.8, maxsd * 0.8, "RDI", srt = 315)
      segments(cos(medtickangles) * maxsd, sin(medtickangles) * 
                 maxsd, cos(medtickangles) * 0.98 * maxsd, sin(medtickangles) * 
                 0.98 * maxsd)
      segments(cos(smltickangles) * maxsd, sin(smltickangles) * 
                 maxsd, cos(smltickangles) * 0.99 * maxsd, sin(smltickangles) * 
                 0.99 * maxsd)
    }
    else {
      # x <- ref
      # y <- model
      # R <- cor(x, y, use = "pairwise.complete.obs")
      # E <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
      # xprime <- x - mean(x, na.rm = TRUE)
      # yprime <- y - mean(y, na.rm = TRUE)
      # sumofsquares <- (xprime - yprime)^2
      # Eprime <- sqrt(sum(sumofsquares)/length(complete.cases(x)))
      # E2 <- E^2 + Eprime^2
      # if (add == FALSE) {
      #   maxray <- 1.5 * max(sd.f, sd.r)
      #   plot(c(-maxray, maxray), c(0, maxray), type = "n", 
      #        asp = 1, bty = "n", xaxt = "n", yaxt = "n", 
      #        xlab = xlab, ylab = ylab, main = main, cex = cex.axis)
      #   discrete <- seq(180, 0, by = -1)
      #   listepoints <- NULL
      #   for (i in discrete) {
      #     listepoints <- cbind(listepoints, maxray * 
      #                            cos(i * pi/180), maxray * sin(i * pi/180))
      #   }
      #   listepoints <- matrix(listepoints, 2, length(listepoints)/2)
      #   listepoints <- t(listepoints)
      #   lines(listepoints[, 1], listepoints[, 2])
      #   lines(c(-maxray, maxray), c(0, 0))
      #   lines(c(0, 0), c(0, maxray))
      #   for (i in grad.corr.lines) {
      #     lines(c(0, maxray * i), c(0, maxray * sqrt(1 - 
      #                                                  i^2)), lty = 3)
      #     lines(c(0, -maxray * i), c(0, maxray * sqrt(1 - 
      #                                                   i^2)), lty = 3)
      #   }
      #   for (i in grad.corr.full) {
      #     text(1.05 * maxray * i, 1.05 * maxray * sqrt(1 - 
      #                                                    i^2), i, cex = 0.6)
      #     text(-1.05 * maxray * i, 1.05 * maxray * sqrt(1 - 
      #                                                     i^2), -i, cex = 0.6)
      #   }
      #   seq.sd <- seq.int(0, 2 * maxray, by = (maxray/10))[-1]
      #   for (i in seq.sd) {
      #     xcircle <- sd.r + (cos(discrete * pi/180) * 
      #                          i)
      #     ycircle <- sin(discrete * pi/180) * i
      #     for (j in 1:length(xcircle)) {
      #       if ((xcircle[j]^2 + ycircle[j]^2) < (maxray^2)) {
      #         points(xcircle[j], ycircle[j], col = "darkgreen", 
      #                pch = ".")
      #         if (j == 10) 
      #           text(xcircle[j], ycircle[j], signif(i, 
      #                                               2), cex = 0.5, col = "darkgreen")
      #       }
      #     }
      #   }
      #   seq.sd <- seq.int(0, maxray, length.out = 5)
      #   for (i in seq.sd) {
      #     xcircle <- (cos(discrete * pi/180) * i)
      #     ycircle <- sin(discrete * pi/180) * i
      #     if (i) 
      #       lines(xcircle, ycircle, lty = 3, col = "blue")
      #     text(min(xcircle), -0.03 * maxray, signif(i, 
      #                                               2), cex = 0.5, col = "blue")
      #     text(max(xcircle), -0.03 * maxray, signif(i, 
      #                                               2), cex = 0.5, col = "blue")
      #   }
      #   text(0, -0.08 * maxray, "Standard Deviation", 
      #        cex = 0.7, col = "blue")
      #   text(0, -0.12 * maxray, "Centered RMS Difference", 
      #        cex = 0.7, col = "darkgreen")
      #   points(sd.r, 0, pch = 22, bg = "darkgreen", cex = 1.1)
      #   text(0, 1.1 * maxray, "Correlation Coefficient", 
      #        cex = 0.7)
      # }
      # S <- (2 * (1 + R))/(sd.f + (1/sd.f))^2
    }
  }
  points(sd.f * R, sd.f * sin(acos(R)), pch = pch, col = col,
         cex = pcex)
  
  if(! is.null(text)) {
    text(sd.f * R, sd.f * sin(acos(R)),  #the line to add
         labels=text, cex = pcex, pos=3) #You can change the pos argument to your liking
  }
  
  invisible(oldpar)
}


#' Plot a heatmap to demonstrate the hypothetical time delay between any two genes 
#' 
#' This plot function uses the estimated turning_point for each gene in the CDS to estimate the time delay between any gene pairs. 
#' It then plot the time delay using a heatmap. 
#' 
#' @param cds CellDataSet for the experiment
#' @param gene_pairs_mat A matrix where each row is the gene pair and the first column is the hypothetical source or regulator while the second column represents the hypothetical target. The name in this matrix should match the name in the gene_short_name column of the cds object. 
#' @param n_row number of columns used to layout the faceted cluster panels
#' @param n_col number of columns used to layout the faceted cluster panels
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
# function to plot the time delay between each gene: 
plot_time_delay_heatmap <- function(cds, use_gene_short_name = TRUE, cluster_row = TRUE, cluster_cols = TRUE) {
  if(!("turning_point" %in% colnames(fData(cds)))) 
    stop('Please first run the xxx function before running this function')
  
  turning_point <- fData(cds)[, c('gene_short_name', 'turning_point')] 
  
  if(use_gene_short_name) {
    if(any(duplicated(turning_point$gene_short_name)))
      stop('gene_short_name column in fData has duplicated points!')
    
    row.names(turning_point) <- turning_point$gene_short_name
  }
  
  time_delay_mat <- matrix(nrow = nrow(turning_point), ncol = nrow(turning_point), 
                           dimnames = list(row.names(turning_point), row.names(turning_point)))
  for(gene_i in row.names(turning_point)) {
    for(gene_j in row.names(turning_point)) {
      time_delay_mat[gene_i, gene_j] <- abs(abs(turning_point[gene_i, 'turning_point']) - abs(turning_point[gene_j, 'turning_point']))
    }
  }
  
  if(any(!is.finite(time_delay_mat))) {
    message('There is NA values in turining points calculated, the time delay is set to 0 by default')
    
    time_delay_mat[!is.finite(time_delay_mat)] <- 0
  }
  
  ph_res <- pheatmap::pheatmap(time_delay_mat[, ], #ph$tree_row$order
                               useRaster = T,
                               cluster_cols=cluster_row, 
                               cluster_rows=cluster_row, 
                               show_rownames=T, 
                               show_colnames=T, 
                               #scale="row",
                               # clustering_distance_rows=row_dist, #row_dist
                               clustering_method = "ward.D2", #
                               # cutree_rows=num_clusters,
                               # cutree_cols = 2,
                               # annotation_row=annotation_row,
                               # annotation_col=annotation_col,
                               # annotation_colors=annotation_colors,
                               # gaps_col = col_gap_ind,
                               treeheight_row = 20, 
                               # breaks=bks,
                               fontsize = 6,
                               # color=hmcols, 
                               silent=TRUE)
  
  grid::grid.rect(gp=grid::gpar("fill", col=NA))
  grid::grid.draw(ph_res$gtable)
}

#' Plot the causal network inferred by Scribe 
#' 
#' This plot function uses state-of-art graph visualization tool (including Hive plot, hiearchy plot and arc diagram) to visualize the causal network. 
#' 
#' @param graph A igraph object describing inferred causal network 
#' @param type A character describing which type of the plot will be used for making the network plot  
#' @param layout A layout function which will only be used when we set type = 'igraph'. 
#' @return a plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung() 
#' gene_pairs_mat <- matrix(c('H19', 'Ccnd2', 'Ccnd2', 'Scnn1g'), ncol = 2)
#' plot_causal_network(lung, gene_pairs_mat)
#' }
#' @references Hive plot, hiearchy plot and arc diagram plot, etc. 
#' @export
#' 
# plot_time_delay_heatmap(lung, use_gene_short_name = F)
plot_causal_network <- function(graph, type = c('igraph', 'hiearchy', 'hive', 'arcdiagram'), layout = NULL, ...) {
  if(class(graph) != 'igraph' | length(E(graph)) == 0)
    stop('Please pass an igraph object for the graph object and make sure it has at least one edge')
  
  if(type == 'igraph') {
    if(!is.null(layout)) {
      if(!is.function(layout)) {
        stop("Please make sure layout is a supported FUNCTION from igraph or customized layout FUNCTION")
      } else {
        layout_coord = layout(graph)
      }
    }
    
    plot(graph, layout = layout_coord) #, vertex.size=2, vertex.label=NA, vertex.color = 'black'        
  }
  else if(type == 'hiearchy') {
    res <- level.plot(graph) # create the hiearchical plot 
    
    cus_layout <- res$layout
    master_regulator_id <- cus_layout[, 2] %in% min(unique(cus_layout[, 2])) #1:2 #
    direct_target_id <- cus_layout[, 2] %in% unique(cus_layout[, 2])[order(unique(cus_layout[, 2])) == 2]
    secondary_target_id <- cus_layout[, 2] %in% max(unique(cus_layout[, 2]))
    
    cus_layout[master_regulator_id, 1]
    cus_layout[master_regulator_id, 1] <- c(-10, 10)
    cus_layout[direct_target_id, 1] <- cus_layout[direct_target_id, 1] * 1.5 #order(cus_layout[direct_target_id, 1]) * 2  - length(cus_layout[direct_target_id, 1]) * 20
    cus_layout[secondary_target_id, 1] <- cus_layout[secondary_target_id, 1] * 1#order(cus_layout[secondary_target_id, 1]) * 10 - length(cus_layout[secondary_target_id, 1])  * 5
    
    v_size <- rep(res$vertex.size, nrow(cus_layout))
    v_size[master_regulator_id] <- 12
    v_size[direct_target_id] <- 4#v_size[direct_target_id] * 1.5
    v_size[secondary_target_id] <- 2.5
    
    cus_layout[master_regulator_id, 2] <- 0
    cus_layout[direct_target_id, 2] <- -1
    cus_layout[secondary_target_id, 2] <- -2
    
    res$vertex.label.cex <- rep(0.25, length(v_size))
    res$vertex.label.cex[1:2] <- 1
    res$layout <- cus_layout
    res$bg <- 'white'
    res$vertex.size <- v_size
    
    res$vertex.label <- V(res$g)$name
    res$lab.color <- 'black'
    
    res$vertex.color[1:2] <- c('#BCA0CC')
    res$vertex.color[direct_target_id] <- '#77CCD2'
    res$vertex.color[res$vertex.color == 'orange2'] <- '#7EB044'
    
    secondary_layer_degree <- degree(res$g, v = V(res$g)$name[cus_layout[, 2] == -1], mode = c("out"),
                                     loops = TRUE, normalized = FALSE)
    res$vertex.color[which(cus_layout[, 2] == -1)[secondary_layer_degree == 0]] <- '#F3756C'
    # secondary_layer <- degree(res$g, v = V(res$g)$name[cus_layout[, 2] == -1], mode = c("out"),
    #                           loops = TRUE, normalized = FALSE)
    
    res$vertex.frame.color <- res$vertex.color
    res$vertex.label <- c(res$vertex.label[1:2], rep('', length(res$vertex.label) - 2))
    
    netbiov::plot.netbiov(res)
  }
  else if(type == 'hive') {
    
    # use the hiveR package 
  }
  else if(type == 'arcdiagram') {
    edge <- get.edgelist(graph)
    
    arcplot(edge, ordering=sort(unique(c(edge[, 1], edge[, 2]))), 
            horizontal=TRUE,
            #labels=paste("node",1:10,sep="-"),
            #lwd.arcs=4*runif(10,.5,2), 
            col.arcs=hsv(runif(9,0.6,0.8),alpha=0.4),
            show.nodes=TRUE, pch.nodes=21, cex.nodes=runif(10,1,3), 
            col.nodes="gray80", bg.nodes="gray90", lwd.nodes=2)
    
    # # create the figure according the hiearchical plot 
    # arcplot(edge, ordering=unique(c(edge[, 1], edge[, 2])), horizontal=TRUE,
    #         #labels=paste("node",1:10,sep="-"),
    #         #lwd.arcs=4*runif(10,.5,2), 
    #         col.arcs=hsv(abs(cus_layout) / max(abs(cus_layout)),alpha=0.4),
    #         show.nodes=TRUE, pch.nodes=21, cex.nodes=runif(10,1,3), 
    #         col.nodes="gray80", bg.nodes="gray90", lwd.nodes=2)
  }
}

#' Plot the relationship between the predicted value (\hat{x}(t) | M_y) and the observed x(t). 
#' 
#' This plotting function uses the simplex projection from the shadow manifold to estimate the source gene's value and then 
#' plot against the observe value. If there is a true causality from x to y, the predicted point should match up with 
#' the observed value pretty well. It also supports to plot the predicted value and the weighted target values based on 
#' the nearest neighbors used in the simplex projection. This two different plots are controlled by the type argument, 
#' where "Prediction" and "Relationship" corresponds to first and second plot respectively. This function accepts a matrix 
#' where each row is the gene pair and the first column is the hypothetical source or regulator while the second column 
#' represents the hypothetical target. The name in this matrix should match the name in the gene_short_name column of the 
#' cds object.
#' 
#' @param cds CellDataSet for the experiment
#' @param gene_pairs_mat A matrix where each row is the gene pair and the first column is the hypothetical source or regulator while the second column represents the hypothetical target. The name in this matrix should match the name in the gene_short_name column of the cds object. 
#' @param tau the lag to use for time delay embedding (same as in the ccm function from the rEDM package)
#' @param E the embedding dimensions to use for time delay embedding
#' @param Type 
#' @param verbose whether or not to print the warning messages to the R console
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
plot_cross_map <-function(cds_subset, gene_pairs_mat, tau = 1, E = 2, Type = c("Prediction", "Relationship"), verbose = FALSE, return_all = FALSE) 
{
  gene_pairs_mat <- as.matrix(gene_pairs_mat)
  all_genes_in_pair <- as.vector(gene_pairs_mat)
  if(! all(all_genes_in_pair %in% fData(cds_subset)$gene_short_name)) {
    stop("cds_subset doesn't include all genes in gene_pairs_mat Make sure all genes are included in gene_short_name column of the cds")
  }
  
  cds_subset <- cds_subset[row.names(subset(fData(cds_subset), gene_short_name %in% all_genes_in_pair)), ]
  
  Time <- ncol(cds_subset)
  n_shadow <- Time - (E - 1)*tau # maximal number of points for shadow manifold 
  
  flat_res <- as.data.frame(matrix(0, nrow = n_shadow * 2 * nrow(gene_pairs_mat), ncol = 4))
  colnames(flat_res) <- c("Observations", "Predictions", "Weighted_source", "Pair")
  
  id <- 0
  for(gene_pairs_ind in 1:nrow(gene_pairs_mat))  {
    # gene_pairs_ind <- 1
    if(verbose)
      message("current gene pair is ", gene_pairs_mat[gene_pairs_ind, 1], " -> ",  gene_pairs_mat[gene_pairs_ind, 2])
    
    gene_pairs <- gene_pairs_mat[gene_pairs_ind, ]
    X <- as.numeric(exprs(cds_subset)[row.names(subset(fData(cds_subset), gene_short_name %in% gene_pairs[1])), ])# as.numeric(exprs(cds_subset)[gene_pairs[1], ])
    Y <- as.numeric(exprs(cds_subset)[row.names(subset(fData(cds_subset), gene_short_name %in% gene_pairs[2])), ])
    
    laggedX<-gen_lagged_vec(X=X,E=E,tau=tau)
    laggedY<-gen_lagged_vec(X=Y,E=E,tau=tau)
    dists_x_all <- as.matrix(dist(laggedX$x))
    dists_y_all <- as.matrix(dist(laggedY$x))
    max_shadow <- dim(dists_x_all)[1]
    # cor_x <- numeric(Time)
    # cor_y <- numeric(Time)
    
    # for(L in (5+(E-1)*tau):Time) # why 5 + (E - 1) * tau here? 
    # {   
    L = Time
    n_shadow <- L-(E-1)*tau
    pred_x <- numeric(n_shadow)
    pred_y <- numeric(n_shadow)
    actual_x <- numeric(n_shadow)
    actual_y <- numeric(n_shadow) 
    
    source_x <- numeric(n_shadow)
    source_y <- numeric(n_shadow)
    
    # cor_x_tmp <- numeric(reps)
    # cor_y_tmp <- numeric(reps)
    
    start <- 1
    # if(randomWindow)
    #     start <- sample(1:(max_shadow-n_shadow+1),1)
    window <- start:(start+n_shadow-1)
    
    dists_x <- dists_x_all[window,window]
    dists_y <- dists_y_all[window,window]
    
    x=laggedX$x[window,]
    y=laggedY$x[window,]
    nn_x <- get_nn(x,dists_x)
    nn_y <- get_nn(y,dists_y)
    
    
    for(i in 1:n_shadow)
    {
      ww <- calc_weights(x=x,dists=dists_x,index=i,nn=nn_x)
      pred_y[i] <- predict_y(y=y,wt=ww)
      actual_y[i] <- y[i,1]
      source_x[i] <- predict_y(y=x,wt=ww)
      
      ww <- calc_weights(x=y,dists=dists_y,index=i,nn=nn_y)
      pred_x[i] <- predict_y(y=x,wt=ww)
      actual_x[i] <- x[i,1] 
      source_y[i] <- predict_y(y=y,wt=ww)
    }      
    
    # cor_x[L] <- cor(pred_x,actual_x)
    # cor_y[L] <- cor(pred_y,actual_y)
    
    # cat(L,'\t',cor_x[L],'\t',cor_y[L],'\n')
    
    data <- data.frame(Observations = c(actual_x, actual_y), Predictions = c(pred_x, pred_y), 
                       Weighted_source = c(source_y, source_x), 
                       Pair = c(rep(paste0(gene_pairs[2], ' xmap ', gene_pairs[1]), n_shadow), 
                                rep(paste0(gene_pairs[1], ' xmap ', gene_pairs[2]), n_shadow)), stringsAsFactors = FALSE)
    
    flat_res[(id + 1):(id + 2 * n_shadow), ] <- data
    
    id <- id + 2 * n_shadow
  }
  
  flat_res$Pair <- factor(flat_res$Pair, levels = unique(flat_res$Pair))
  
  if(Type == 'Prediction')
    g <- qplot(Observations, Predictions, data = flat_res, geom = 'point', log = 'xy', size = I(1), size = 2) +
    xlab('') + ylab('') + geom_point(color = I("black"), size = I(0.85)) + geom_smooth(method = 'lm') + 
    monocle:::monocle_theme_opts() + facet_wrap(~Pair, ncol = 2, scales = 'free') + geom_abline(color = 'red')
  else if(Type == 'Relationship')
    g <- qplot(Predictions, Weighted_source, data = res$res, geom = 'point', log = 'xy', size = I(1), size = 2) +
    xlab('') + ylab('') + geom_point(color = I("black"), size = I(0.85)) + geom_smooth(method = 'loess') + 
    monocle:::monocle_theme_opts() + facet_wrap(~Pair, ncol = 2, scales = 'free') + geom_abline(color = 'red')
  
  if(return_all) {
    return(list(res = flat_res, g = g))
  } else {
    return(g)
  }
}

meshgrid <- function (x, y = x) {
  if (!is.numeric(x) || !is.numeric(y)) 
    stop("Arguments 'x' and 'y' must be numeric vectors.")
  x <- c(x)
  y <- c(y)
  n <- length(x)
  m <- length(y)
  X <- matrix(rep(x, each = m), nrow = m, ncol = n)
  Y <- matrix(rep(y, times = n), nrow = m, ncol = n)
  return(list(X = X, Y = Y))
}

# add the functionality to normalize either row or column for the data 
#' Plot the heatmap for the expected value y(t) given x(t - d) and y(t - 1). 
#' 
#' This plotting function tries to intuitively visualize the informatioin transfer from x(t - d) to y(t) given y(t)'s previous state y(t - 1).
#' Firstly, we divide the expression space for x(t - d) to y(t - 1) based on grid_num and then we estimate the k-nearest neighbor for each of the
#' grid. We then use a Gaussian kernel to estimate the expected value for y(t). It is then displayed in two dimension with x(t - d) and y(t - 1)
#' as two axis and the color represents the value of the expected of y(t) This function accepts a matrix where each row is the gene pair and the 
#' first column is the hypothetical source or regulator while the second column represents the hypothetical target. The name in this matrix should 
#' match the name in the gene_short_name column of the cds object.
#' 
#' @param cds_subset CellDataSet for the experiment
#' @param gene_pairs_mat A matrix where each row is the gene pair and the first column is the hypothetical source or regulator while the second column represents the hypothetical target. The name in this matrix should match the name in the gene_short_name column of the cds object. 
#' @param log A logic argument to determine whether or not we should log the data before we perform density estimation.
#' @param d The time delay between the source gene and the target gene. 
#' @param k The number of the nearest neighbors used in calculating the expectation value of y(t). 
#' @param grid_num The number of grids used in estimating density. 
#' @param n_row number of columns used to layout the faceted cluster panels
#' @param n_col number of columns used to layout the faceted cluster panels
#' @param scales 
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
# create the heatmap for the rdi value: 
plot_rdi_pairs_heatmap <- function(cds_subset, gene_pairs_mat, 
                                   log = FALSE, 
                                   d = 1, 
                                   k = 5, 
                                   grid_num = 25,
                                   n_row = NULL,
                                   n_col = 1,
                                   scales = "free",
                                   verbose = FALSE) {
  
  gene_pairs_mat <- as.matrix(gene_pairs_mat)
  all_genes_in_pair <- as.vector(gene_pairs_mat)
  if(! all(all_genes_in_pair %in% fData(cds_subset)$gene_short_name)) {
    stop("cds_subset doesn't include all genes in gene_pairs_mat Make sure all genes are included in gene_short_name column of the cds")
  }
  
  cds_subset <- cds_subset[row.names(subset(fData(cds_subset), gene_short_name %in% all_genes_in_pair)), ]
  
  flat_res <- as.data.frame(matrix(0, nrow = grid_num^2  * nrow(gene_pairs_mat), ncol = 4))
  colnames(flat_res) <- c("x", "z", "expected_y", "pair")
  # xy <- data.frame()
  
  id <- 0
  for(gene_pairs_ind in 1:nrow(gene_pairs_mat))  {
    # gene_pairs_ind <- 1
    if(verbose)
      message("current gene pair is ", gene_pairs_mat[gene_pairs_ind, 1], " -> ",  gene_pairs_mat[gene_pairs_ind, 2])
    
    gene_pairs <- gene_pairs_mat[gene_pairs_ind, ]
    f_ini_ind <- grid_num^2 * id #flat_res (normalized density results)
    
    gene_pair_name <- paste(gene_pairs[1], gene_pairs[2], sep = ' -> ')
    
    x <- matrix(exprs(cds_subset)[row.names(subset(fData(cds_subset), gene_short_name %in% gene_pairs[1])), ], ncol = 1)
    y_ori <- matrix(exprs(cds_subset)[row.names(subset(fData(cds_subset), gene_short_name %in% gene_pairs[2])), ], ncol = 1)
    
    if(log) {
      x <- log(x + 1)
      y_ori <- log(y_ori + 1)
    }
    
    if(d != 0) {
      x  <- matrix(x[1:(nrow(x) - d), ], ncol = 1)
      y <- matrix(y_ori[-(1:d), ], ncol = 1)
      z <- y_ori[d:(nrow(y_ori) - 1), ]
    }
    else {
      x  <- matrix(x, ncol = 1)
      y <- matrix(y_ori, ncol = 1)
      z <- y_ori
    }
    
    if(length(unique(x)) < grid_num | length(unique(y)) < grid_num) {
      # stop(paste0("Genes ", gene_pairs, "are required to express in at least ", grid_num, " cells"))
      # x <- x + rnorm(length(x), sd = 1e-10)
    }
    
    data <- data.frame(x = x, y = y, z = z)
    
    rng_x <- range(x)
    rng_z <- range(z)
    x_meshgrid <- seq(rng_x[1], rng_x[2], length.out = grid_num)
    z_meshgrid <- seq(rng_z[1], rng_z[2], length.out = grid_num)
    
    xz_meshgrid <- meshgrid(x_meshgrid, z_meshgrid)
    xz_query <- matrix(c(as.vector(xz_meshgrid$X), as.vector(xz_meshgrid$Y)), ncol = 2)
    
    knn_res <- RANN::nn2(data[, 1:2], query = xz_query, k = k + 1)
    
    dist_mat <- apply(knn_res$nn.idx, 1, function(x) {
      abs(data[x[-1], 2] - data[x[1], 2])
    })
    dist_mat <- knn_res$nn.dists
    
    for(index in 1:nrow(dist_mat)) {
      u <- exp(-dist_mat[index, -1]/ min(dist_mat[index, -1]))
      w <- u / sum(u)
      
      subset_dat <- data[knn_res$nn.idx[index, -1], 2]
      tmp <- sum(w * subset_dat)
      flat_res[f_ini_ind + index, ] <- c(xz_query[index, ], tmp, gene_pair_name)
    }
    
    # row-scale the data 
    vals <- as.numeric(flat_res[c(f_ini_ind + 1):c(f_ini_ind + nrow(dist_mat)), 3])
    max_val <- max(vals, na.rm = T) 
    if(!is.finite(max_val))
      max_val <- 1e10
    
    print(max_val)
    for(index in 1:nrow(dist_mat)) {
      flat_res[f_ini_ind + index, 3] <- as.numeric(flat_res[f_ini_ind + index, 3])  / max_val
    }
    
    id <- id + 1
  }
  
  flat_res[, 1:3] <- as.matrix(flat_res[, 1:3])
  flat_res[, 4] <- factor(flat_res[, 4], levels = unique(flat_res[, 4]))
  ggplot(aes(as.numeric(x), as.numeric(z)), data = flat_res) +  geom_raster(aes(fill = as.numeric(expected_y))) + 
    scale_fill_gradientn(expression(paste("Expection (", "y", phantom()[{
      paste("t")
    }], ")", "")), colours = terrain.colors(10)) + 
    # geom_rug(aes(as.numeric(x), as.numeric(y)), data = xy, col="darkred",alpha=.1) + 
    # geom_path(aes(as.numeric(x), as.numeric(y)), data  = ridge_curve, color = 'red') + 
    facet_wrap(~pair, scales = scales, nrow = n_row, ncol = n_col) + 
    xlab( expression(paste("Source (", "x", phantom()[{
      paste("t", phantom() - phantom(), d)
    }], ")", ""))) + 
    ylab(expression(paste("Target (", "y", phantom()[{
      paste("t", phantom() - phantom(), "1")
    }], ")", ""))) + monocle:::monocle_theme_opts()
}

#' Plot the combinatorial influence of two genes x, y to the target z. 
#' 
#' This plotting function tries to intuitively visualize the influence from genes x and y to the target z. 
#' Firstly, we divide the expression space for x and y based on grid_num and then we estimate the k-nearest neighbor for each of the
#' grid. We then use a Gaussian kernel to estimate the expected value for z. It is then displayed in two dimension with x and y
#' as two axis and the color represents the value of the expected of z. This function accepts a matrix where each row is the gene pair
#' and the target genes for this pair. The first column is the first hypothetical source or regulator, the second column represents 
#' the second hypothetical target while the third column represents the hypothetical target gene. The name in this matrix should match 
#' the name in the gene_short_name column of the cds object.
#'  
#' @param cds CellDataSet for the experiment
#' @param gene_pairs_mat A matrix where each row is the gene pair and the first column is the hypothetical source or regulator while the second column represents the hypothetical target. The name in this matrix should match the name in the gene_short_name column of the cds object. 
#' @param n_row number of columns used to layout the faceted cluster panels
#' @param n_col number of columns used to layout the faceted cluster panels
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
plot_comb_logic_heatmap <- function(cds_subset, gene_pairs_target_mat, 
                                    log = FALSE, 
                                    d = 1,
                                    grid_num = 25,
                                    n_row = NULL,
                                    n_col = 1,
                                    normalized = TRUE,
                                    scales = "free",
                                    verbose = FALSE) {
  
  gene_pairs_target_mat <- as.matrix(gene_pairs_target_mat)
  all_genes_in_pair <- as.vector(gene_pairs_target_mat)
  if(! all(all_genes_in_pair %in% fData(cds_subset)$gene_short_name)) {
    stop("cds_subset doesn't include all genes in gene_pairs_target_mat Make sure all genes are included in gene_short_name column of the cds")
  }
  
  cds_subset <- cds_subset[row.names(subset(fData(cds_subset), gene_short_name %in% all_genes_in_pair)), ]
  
  flat_res <- as.data.frame(matrix(0, nrow = grid_num^2  * nrow(gene_pairs_target_mat), ncol = 4))
  colnames(flat_res) <- c("x", "y", "expected_z", "pair")
  # xy <- data.frame()
  
  id <- 0
  for(gene_pairs_ind in 1:nrow(gene_pairs_target_mat))  {
    # gene_pairs_ind <- 1
    if(verbose)
      message("current gene pair is ",  paste0(gene_pairs[1], " : ", gene_pairs[2], " -> ", gene_pairs[3]))
    
    gene_pairs <- gene_pairs_target_mat[gene_pairs_ind, ]
    f_ini_ind <- grid_num^2 * id #flat_res (normalized density results)
    
    gene_pair_name <- paste0(gene_pairs[1], " : ", gene_pairs[2], " -> ", gene_pairs[3])
    
    x <- matrix(exprs(cds_subset)[row.names(subset(fData(cds_subset), gene_short_name %in% gene_pairs[1])), ], ncol = 1)
    y <- matrix(exprs(cds_subset)[row.names(subset(fData(cds_subset), gene_short_name %in% gene_pairs[2])), ], ncol = 1)
    z <- matrix(exprs(cds_subset)[row.names(subset(fData(cds_subset), gene_short_name %in% gene_pairs[3])), ], ncol = 1)
    
    if(log) {
      x <- log(x + 1)
      y <- log(y + 1)
      z <- log(z + 1)
    }
    
    data <- data.frame(x = x, y = y, z = z)
    
    rng_x <- range(x)
    rng_y <- range(y)
    x_meshgrid <- seq(rng_x[1], rng_x[2], length.out = grid_num)
    y_meshgrid <- seq(rng_y[1], rng_y[2], length.out = grid_num)
    
    xy_meshgrid <- meshgrid(x_meshgrid, y_meshgrid)
    xy_query <- matrix(c(as.vector(xy_meshgrid$X), as.vector(xy_meshgrid$Y)), ncol = 2)
    
    knn_res <- RANN::nn2(data[, 1:2], query = xy_query, k = k + 1)
    
    dist_mat <- apply(knn_res$nn.idx, 1, function(x) {
      abs(data[x[-1], 3] - data[x[1], 3])
    })
    dist_mat <- knn_res$nn.dists
    
    for(index in 1:nrow(dist_mat)) {
      u <- exp(-dist_mat[index, -1]/ min(dist_mat[index, -1]))
      w <- u / sum(u)
      
      subset_dat <- data[knn_res$nn.idx[index, -1], 3]
      tmp <- sum(w * subset_dat)
      flat_res[f_ini_ind + index, ] <- c(xy_query[index, ], tmp, gene_pair_name)
    }
    
    # row-scale the data 
    if(normalized) {
      vals <- as.numeric(flat_res[c(f_ini_ind + 1):c(f_ini_ind + nrow(dist_mat)), 3])
      max_val <- max(vals, na.rm = T) 
      if(!is.finite(max_val))
        max_val <- 1e10
      
      print(max_val)
      for(index in 1:nrow(dist_mat)) {
        flat_res[f_ini_ind + index, 3] <- as.numeric(flat_res[f_ini_ind + index, 3])  / max_val
      }
    }
    
    id <- id + 1
  }
  
  flat_res[, 1:3] <- as.matrix(flat_res[, 1:3])
  flat_res[, 4] <- factor(flat_res[, 4], levels = unique(flat_res[, 4]))
  ggplot(aes(as.numeric(x), as.numeric(y)), data = flat_res) +  geom_raster(aes(fill = as.numeric(expected_z))) + 
    scale_fill_gradientn(expression(paste("Expection (", "z", phantom()[{
      paste("t")
    }], ")", "")), colours = terrain.colors(10)) + 
    # geom_rug(aes(as.numeric(x), as.numeric(y)), data = xy, col="darkred",alpha=.1) + 
    # geom_path(aes(as.numeric(x), as.numeric(y)), data  = ridge_curve, color = 'red') + 
    facet_wrap(~pair, scales = scales, nrow = n_row, ncol = n_col) + 
    xlab( expression(paste("Source (", "x", ")", ""))) + 
    ylab(expression(paste("Source (", "y", ")", ""))) + monocle:::monocle_theme_opts()
}
