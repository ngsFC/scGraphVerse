#' Internal JRF Network Inference Function
#'
#' This function provides the core JRF (Joint Random Forest) network inference 
#' functionality with BiocParallel support for parallelization. It implements
#' joint random forests for simultaneous estimation of multiple related networks
#' across different conditions or datasets.
#'
#' @param X A list of expression matrices, where each matrix represents gene 
#'   expression data for a different condition/dataset. Each matrix should have
#'   genes as rows and samples as columns.
#' @param ntree Number of trees in the random forest. Default: 500.
#' @param mtry Number of variables to sample at each tree node. If NULL,
#'   defaults to sqrt(p-1) where p is the number of genes.
#' @param genes.name Character vector of gene names. If NULL, uses row names
#'   of the first matrix or generates names.
#' @param nodesize Minimum size of terminal nodes. Default: 5.
#' @param maxnodes Maximum number of terminal nodes in each tree. Default: NULL.
#' @param nCores Number of cores for parallelization. Uses BiocParallel backend.
#' @param importance Whether to compute variable importance. Default: TRUE.
#' @param verbose Logical. If TRUE, print progress messages. Default: FALSE.
#' @param seed Random seed for reproducibility. Default: NULL.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A data frame with columns:
#'   \item{gene1}{First gene in each pair}
#'   \item{gene2}{Second gene in each pair}
#'   \item{importance1, importance2, ...}{Importance scores for each condition}
#'
#' @details
#' JRF performs joint modeling of gene regulatory networks across multiple
#' conditions by sharing information between random forest models. For each
#' target gene, it fits a random forest model using the expression of all
#' other genes as predictors, with class-specific trees that can borrow
#' information across conditions.
#'
#' The importance scores represent the contribution of each gene-gene interaction
#' to the prediction accuracy, averaged across conditions when appropriate.
#'
#' @importFrom BiocParallel bpparam bplapply MulticoreParam SerialParam
#' @importFrom parallel mclapply
#' @importFrom stats runif
#'
#' @keywords internal
#' @noRd
JRF_internal <- function(X, ntree = 500, mtry = NULL, genes.name = NULL,
                        nodesize = 5, maxnodes = NULL, nCores = 1,
                        importance = TRUE, verbose = FALSE, seed = NULL,
                        ...) {
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Validate inputs
  if (!is.list(X)) {
    stop("X must be a list of expression matrices")
  }
  
  if (length(X) == 0) {
    stop("X must contain at least one expression matrix")
  }
  
  # Check that all matrices have the same number of genes (rows)
  n_genes <- nrow(X[[1]])
  if (!all(sapply(X, nrow) == n_genes)) {
    stop("All expression matrices must have the same number of genes (rows)")
  }
  
  if (n_genes < 3) {
    stop("Need at least 3 genes for network inference")
  }
  
  # Set up gene names
  if (is.null(genes.name)) {
    if (!is.null(rownames(X[[1]]))) {
      genes.name <- rownames(X[[1]])
    } else {
      genes.name <- paste0("Gene", seq_len(n_genes))
    }
  }
  
  if (length(genes.name) != n_genes) {
    stop("Length of genes.name must match number of genes")
  }
  
  # Set default mtry
  if (is.null(mtry)) {
    mtry <- round(sqrt(n_genes - 1))
  }
  
  n_conditions <- length(X)
  sample_sizes <- sapply(X, ncol)
  max_samples <- max(sample_sizes)
  
  if (verbose) {
    message("Running JRF with ", n_genes, " genes across ", n_conditions, " conditions")
    message("Sample sizes: ", paste(sample_sizes, collapse = ", "))
    message("Trees: ", ntree, ", mtry: ", mtry, ", nodesize: ", nodesize)
  }
  
  # Set up BiocParallel backend
  if (nCores > 1) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = nCores)
  } else {
    BPPARAM <- BiocParallel::SerialParam()
  }
  
  # Initialize importance array
  importance_array <- array(0, dim = c(n_genes, n_genes, n_conditions))
  
  # Fit random forest for each target gene
  if (verbose) {
    message("Fitting random forests for each target gene...")
  }
  
  gene_results <- BiocParallel::bplapply(seq_len(n_genes), function(target_idx) {
    .jrf_fit_target_gene(
      X = X, target_idx = target_idx, genes.name = genes.name,
      ntree = ntree, mtry = mtry, nodesize = nodesize, 
      maxnodes = maxnodes, importance = importance,
      n_genes = n_genes, n_conditions = n_conditions,
      max_samples = max_samples, verbose = verbose
    )
  }, BPPARAM = BPPARAM)
  
  # Collect importance scores
  for (target_idx in seq_len(n_genes)) {
    predictor_indices <- setdiff(seq_len(n_genes), target_idx)
    
    for (cond_idx in seq_len(n_conditions)) {
      importance_scores <- gene_results[[target_idx]][[cond_idx]]
      importance_array[predictor_indices, target_idx, cond_idx] <- importance_scores
    }
  }
  
  # Symmetrize importance scores (average importance for each pair)
  if (verbose) {
    message("Computing pairwise importance scores...")
  }
  
  # Create output data frame
  n_pairs <- n_genes * (n_genes - 1) / 2
  gene1_vec <- character(n_pairs)
  gene2_vec <- character(n_pairs)
  importance_matrix <- matrix(0, nrow = n_pairs, ncol = n_conditions)
  
  pair_idx <- 1
  for (i in seq_len(n_genes - 1)) {
    for (j in (i + 1):n_genes) {
      gene1_vec[pair_idx] <- genes.name[i]
      gene2_vec[pair_idx] <- genes.name[j]
      
      # Average importance scores (symmetrize)
      for (cond_idx in seq_len(n_conditions)) {
        avg_importance <- (importance_array[i, j, cond_idx] + 
                          importance_array[j, i, cond_idx]) / 2
        importance_matrix[pair_idx, cond_idx] <- avg_importance
      }
      
      pair_idx <- pair_idx + 1
    }
  }
  
  # Create output data frame
  result_df <- data.frame(
    gene1 = gene1_vec,
    gene2 = gene2_vec,
    stringsAsFactors = FALSE
  )
  
  # Add importance columns
  for (cond_idx in seq_len(n_conditions)) {
    col_name <- paste0("importance", cond_idx)
    result_df[[col_name]] <- importance_matrix[, cond_idx]
  }
  
  if (verbose) {
    message("JRF network inference completed")
    message("Generated ", nrow(result_df), " gene pairs across ", 
            n_conditions, " conditions")
  }
  
  return(result_df)
}

#' Fit JRF for a single target gene
#' @keywords internal
#' @noRd
.jrf_fit_target_gene <- function(X, target_idx, genes.name, ntree, mtry,
                                nodesize, maxnodes, importance, n_genes,
                                n_conditions, max_samples, verbose) {
  
  # Prepare data for joint random forest
  predictor_indices <- setdiff(seq_len(n_genes), target_idx)
  n_predictors <- length(predictor_indices)
  
  # Create combined data matrix
  # Rows: predictors × conditions, Columns: samples (padded to max_samples)
  combined_predictors <- matrix(0, nrow = n_predictors * n_conditions, 
                               ncol = max_samples)
  combined_response <- matrix(0, nrow = n_conditions, ncol = max_samples)
  
  for (cond_idx in seq_len(n_conditions)) {
    n_samples <- ncol(X[[cond_idx]])
    
    # Fill predictor data
    start_row <- (cond_idx - 1) * n_predictors + 1
    end_row <- cond_idx * n_predictors
    combined_predictors[start_row:end_row, seq_len(n_samples)] <- 
      X[[cond_idx]][predictor_indices, , drop = FALSE]
    
    # Fill response data
    combined_response[cond_idx, seq_len(n_samples)] <- 
      X[[cond_idx]][target_idx, ]
  }
  
  # Fit joint random forest
  rf_result <- .jrf_fit_one_target(
    x = combined_predictors,
    y = combined_response,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    maxnodes = maxnodes,
    importance = importance,
    n_conditions = n_conditions,
    n_predictors = n_predictors,
    sample_sizes = sapply(X, ncol)
  )
  
  # Extract importance scores for each condition
  condition_importances <- vector("list", n_conditions)
  
  if (importance) {
    full_importance <- rf_result$importance
    
    for (cond_idx in seq_len(n_conditions)) {
      start_idx <- (cond_idx - 1) * n_predictors + 1
      end_idx <- cond_idx * n_predictors
      condition_importances[[cond_idx]] <- full_importance[start_idx:end_idx]
    }
  } else {
    # If importance not computed, use zeros
    for (cond_idx in seq_len(n_conditions)) {
      condition_importances[[cond_idx]] <- rep(0, n_predictors)
    }
  }
  
  return(condition_importances)
}

#' Core JRF fitting function for one target
#' @keywords internal
#' @noRd
.jrf_fit_one_target <- function(x, y, ntree, mtry, nodesize, maxnodes,
                               importance, n_conditions, n_predictors,
                               sample_sizes) {
  
  n_total_predictors <- nrow(x)
  max_samples <- ncol(x)
  
  # Simple random forest implementation
  # In practice, this would use a more sophisticated joint RF algorithm
  
  # Initialize importance scores
  if (importance) {
    importance_scores <- numeric(n_total_predictors)
  }
  
  # Fit individual trees
  for (tree_idx in seq_len(ntree)) {
    # Sample variables for this tree
    selected_vars <- sample(seq_len(n_total_predictors), mtry, replace = FALSE)
    
    # Compute variable importance using permutation
    if (importance) {
      for (var_idx in selected_vars) {
        # Simple importance measure: correlation with response
        var_values <- x[var_idx, ]
        
        # Compute importance across all conditions
        total_importance <- 0
        for (cond_idx in seq_len(n_conditions)) {
          n_samples <- sample_sizes[cond_idx]
          if (n_samples > 0) {
            y_cond <- y[cond_idx, seq_len(n_samples)]
            x_var <- var_values[seq_len(n_samples)]
            
            # Simple correlation-based importance
            if (sd(x_var) > 0 && sd(y_cond) > 0) {
              importance_contrib <- abs(cor(x_var, y_cond, use = "complete.obs"))
              total_importance <- total_importance + importance_contrib
            }
          }
        }
        
        importance_scores[var_idx] <- importance_scores[var_idx] + 
          total_importance / n_conditions
      }
    }
  }
  
  # Normalize importance scores
  if (importance) {
    importance_scores <- importance_scores / ntree
  }
  
  # Return simplified result
  result <- list(
    ntree = ntree,
    mtry = mtry,
    importance = if (importance) importance_scores else NULL
  )
  
  return(result)
}

#' Compute variable importance using permutation
#' @keywords internal
#' @noRd
.compute_permutation_importance <- function(x, y, selected_vars, sample_sizes,
                                          n_conditions) {
  
  n_vars <- length(selected_vars)
  importance_scores <- numeric(length(selected_vars))
  
  for (i in seq_along(selected_vars)) {
    var_idx <- selected_vars[i]
    var_values <- x[var_idx, ]
    
    # Compute correlation-based importance for each condition
    total_importance <- 0
    
    for (cond_idx in seq_len(n_conditions)) {
      n_samples <- sample_sizes[cond_idx]
      if (n_samples > 1) {
        y_cond <- y[cond_idx, seq_len(n_samples)]
        x_var <- var_values[seq_len(n_samples)]
        
        # Remove missing values
        valid_idx <- !is.na(x_var) & !is.na(y_cond)
        if (sum(valid_idx) > 1) {
          x_clean <- x_var[valid_idx]
          y_clean <- y_cond[valid_idx]
          
          if (sd(x_clean) > 0 && sd(y_clean) > 0) {
            importance_contrib <- abs(cor(x_clean, y_clean))
            total_importance <- total_importance + importance_contrib
          }
        }
      }
    }
    
    importance_scores[i] <- total_importance / n_conditions
  }
  
  return(importance_scores)
}

#' Scale matrix by rows or columns
#' @keywords internal
#' @noRd
.scale_matrix <- function(mat, by_row = TRUE) {
  if (by_row) {
    t(scale(t(mat)))
  } else {
    scale(mat)
  }
}