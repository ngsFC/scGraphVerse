#' Joint Random Forest Internal Implementation
#'
#' Internal implementation of Joint Random Forest for gene regulatory network
#' inference, rebuilt from original sources with BiocParallel integration.
#'
#' @param X List of expression matrices (genes x cells) for each condition
#' @param ntree Number of trees in the random forest
#' @param mtry Number of variables randomly sampled at each split
#' @param genes.name Character vector of gene names
#' @param nCores Number of cores for parallelization
#'
#' @return Data frame with gene pairs and importance scores
#'
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam
#' @importFrom randomForest randomForest
#' @keywords internal
#' @noRd
JRF_internal <- function(X, ntree = 500, mtry = NULL, genes.name, nCores = 1) {
    # Check dependencies
    if (!requireNamespace("randomForest", quietly = TRUE)) {
        stop("randomForest package is required but not installed.\n",
             "Install with: install.packages('randomForest')")
    }
    
    # Setup parallelization
    if (nCores > 1) {
        bp_param <- BiocParallel::MulticoreParam(workers = nCores)
    } else {
        bp_param <- BiocParallel::SerialParam()
    }
    
    # Validate inputs
    nclasses <- length(X)
    if (nclasses < 1) stop("At least one expression matrix required")
    
    p <- nrow(X[[1]])
    if (is.null(mtry)) mtry <- max(1, floor(sqrt(p - 1)))
    
    # Check gene names
    if (length(genes.name) != p) {
        stop("Length of genes.name must match number of genes in expression matrices")
    }
    
    # Parallel processing of target genes
    target_results <- BiocParallel::bplapply(seq_len(p), function(target_idx) {
        .process_target_gene(X, target_idx, mtry, ntree, nclasses, p)
    }, BPPARAM = bp_param)
    
    # Combine results into importance matrix
    imp <- array(0, dim = c(p, p, nclasses))
    for (j in seq_len(p)) {
        if (!is.null(target_results[[j]])) {
            imp[-j, j, ] <- target_results[[j]]
        }
    }
    
    # Create final result matrix
    .create_jrf_result(imp, genes.name, p, nclasses)
}

#' Process single target gene for JRF
#' @keywords internal
#' @noRd
.process_target_gene <- function(X, target_idx, mtry, ntree, nclasses, p) {
    tryCatch({
        # Prepare data for current target gene
        sampsize <- vapply(X, ncol, integer(1))
        max_samples <- max(sampsize)
        
        # Create predictor matrix
        predictor_data <- matrix(0, nrow = max_samples, ncol = (p - 1) * nclasses)
        response_data <- numeric(max_samples)
        sample_weights <- numeric(max_samples)
        
        col_idx <- 1
        for (class_idx in seq_len(nclasses)) {
            n_samples <- sampsize[class_idx]
            if (n_samples > 0) {
                # Extract predictors (all genes except target)
                predictors <- t(X[[class_idx]][-target_idx, , drop = FALSE])
                predictor_data[seq_len(n_samples), 
                             seq(col_idx, col_idx + p - 2)] <- predictors
                
                # Extract response (target gene)
                response_data[seq_len(n_samples)] <- X[[class_idx]][target_idx, ]
                
                # Set sample weights for this class
                sample_weights[seq_len(n_samples)] <- 1.0
            }
            col_idx <- col_idx + (p - 1)
        }
        
        # Remove rows with all zeros
        valid_samples <- which(rowSums(abs(predictor_data)) > 0)
        if (length(valid_samples) < 10) return(NULL)
        
        predictor_data <- predictor_data[valid_samples, , drop = FALSE]
        response_data <- response_data[valid_samples]
        
        # Fit random forest
        rf_model <- randomForest::randomForest(
            x = predictor_data,
            y = response_data,
            ntree = ntree,
            mtry = min(mtry, ncol(predictor_data)),
            importance = TRUE,
            nodesize = 5
        )
        
        # Extract importance scores for each class
        importance_scores <- rf_model$importance[, "IncNodePurity"]
        
        # Reshape importance by class
        class_importance <- matrix(0, nrow = p - 1, ncol = nclasses)
        for (class_idx in seq_len(nclasses)) {
            start_col <- (class_idx - 1) * (p - 1) + 1
            end_col <- class_idx * (p - 1)
            if (end_col <= length(importance_scores)) {
                class_importance[, class_idx] <- importance_scores[start_col:end_col]
            }
        }
        
        return(class_importance)
        
    }, error = function(e) {
        warning("Error processing target gene ", target_idx, ": ", e$message)
        return(NULL)
    })
}

#' Create JRF result data frame
#' @keywords internal
#' @noRd
.create_jrf_result <- function(imp, genes.name, p, nclasses) {
    # Create gene pair combinations
    gene_pairs <- combn(seq_len(p), 2, simplify = FALSE)
    n_pairs <- length(gene_pairs)
    
    result_df <- data.frame(
        gene1 = character(n_pairs),
        gene2 = character(n_pairs),
        stringsAsFactors = FALSE
    )
    
    # Fill gene names
    for (i in seq_len(n_pairs)) {
        pair <- gene_pairs[[i]]
        result_df$gene1[i] <- genes.name[pair[1]]
        result_df$gene2[i] <- genes.name[pair[2]]
    }
    
    # Add importance scores for each class
    for (class_idx in seq_len(nclasses)) {
        importance_col <- numeric(n_pairs)
        
        for (i in seq_len(n_pairs)) {
            pair <- gene_pairs[[i]]
            gene1_idx <- pair[1]
            gene2_idx <- pair[2]
            
            # Average bidirectional importance
            imp1_to_2 <- imp[gene1_idx, gene2_idx, class_idx]
            imp2_to_1 <- imp[gene2_idx, gene1_idx, class_idx]
            importance_col[i] <- (imp1_to_2 + imp2_to_1) / 2
        }
        
        result_df[[paste0("importance_class_", class_idx)]] <- importance_col
    }
    
    return(result_df)
}