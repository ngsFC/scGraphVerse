#' Joint Random Forest Internal Implementation
#'
#' Internal implementation matching original CRAN JRF package for gene regulatory
#' network inference, with BiocParallel integration for parallelization.
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
        stop(
            "randomForest package is required but not installed.",
            "Install with: install.packages('randomForest')"
        )
    }

    # Setup parallelization
    if (nCores > 1) {
        bp_param <- BiocParallel::MulticoreParam(workers = nCores)
    } else {
        bp_param <- BiocParallel::SerialParam()
    }

    # Validate inputs - following original JRF structure
    nclasses <- length(X)
    if (nclasses < 1) stop("At least one expression matrix required")

    # Get sample sizes for each class
    sampsize <- vapply(X, ncol, integer(1))
    tot <- max(sampsize)
    p <- nrow(X[[1]])
    
    if (is.null(mtry)) {
        mtry <- if (p > 3) max(floor(p/3), 1) else floor(sqrt(p))
    }

    # Check gene names
    if (length(genes.name) != p) {
        stop(
            "Length of genes.name must match number of genes in ",
            "expression matrices"
        )
    }

    # Initialize importance array
    imp <- array(0, dim = c(p, p, nclasses))

    # Process each target gene (following original algorithm)
    target_results <- BiocParallel::bplapply(seq_len(p), function(j) {
        .process_target_gene_original(X, j, mtry, ntree, nclasses, p, sampsize, tot)
    }, BPPARAM = bp_param)

    # Fill importance matrix
    for (j in seq_len(p)) {
        if (!is.null(target_results[[j]])) {
            imp[-j, j, ] <- target_results[[j]]
        }
    }

    # Create final result following original format
    .create_jrf_result_original(imp, genes.name, p, nclasses)
}

#' Process single target gene following original JRF algorithm
#' @keywords internal
#' @noRd
.process_target_gene_original <- function(X, j, mtry, ntree, nclasses, p, sampsize, tot) {
    tryCatch({
        # Following original JRF algorithm exactly:
        # Create covariate matrix and response following line 494-500 in original
        covar <- matrix(0, (p - 1) * nclasses, tot)
        y <- matrix(0, nclasses, tot)
        
        for (c in seq_len(nclasses)) {
            # Extract target gene expression for class c
            y[c, seq_len(sampsize[c])] <- as.matrix(X[[c]][j, ])
            
            # Extract predictor genes (all except target j) for class c
            covar[
                seq((c - 1) * (p - 1) + 1, c * (p - 1)),
                seq_len(sampsize[c])
            ] <- X[[c]][-j, ]
        }
        
        # Call modified randomForest following original line 502
        rf_model <- randomForest::randomForest(
            x = t(covar),
            y = as.vector(t(y)),
            ntree = ntree,
            mtry = mtry,
            importance = TRUE,
            nodesize = 5
        )
        
        # Extract importance scores using custom importance function logic
        # Following original line 504: importance(jrf.out, scale=FALSE)
        imp_scores <- rf_model$importance[, "IncNodePurity"]
        
        # Reshape importance scores by class (following line 504)
        class_importance <- matrix(0, nrow = p - 1, ncol = nclasses)
        for (s in seq_len(nclasses)) {
            start_idx <- (p - 1) * (s - 1) + 1
            end_idx <- (p - 1) * (s - 1) + (p - 1)
            class_importance[, s] <- imp_scores[start_idx:end_idx]
        }
        
        return(class_importance)
    }, error = function(e) {
        return(NULL)
    })
}

#' Create JRF result data frame following original format
#' @keywords internal
#' @noRd
.create_jrf_result_original <- function(imp, genes.name, p, nclasses) {
    # Following original algorithm lines 508-516:
    # Derive importance score for each interaction and create output format
    
    # Create lower triangular index vectors following lines 485-488
    vec1 <- matrix(rep(genes.name, p), p, p)
    vec2 <- t(vec1)
    vec1 <- vec1[lower.tri(vec1, diag = FALSE)]
    vec2 <- vec2[lower.tri(vec2, diag = FALSE)]
    
    imp.final <- matrix(0, p * (p - 1) / 2, nclasses)
    
    # Derive importance score for each interaction (lines 509-512)
    for (s in seq_len(nclasses)) {
        imp.s <- imp[, , s]
        t.imp <- t(imp.s)
        # Average bidirectional importance (gene1->gene2 and gene2->gene1)
        imp.final[, s] <- (imp.s[lower.tri(imp.s, diag = FALSE)] + 
                          t.imp[lower.tri(t.imp, diag = FALSE)]) / 2
    }
    
    # Create output following original format (lines 514-516)
    out <- cbind(
        as.character(vec1),
        as.character(vec2),
        as.data.frame(imp.final),
        stringsAsFactors = FALSE
    )
    colnames(out) <- c(
        paste0('gene', 1:2),
        paste0('importance', seq_len(nclasses))
    )
    
    return(out)
}
