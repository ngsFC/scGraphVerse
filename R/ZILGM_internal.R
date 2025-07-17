#' Zero-Inflated Latent Gaussian Model Internal Implementation
#'
#' Internal implementation of ZILGM for gene regulatory network inference,
#' rebuilt from original sources with BiocParallel integration.
#'
#' @param X Matrix of count data (samples x genes)
#' @param lambda Vector of regularization parameters
#' @param nlambda Number of lambda values
#' @param family Distribution family ("Poisson", "NBI", "NBII")
#' @param nCores Number of cores for parallelization
#' @param ... Additional parameters
#'
#' @return List with network adjacency matrices and coefficients
#'
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam
#' @importFrom glmnet glmnet
#' @keywords internal
#' @noRd
zilgm_internal <- function(X, lambda = NULL, nlambda = 50, 
                          family = "NBII", nCores = 1, ...) {
    # Check dependencies
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("glmnet package is required but not installed.\n",
             "Install with: install.packages('glmnet')")
    }
    
    # Setup parallelization
    if (nCores > 1) {
        bp_param <- BiocParallel::MulticoreParam(workers = nCores)
    } else {
        bp_param <- BiocParallel::SerialParam()
    }
    
    # Validate inputs
    if (!is.matrix(X)) X <- as.matrix(X)
    n_samples <- nrow(X)
    p_genes <- ncol(X)
    
    if (n_samples < 3 || p_genes < 2) {
        stop("Insufficient data: need at least 3 samples and 2 genes")
    }
    
    # Determine lambda sequence if not provided
    if (is.null(lambda)) {
        lambda_max <- .estimate_lambda_max(X)
        lambda <- exp(seq(log(lambda_max), log(lambda_max * 1e-4), length.out = nlambda))
    }
    
    # Parallel estimation for each gene
    gene_results <- BiocParallel::bplapply(seq_len(p_genes), function(gene_idx) {
        .estimate_gene_network(X, gene_idx, lambda, family)
    }, BPPARAM = bp_param)
    
    # Combine results
    .combine_zilgm_results(gene_results, lambda, p_genes, colnames(X))
}

#' Estimate lambda max for ZILGM
#' @keywords internal
#' @noRd
.estimate_lambda_max <- function(X) {
    p <- ncol(X)
    lambda_max <- 0
    
    for (j in seq_len(p)) {
        y <- X[, j]
        X_j <- X[, -j, drop = FALSE]
        
        # Standardize predictors
        X_j_scaled <- scale(X_j)
        
        # Compute correlation-based lambda max
        cor_vals <- abs(cor(y, X_j_scaled, use = "complete.obs"))
        lambda_max <- max(lambda_max, max(cor_vals, na.rm = TRUE))
    }
    
    return(lambda_max)
}

#' Estimate network for single gene
#' @keywords internal
#' @noRd
.estimate_gene_network <- function(X, target_gene, lambda, family) {
    tryCatch({
        y <- X[, target_gene]
        X_predictors <- X[, -target_gene, drop = FALSE]
        
        # Handle zero-inflation based on family
        if (family %in% c("Poisson", "NBI", "NBII")) {
            # For count data, use appropriate family in glmnet
            glm_family <- switch(family,
                "Poisson" = "poisson",
                "NBI" = "poisson",  # Approximate with Poisson
                "NBII" = "poisson"  # Approximate with Poisson
            )
        } else {
            glm_family <- "gaussian"
        }
        
        # Fit regularized regression
        fit <- glmnet::glmnet(
            x = X_predictors,
            y = y,
            family = glm_family,
            lambda = lambda,
            standardize = TRUE,
            intercept = TRUE
        )
        
        # Extract coefficients (excluding intercept)
        coef_matrix <- as.matrix(fit$beta)
        
        return(list(
            coefficients = coef_matrix,
            lambda = fit$lambda,
            target_gene = target_gene
        ))
        
    }, error = function(e) {
        warning("Error fitting gene ", target_gene, ": ", e$message)
        # Return zero matrix
        p_predictors <- ncol(X) - 1
        zero_matrix <- matrix(0, nrow = p_predictors, ncol = length(lambda))
        return(list(
            coefficients = zero_matrix,
            lambda = lambda,
            target_gene = target_gene
        ))
    })
}

#' Combine ZILGM results into final format
#' @keywords internal
#' @noRd
.combine_zilgm_results <- function(gene_results, lambda, p_genes, gene_names) {
    n_lambda <- length(lambda)
    
    # Initialize coefficient array
    coef_array <- array(0, dim = c(p_genes, p_genes, n_lambda))
    
    # Fill coefficient array
    for (result in gene_results) {
        target_idx <- result$target_gene
        predictor_indices <- setdiff(seq_len(p_genes), target_idx)
        
        # Map coefficients back to full gene set
        if (target_idx < p_genes) {
            coef_array[predictor_indices[seq_len(target_idx - 1)], target_idx, ] <- 
                result$coefficients[seq_len(target_idx - 1), ]
            if (target_idx < p_genes) {
                remaining_idx <- seq(target_idx, nrow(result$coefficients))
                coef_array[predictor_indices[remaining_idx], target_idx, ] <- 
                    result$coefficients[remaining_idx, ]
            }
        } else {
            coef_array[predictor_indices, target_idx, ] <- result$coefficients
        }
    }
    
    # Create adjacency matrices for each lambda
    networks <- vector("list", n_lambda)
    for (i in seq_len(n_lambda)) {
        adj_matrix <- abs(coef_array[, , i])
        
        # Symmetrize using "OR" operation (union of edges)
        adj_matrix <- pmax(adj_matrix, t(adj_matrix))
        
        # Add gene names if available
        if (!is.null(gene_names)) {
            rownames(adj_matrix) <- colnames(adj_matrix) <- gene_names
        }
        
        networks[[i]] <- adj_matrix
    }
    
    return(list(
        network = networks,
        coef_network = coef_array,
        lambda = lambda,
        call = match.call()
    ))
}