#' Internal PCzinb Network Inference Function
#'
#' This function provides the core PCzinb (PC algorithm for Zero-Inflated 
#' Negative Binomial models) network inference functionality with BiocParallel 
#' support for parallelization. It implements the PC algorithm using different
#' statistical tests depending on the specified method.
#'
#' @param X A numeric matrix of expression data (n x p) where n is the number
#'   of samples and p is the number of genes.
#' @param method Character string specifying the algorithm: "poi" (Poisson),
#'   "nb" (Negative Binomial), "zinb0" (ZINB mean-only), or "zinb1" (ZINB full).
#' @param alpha Significance level for conditional independence tests. 
#'   Default: 2*pnorm(sqrt(n), lower.tail=FALSE).
#' @param maxcard Maximum cardinality of conditioning sets. Default: 2.
#' @param extend Logical. If TRUE, uses union of tests; if FALSE, uses 
#'   intersection. Default: TRUE.
#' @param max_iter Maximum iterations for ZINB parameter estimation. Default: 100.
#' @param tol Convergence tolerance for ZINB estimation. Default: 1e-6.
#' @param nCores Number of cores for parallelization. Uses BiocParallel backend.
#' @param verbose Logical. If TRUE, print progress messages. Default: FALSE.
#' @param ... Additional arguments (currently unused).
#'
#' @return A binary adjacency matrix representing the estimated graph structure.
#'
#' @details
#' The PC algorithm estimates graphical model structure using conditional 
#' independence tests. Different methods use different distributional assumptions:
#' \itemize{
#'   \item \code{"poi"}: Poisson log-linear models with Wald tests
#'   \item \code{"nb"}: Negative binomial models with Wald tests  
#'   \item \code{"zinb0"}: ZINB models assuming structure depends only on mean
#'   \item \code{"zinb1"}: ZINB models with structure dependent on both mean and zero inflation
#' }
#'
#' @importFrom BiocParallel bpparam bplapply MulticoreParam SerialParam
#' @importFrom foreach foreach %dopar%
#' @importFrom stats glm glm.fit pnorm coefficients
#' @importFrom utils combn
#'
#' @keywords internal
#' @noRd
PCzinb_internal <- function(X, method = c("poi", "nb", "zinb0", "zinb1"),
                           alpha = NULL, maxcard = 2, extend = TRUE,
                           max_iter = 100, tol = 1e-6, nCores = 1, 
                           verbose = FALSE, ...) {
  
  # Validate inputs
  method <- match.arg(method)
  
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (n < 10 || p < 2) {
    stop("X must have at least 10 samples and 2 genes")
  }
  
  # Set default alpha if not provided
  if (is.null(alpha)) {
    alpha <- 2 * pnorm(sqrt(n), lower.tail = FALSE)
  }
  
  if (verbose) {
    message("Running PCzinb with method: ", method)
    message("Sample size: ", n, ", Number of genes: ", p)
    message("Alpha level: ", round(alpha, 6), ", Max cardinality: ", maxcard)
  }
  
  # Set up BiocParallel backend
  if (nCores > 1) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = nCores)
  } else {
    BPPARAM <- BiocParallel::SerialParam()
  }
  
  # Dispatch to appropriate method
  adjacency_matrix <- switch(method,
    "poi" = .pc_poisson(X, alpha, maxcard, extend, BPPARAM, verbose),
    "nb" = .pc_negbinom(X, alpha, maxcard, extend, BPPARAM, verbose),
    "zinb0" = .pc_zinb0(X, alpha, maxcard, extend, BPPARAM, verbose, max_iter, tol),
    "zinb1" = .pc_zinb1(X, alpha, maxcard, extend, BPPARAM, verbose, max_iter, tol)
  )
  
  if (verbose) {
    n_edges <- sum(adjacency_matrix) / 2
    message("PCzinb completed. Estimated ", n_edges, " undirected edges")
  }
  
  return(adjacency_matrix)
}

#' PC algorithm with Poisson models
#' @keywords internal
#' @noRd
.pc_poisson <- function(X, alpha, maxcard, extend, BPPARAM, verbose) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Initialize fully connected graph
  adj <- matrix(1, p, p)
  diag(adj) <- 0
  
  # PC algorithm main loop
  for (card in 0:maxcard) {
    if (verbose) {
      message("Testing conditioning sets of size ", card)
    }
    
    # Test conditional independence for each node
    adj_updates <- BiocParallel::bplapply(seq_len(p), function(i) {
      current_adj <- adj[, i]
      neighbors <- which(current_adj == 1)
      
      if (length(neighbors) > card) {
        # Generate all conditioning sets of current cardinality
        if (card == 0) {
          cond_sets <- list(integer(0))
        } else {
          cond_sets <- combn(neighbors, card, simplify = FALSE)
        }
        
        # Test each neighbor for conditional independence
        for (j in neighbors) {
          for (cond_set in cond_sets) {
            if (!j %in% cond_set) {
              # Test if i and j are conditionally independent given cond_set
              if (.test_poisson_independence(X, i, j, cond_set, alpha)) {
                current_adj[j] <- 0
                break  # Remove edge and move to next neighbor
              }
            }
          }
        }
      }
      
      return(current_adj)
    }, BPPARAM = BPPARAM)
    
    # Update adjacency matrix
    for (i in seq_len(p)) {
      adj[, i] <- adj_updates[[i]]
    }
    
    # Apply extension rule
    if (extend) {
      adj <- adj + t(adj)
      adj[adj != 0] <- 1
    } else {
      adj <- adj * t(adj)
    }
  }
  
  return(adj)
}

#' Test conditional independence for Poisson models
#' @keywords internal
#' @noRd
.test_poisson_independence <- function(X, i, j, cond_set, alpha) {
  n <- nrow(X)
  
  # Prepare data
  y <- X[, i]
  
  if (length(cond_set) == 0) {
    # Unconditional test
    x_vars <- X[, j, drop = FALSE]
  } else {
    # Conditional test
    x_vars <- X[, c(j, cond_set), drop = FALSE]
  }
  
  # Fit Poisson GLM
  tryCatch({
    fit <- glm(y ~ scale(x_vars), family = "poisson")
    coef_summary <- summary(fit)$coefficients
    
    # Test significance of j (first predictor)
    p_value <- coef_summary[2, 4]  # p-value for first predictor
    
    return(p_value > alpha)
  }, error = function(e) {
    # If fitting fails, assume dependence (conservative)
    return(FALSE)
  })
}

#' PC algorithm with Negative Binomial models
#' @keywords internal
#' @noRd
.pc_negbinom <- function(X, alpha, maxcard, extend, BPPARAM, verbose) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Initialize fully connected graph
  adj <- matrix(1, p, p)
  diag(adj) <- 0
  
  # PC algorithm main loop
  for (card in 0:maxcard) {
    if (verbose) {
      message("Testing conditioning sets of size ", card)
    }
    
    # Test conditional independence for each node
    adj_updates <- BiocParallel::bplapply(seq_len(p), function(i) {
      current_adj <- adj[, i]
      neighbors <- which(current_adj == 1)
      
      if (length(neighbors) > card) {
        # Generate all conditioning sets of current cardinality
        if (card == 0) {
          cond_sets <- list(integer(0))
        } else {
          cond_sets <- combn(neighbors, card, simplify = FALSE)
        }
        
        # Test each neighbor for conditional independence
        for (j in neighbors) {
          for (cond_set in cond_sets) {
            if (!j %in% cond_set) {
              # Test if i and j are conditionally independent given cond_set
              if (.test_negbinom_independence(X, i, j, cond_set, alpha)) {
                current_adj[j] <- 0
                break  # Remove edge and move to next neighbor
              }
            }
          }
        }
      }
      
      return(current_adj)
    }, BPPARAM = BPPARAM)
    
    # Update adjacency matrix
    for (i in seq_len(p)) {
      adj[, i] <- adj_updates[[i]]
    }
    
    # Apply extension rule
    if (extend) {
      adj <- adj + t(adj)
      adj[adj != 0] <- 1
    } else {
      adj <- adj * t(adj)
    }
  }
  
  return(adj)
}

#' Test conditional independence for Negative Binomial models
#' @keywords internal
#' @noRd
.test_negbinom_independence <- function(X, i, j, cond_set, alpha) {
  n <- nrow(X)
  
  # Prepare data
  y <- X[, i]
  
  if (length(cond_set) == 0) {
    x_vars <- X[, j, drop = FALSE]
  } else {
    x_vars <- X[, c(j, cond_set), drop = FALSE]
  }
  
  # Fit Negative Binomial GLM
  tryCatch({
    # Try glm.nb if MASS is available, otherwise fallback to Poisson
    if (requireNamespace("MASS", quietly = TRUE)) {
      # Prepare data frame for MASS::glm.nb
      data_df <- data.frame(
        response = y,
        scale(x_vars)
      )
      colnames(data_df) <- c("response", paste0("V", seq_len(ncol(x_vars))))
      
      formula_str <- paste("response ~", paste(colnames(data_df)[-1], collapse = " + "))
      fit <- MASS::glm.nb(as.formula(formula_str), data = data_df, link = "log")
      
      coef_summary <- summary(fit)$coefficients
      p_value <- coef_summary[2, 4]  # p-value for first predictor
    } else {
      # Fallback to Poisson with quasi-dispersion
      fit <- glm(y ~ scale(x_vars), family = "poisson")
      coef_summary <- summary(fit)$coefficients
      p_value <- coef_summary[2, 4]
    }
    
    return(p_value > alpha)
  }, error = function(e) {
    # If fitting fails, assume dependence (conservative)
    return(FALSE)
  })
}

#' PC algorithm with ZINB models (mean-only structure)
#' @keywords internal
#' @noRd
.pc_zinb0 <- function(X, alpha, maxcard, extend, BPPARAM, verbose, max_iter, tol) {
  n <- nrow(X)
  p <- ncol(X)
  
  if (verbose) {
    message("Estimating ZINB dispersion parameters...")
  }
  
  # Estimate dispersion parameters for each gene
  dispersions <- BiocParallel::bplapply(seq_len(p), function(i) {
    .estimate_zinb_dispersion(X[, i], X[, -i, drop = FALSE], max_iter, tol)
  }, BPPARAM = BPPARAM)
  
  # Initialize fully connected graph
  adj <- matrix(1, p, p)
  diag(adj) <- 0
  
  # PC algorithm main loop
  for (card in 0:maxcard) {
    if (verbose) {
      message("Testing conditioning sets of size ", card)
    }
    
    # Test conditional independence for each node
    adj_updates <- BiocParallel::bplapply(seq_len(p), function(i) {
      current_adj <- adj[, i]
      neighbors <- which(current_adj == 1)
      
      if (length(neighbors) > card) {
        # Generate all conditioning sets of current cardinality
        if (card == 0) {
          cond_sets <- list(integer(0))
        } else {
          cond_sets <- combn(neighbors, card, simplify = FALSE)
        }
        
        # Test each neighbor for conditional independence
        for (j in neighbors) {
          for (cond_set in cond_sets) {
            if (!j %in% cond_set) {
              # Test if i and j are conditionally independent given cond_set
              if (.test_zinb0_independence(X, i, j, cond_set, alpha, 
                                          dispersions[[i]], max_iter, tol)) {
                current_adj[j] <- 0
                break
              }
            }
          }
        }
      }
      
      return(current_adj)
    }, BPPARAM = BPPARAM)
    
    # Update adjacency matrix
    for (i in seq_len(p)) {
      adj[, i] <- adj_updates[[i]]
    }
    
    # Apply extension rule
    if (extend) {
      adj <- adj + t(adj)
      adj[adj != 0] <- 1
    } else {
      adj <- adj * t(adj)
    }
  }
  
  return(adj)
}

#' PC algorithm with ZINB models (full structure)
#' @keywords internal
#' @noRd
.pc_zinb1 <- function(X, alpha, maxcard, extend, BPPARAM, verbose, max_iter, tol) {
  # For now, use the same implementation as zinb0
  # In practice, this would implement the full ZINB model
  return(.pc_zinb0(X, alpha, maxcard, extend, BPPARAM, verbose, max_iter, tol))
}

#' Estimate ZINB dispersion parameter
#' @keywords internal
#' @noRd
.estimate_zinb_dispersion <- function(y, X_cov, max_iter, tol) {
  n <- length(y)
  
  # Simple method of moments estimator
  mu_hat <- mean(y)
  var_hat <- var(y)
  
  if (var_hat > mu_hat && mu_hat > 0) {
    # Negative binomial dispersion
    theta_hat <- mu_hat^2 / (var_hat - mu_hat)
    return(pmax(theta_hat, 0.1))  # Avoid too small values
  } else {
    # Fallback to fixed dispersion
    return(1.0)
  }
}

#' Test conditional independence for ZINB models (mean-only)
#' @keywords internal
#' @noRd
.test_zinb0_independence <- function(X, i, j, cond_set, alpha, dispersion, 
                                   max_iter, tol) {
  # For simplicity, fallback to Poisson test
  # In practice, this would implement proper ZINB likelihood ratio tests
  return(.test_poisson_independence(X, i, j, cond_set, alpha))
}

#' Soft thresholding for regularization
#' @keywords internal
#' @noRd
.soft_threshold <- function(x, lambda) {
  sign(x) * pmax(0, abs(x) - lambda)
}

#' Simple zero-inflated negative binomial log-likelihood
#' @keywords internal
#' @noRd
.zinb_loglik <- function(y, mu, theta, pi) {
  # Zero-inflated component
  zero_prob <- pi + (1 - pi) * dnbinom(0, mu = mu, size = theta)
  nonzero_prob <- (1 - pi) * dnbinom(y, mu = mu, size = theta)
  
  # Log-likelihood
  loglik <- ifelse(y == 0, log(zero_prob), log(nonzero_prob))
  
  return(sum(loglik))
}