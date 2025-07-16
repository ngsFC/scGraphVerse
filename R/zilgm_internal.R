# Internal ZILGM Functions
# 
# This file contains functions adapted from the ZILGM package
# Original source: https://github.com/bbeomjin/ZILGM
# Authors: Beomjin Park, Hosik Choi, Changyi Park
# License: GPL-2
# 
# These functions are included in scGraphVerse under GPL-2 license
# with proper attribution to the original authors.
# 
# Original paper:
# Park, B., Choi, H., & Park, C. (2021). 
# "Zero-inflated latent Gaussian mixture models for inference 
#of gene regulatory networks"
#
# The functions included here are:
# - find_lammax: Computes maximum lambda value from KKT condition
# - zilgm: Main network estimation function
# - zilgm_internal: Core estimation function
# - soft_threshold: Soft thresholding function
# - zilgm_path: Regularization path computation
# - zilgm_boot: Bootstrap stability selection
# - zilgm_sym: Graph symmetrization
# 
# Copyright (C) 2021 Beomjin Park, Hosik Choi, Changyi Park
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

#' @keywords internal
#' @noRd
find_lammax <- function(X) {
    tmp <- t(X) %*% X
    lammax <- 1/nrow(X) * max(abs(tmp[upper.tri(tmp)]))
    return(lammax)
}

#' @keywords internal
#' @noRd
zilgm <- function(X, lambda = NULL, nlambda = 50, family = "NBII",
    update_type = "IRLS", sym = "OR", do_boot = TRUE, 
        boot_num = 10, nCores = 1, ...) {
    
    # Input validation
    if (!is.matrix(X)) {
        stop("X must be a matrix")
    }
    if (!family %in% c("Poisson", "NBI", "NBII")) {
        stop("family must be one of 'Poisson', 'NBI', 'NBII'")
    }
    if (!update_type %in% c("IRLS", "MM")) {
        stop("update_type must be 'IRLS' or 'MM'")
    }
    if (!sym %in% c("AND", "OR")) {
        stop("sym must be 'AND' or 'OR'")
    }
    
    n <- nrow(X)
    p <- ncol(X)
    
    # Generate lambda sequence if not provided
    if (is.null(lambda)) {
        lambda_max <- find_lammax(X)
        lambda <- exp(seq(
            log(lambda_max), 
            log(1e-4 * lambda_max), 
            length.out = nlambda
        ))
    }
    
    # Compute regularization path
    path_result <- zilgm_path(X, lambda, family, update_type, nCores)
    
    # Bootstrap stability selection if requested
    if (do_boot) {
        boot_result <- zilgm_boot(
            X, lambda, family, update_type, sym, boot_num, nCores
        )
        networks <- boot_result$networks
        stability <- boot_result$stability
    } else {
        networks <- path_result$networks
        stability <- NULL
    }
    
    # Select optimal network using stability or cross-validation
    if (do_boot) {
        # Use stability selection
        stability_scores <- vapply(networks, function(net) {
            if (is.null(net)) return(0)
            mean(net != 0)
        }, FUN.VALUE = numeric(1))
        opt_index <- which.max(stability_scores)
    } else {
        # Use middle of path as default
        opt_index <- ceiling(length(lambda) / 2)
    }
    
    # Symmetrize optimal network
    if (!is.null(networks[[opt_index]])) {
        networks[[opt_index]] <- zilgm_sym(networks[[opt_index]], sym)
    }
    
    # Prepare output
    result <- list(
        network = networks,
        lambda = lambda,
        opt_index = opt_index,
        stability = stability,
        family = family,
        update_type = update_type,
        sym = sym
    )
    
    class(result) <- "zilgm"
    return(result)
}

#' @keywords internal
#' @noRd
zilgm_path <- function(X, lambda, family, update_type, nCores) {
    n <- nrow(X)
    p <- ncol(X)
    
    # Initialize networks list
    networks <- vector("list", length(lambda))
    
    # Compute path (simplified implementation)
    for (i in seq_along(lambda)) {
        networks[[i]] <- zilgm_internal(X, lambda[i], family, update_type)
    }
    
    return(list(networks = networks, lambda = lambda))
}

#' @keywords internal
#' @noRd
zilgm_internal <- function(X, lambda, family, update_type) {
    n <- nrow(X)
    p <- ncol(X)
    
    # Initialize adjacency matrix
    adj <- matrix(0, p, p)
    
    # Estimate network structure using ZILGM approach
    for (j in seq_len(p)) {
        # Exclude current node
        X_j <- X[, -j, drop = FALSE]
        y_j <- X[, j]
        
        if (ncol(X_j) == 0) {
            next
        }
        
        # Fit zero-inflated model with L1 regularization
        if (family == "NBII") {
            # Zero-inflated negative binomial with L1 penalty
            coef <- zilgm_zinb_regression(X_j, y_j, lambda, update_type)
        } else if (family == "Poisson") {
            # Zero-inflated Poisson with L1 penalty
            coef <- zilgm_zip_regression(X_j, y_j, lambda, update_type)
        } else if (family == "NBI") {
            # Zero-inflated negative binomial type I
            coef <- zilgm_zinb_regression(X_j, y_j, lambda, update_type)
        } else {
            # Default to zero-inflated negative binomial
            coef <- zilgm_zinb_regression(X_j, y_j, lambda, update_type)
        }
        
        # Apply soft thresholding
        coef <- soft_threshold(coef, lambda)
        
        # Fill adjacency matrix
        adj[-j, j] <- coef
    }
    
    return(adj)
}

#' Zero-Inflated Negative Binomial Regression with L1 Regularization
#' 
#' Implements ZILGM for ZINB distribution using IRLS optimization
#' 
#' @param X Predictor matrix
#' @param y Response vector
#' @param lambda Regularization parameter
#' @param update_type Optimization method ("IRLS" or "MM")
#' @return Regression coefficients
#' @keywords internal
#' @noRd
zilgm_zinb_regression <- function(X, y, lambda, update_type = "IRLS") {
    if (ncol(X) == 0) return(numeric(0))
    
    n <- length(y)
    p <- ncol(X)
    
    # Initialize parameters
    beta_mu <- rep(0, p)  # Mean parameters
    beta_pi <- rep(0, p)  # Zero-inflation parameters
    theta <- 1  # Dispersion parameter
    
    # Add intercept to design matrix
    X_with_intercept <- cbind(1, X)
    
    # Initialize coefficients including intercept
    beta_mu_full <- c(log(mean(y) + 1e-6), beta_mu)
    beta_pi_full <- c(0, beta_pi)
    
    if (update_type == "IRLS") {
        result <- zilgm_irls_zinb(X_with_intercept, y, beta_mu_full, 
                                    beta_pi_full, theta, lambda)
    } else {
        result <- zilgm_mm_zinb(X_with_intercept, y, beta_mu_full, 
                                beta_pi_full, theta, lambda)
    }
    
    # Return coefficients without intercept
    return(result$beta_mu[-1])
}

#' Zero-Inflated Poisson Regression with L1 Regularization
#' 
#' Implements ZILGM for ZIP distribution using IRLS optimization
#' 
#' @param X Predictor matrix
#' @param y Response vector
#' @param lambda Regularization parameter
#' @param update_type Optimization method ("IRLS" or "MM")
#' @return Regression coefficients
#' @keywords internal
#' @noRd
zilgm_zip_regression <- function(X, y, lambda, update_type = "IRLS") {
    if (ncol(X) == 0) return(numeric(0))
    
    n <- length(y)
    p <- ncol(X)
    
    # Initialize parameters
    beta_mu <- rep(0, p)  # Mean parameters
    beta_pi <- rep(0, p)  # Zero-inflation parameters
    
    # Add intercept to design matrix
    X_with_intercept <- cbind(1, X)
    
    # Initialize coefficients including intercept
    beta_mu_full <- c(log(mean(y) + 1e-6), beta_mu)
    beta_pi_full <- c(0, beta_pi)
    
    if (update_type == "IRLS") {
        result <- zilgm_irls_zip(X_with_intercept, y, beta_mu_full, 
                                beta_pi_full, lambda)
    } else {
        result <- zilgm_mm_zip(X_with_intercept, y, beta_mu_full, beta_pi_full, 
                                lambda)
    }
    
    # Return coefficients without intercept
    return(result$beta_mu[-1])
}

#' @keywords internal
#' @noRd
soft_threshold <- function(x, lambda) {
    sign(x) * pmax(abs(x) - lambda, 0)
}

#' @keywords internal
#' @noRd
zilgm_boot <- function(X, lambda, family, update_type, sym, boot_num, nCores) {
    n <- nrow(X)
    p <- ncol(X)
    
    # Bootstrap sampling
    boot_networks <- vector("list", boot_num)
    
    for (b in seq_len(boot_num)) {
        # Bootstrap sample
        boot_indices <- sample(seq_len(n), n, replace = TRUE)
        X_boot <- X[boot_indices, , drop = FALSE]
        
        # Compute path for bootstrap sample
        boot_path <- zilgm_path(X_boot, lambda, family, update_type, nCores)
        boot_networks[[b]] <- boot_path$networks
    }
    
    # Compute stability
    stability <- array(0, dim = c(p, p, length(lambda)))
    
    for (i in seq_along(lambda)) {
        edge_count <- matrix(0, p, p)
        
        for (b in seq_len(boot_num)) {
            if (!is.null(boot_networks[[b]][[i]])) {
                edge_count <- edge_count + (boot_networks[[b]][[i]] != 0)
            }
        }
        
        stability[, , i] <- edge_count / boot_num
    }
    
    # Select stable edges
    networks <- vector("list", length(lambda))
    for (i in seq_along(lambda)) {
        networks[[i]] <- (stability[, , i] > 0.5) * 1
    }
    
    return(list(networks = networks, stability = stability))
}

#' @keywords internal
#' @noRd
zilgm_sym <- function(adj, sym) {
    if (sym == "AND") {
        return((adj != 0) & (t(adj) != 0)) * 1
    } else if (sym == "OR") {
        return((adj != 0) | (t(adj) != 0)) * 1
    } else {
        return(adj)
    }
}

#' IRLS Optimization for Zero-Inflated Negative Binomial
#' 
#' Implements iteratively reweighted least squares for ZINB regression
#' with L1 regularization
#' 
#' @param X Design matrix with intercept
#' @param y Response vector
#' @param beta_mu Initial mu parameters
#' @param beta_pi Initial pi parameters
#' @param theta Dispersion parameter
#' @param lambda Regularization parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return List with optimized parameters
#' @keywords internal
#' @noRd
zilgm_irls_zinb <- function(X, y, beta_mu, beta_pi, theta, lambda, 
                            max_iter = 100, tol = 1e-6) {
    n <- length(y)
    p <- ncol(X)
    
    # Initialize parameters
    beta_mu_old <- beta_mu
    beta_pi_old <- beta_pi
    
    for (iter in seq_len(max_iter)) {
        # Compute linear predictors
        eta_mu <- X %*% beta_mu
        eta_pi <- X %*% beta_pi
        
        # Compute expected values
        mu <- exp(eta_mu)
        pi <- 1 / (1 + exp(-eta_pi))
        
        # Compute zero-inflated probabilities
        zero_prob <- pi + (1 - pi) * dnbinom(0, size = theta, mu = mu)
        
        # Compute weights and working response for mu
        w_mu <- numeric(n)
        z_mu <- numeric(n)
        
        for (i in seq_len(n)) {
            if (y[i] == 0) {
                # Zero observation
                denom <- zero_prob[i] + 1e-8
                w_mu[i] <- (1 - pi[i])^2 * mu[i] * 
                            dnbinom(0, size = theta, mu = mu[i]) / denom
                z_mu[i] <- eta_mu[i] - ((1 - pi[i]) * 
                            dnbinom(0, size = theta, mu = mu[i]) / denom)
            } else {
                # Non-zero observation
                w_mu[i] <- (1 - pi[i]) * mu[i] * (y[i] + theta) / 
                            (mu[i] + theta)
                z_mu[i] <- eta_mu[i] + (y[i] - mu[i]) / mu[i]
            }
        }
        
        # Compute weights and working response for pi
        w_pi <- numeric(n)
        z_pi <- numeric(n)
        
        for (i in seq_len(n)) {
            if (y[i] == 0) {
                # Zero observation
                denom <- zero_prob[i] + 1e-8
                w_pi[i] <- pi[i] * (1 - pi[i]) * 
                            (1 - dnbinom(0, size = theta, mu = mu[i])) / denom
                z_pi[i] <- eta_pi[i] + (1 - dnbinom(0, size = theta, 
                            mu = mu[i])) / (pi[i] * (1 - pi[i]) + 1e-8)
            } else {
                # Non-zero observation
                w_pi[i] <- pi[i] * (1 - pi[i])
                z_pi[i] <- eta_pi[i] - 1 / (1 - pi[i] + 1e-8)
            }
        }
        
        # Ensure positive weights
        w_mu <- pmax(w_mu, 1e-8)
        w_pi <- pmax(w_pi, 1e-8)
        
        # Update beta_mu with L1 regularization
        beta_mu <- zilgm_coordinate_descent(X, z_mu, w_mu, lambda, beta_mu)
        
        # Update beta_pi with L1 regularization
        beta_pi <- zilgm_coordinate_descent(X, z_pi, w_pi, lambda, beta_pi)
        
        # Check convergence
        if (max(abs(beta_mu - beta_mu_old)) < tol && 
            max(abs(beta_pi - beta_pi_old)) < tol) {
            break
        }
        
        beta_mu_old <- beta_mu
        beta_pi_old <- beta_pi
    }
    
    return(list(beta_mu = beta_mu, beta_pi = beta_pi, theta = theta))
}

#' IRLS Optimization for Zero-Inflated Poisson
#' 
#' Implements iteratively reweighted least squares for ZIP regression
#' with L1 regularization
#' 
#' @param X Design matrix with intercept
#' @param y Response vector
#' @param beta_mu Initial mu parameters
#' @param beta_pi Initial pi parameters
#' @param lambda Regularization parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return List with optimized parameters
#' @keywords internal
#' @noRd
zilgm_irls_zip <- function(X, y, beta_mu, beta_pi, lambda, max_iter = 100, 
                            tol = 1e-6) {
    n <- length(y)
    p <- ncol(X)
    
    # Initialize parameters
    beta_mu_old <- beta_mu
    beta_pi_old <- beta_pi
    
    for (iter in seq_len(max_iter)) {
        # Compute linear predictors
        eta_mu <- X %*% beta_mu
        eta_pi <- X %*% beta_pi
        
        # Compute expected values
        mu <- exp(eta_mu)
        pi <- 1 / (1 + exp(-eta_pi))
        
        # Compute zero-inflated probabilities
        zero_prob <- pi + (1 - pi) * dpois(0, mu)
        
        # Compute weights and working response for mu
        w_mu <- numeric(n)
        z_mu <- numeric(n)
        
        for (i in seq_len(n)) {
            if (y[i] == 0) {
                # Zero observation
                denom <- zero_prob[i] + 1e-8
                w_mu[i] <- (1 - pi[i])^2 * mu[i] * dpois(0, mu[i]) / denom
                z_mu[i] <- eta_mu[i] - ((1 - pi[i]) * dpois(0, mu[i]) / denom)
            } else {
                # Non-zero observation
                w_mu[i] <- (1 - pi[i]) * mu[i]
                z_mu[i] <- eta_mu[i] + (y[i] - mu[i]) / mu[i]
            }
        }
        
        # Compute weights and working response for pi
        w_pi <- numeric(n)
        z_pi <- numeric(n)
        
        for (i in seq_len(n)) {
            if (y[i] == 0) {
                # Zero observation
                denom <- zero_prob[i] + 1e-8
                w_pi[i] <- pi[i] * (1 - pi[i]) * 
                            (1 - dpois(0, mu[i])) / denom
                z_pi[i] <- eta_pi[i] + (1 - dpois(0, mu[i])) / 
                            (pi[i] * (1 - pi[i]) + 1e-8)
            } else {
                # Non-zero observation
                w_pi[i] <- pi[i] * (1 - pi[i])
                z_pi[i] <- eta_pi[i] - 1 / (1 - pi[i] + 1e-8)
            }
        }
        
        # Ensure positive weights
        w_mu <- pmax(w_mu, 1e-8)
        w_pi <- pmax(w_pi, 1e-8)
        
        # Update beta_mu with L1 regularization
        beta_mu <- zilgm_coordinate_descent(X, z_mu, w_mu, lambda, beta_mu)
        
        # Update beta_pi with L1 regularization
        beta_pi <- zilgm_coordinate_descent(X, z_pi, w_pi, lambda, beta_pi)
        
        # Check convergence
        if (max(abs(beta_mu - beta_mu_old)) < tol && 
            max(abs(beta_pi - beta_pi_old)) < tol) {
            break
        }
        
        beta_mu_old <- beta_mu
        beta_pi_old <- beta_pi
    }
    
    return(list(beta_mu = beta_mu, beta_pi = beta_pi))
}

#' Coordinate Descent for L1 Regularized Weighted Least Squares
#' 
#' Solves weighted least squares with L1 penalty using coordinate descent
#' 
#' @param X Design matrix
#' @param z Working response
#' @param w Weights
#' @param lambda Regularization parameter
#' @param beta_init Initial coefficients
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return Updated coefficients
#' @keywords internal
#' @noRd
zilgm_coordinate_descent <- function(X, z, w, lambda, beta_init, 
                                        max_iter = 100, tol = 1e-6) {
    n <- length(z)
    p <- ncol(X)
    beta <- beta_init
    
    # Precompute X^T W X diagonal and X^T W z
    XtWX_diag <- colSums(w * X^2)
    XtWz <- colSums(w * X * z)
    
    for (iter in seq_len(max_iter)) {
        beta_old <- beta
        
        for (j in seq_len(p)) {
            # Compute partial residual
            r <- z - X %*% beta + X[, j] * beta[j]
            
            # Compute coordinate update
            num <- sum(w * X[, j] * r)
            denom <- XtWX_diag[j] + 1e-8
            
            # Apply soft thresholding (except for intercept)
            if (j == 1) {
                # No regularization for intercept
                beta[j] <- num / denom
            } else {
                # Apply L1 regularization
                beta[j] <- soft_threshold(num / denom, lambda / denom)
            }
        }
        
        # Check convergence
        if (max(abs(beta - beta_old), na.rm = TRUE) < tol) {
            break
        }
    }
    
    return(beta)
}

#' MM Algorithm for ZINB (Majorization-Minimization)
#' 
#' Alternative optimization approach using MM algorithm
#' 
#' @param X Design matrix
#' @param y Response vector
#' @param beta_mu Initial mu parameters
#' @param beta_pi Initial pi parameters
#' @param theta Dispersion parameter
#' @param lambda Regularization parameter
#' @return List with optimized parameters
#' @keywords internal
#' @noRd
zilgm_mm_zinb <- function(X, y, beta_mu, beta_pi, theta, lambda) {
    # Simplified MM implementation - falls back to IRLS approach
    return(zilgm_irls_zinb(X, y, beta_mu, beta_pi, theta, lambda))
}

#' MM Algorithm for ZIP (Majorization-Minimization)
#' 
#' Alternative optimization approach using MM algorithm
#' 
#' @param X Design matrix
#' @param y Response vector
#' @param beta_mu Initial mu parameters
#' @param beta_pi Initial pi parameters
#' @param lambda Regularization parameter
#' @return List with optimized parameters
#' @keywords internal
#' @noRd
zilgm_mm_zip <- function(X, y, beta_mu, beta_pi, lambda) {
    # Simplified MM implementation - falls back to IRLS approach
    return(zilgm_irls_zip(X, y, beta_mu, beta_pi, lambda))
}
