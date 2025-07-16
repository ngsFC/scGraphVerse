# Internal PCzinb Functions
# 
# This file contains functions adapted from the learn2count package
# Original source: https://github.com/drisso/learn2count
# Authors: Davide Risso, Chiara Romualdi
# License: GPL-2 | GPL-3
# 
# These functions are included in scGraphVerse under GPL license
# with proper attribution to the original authors.
# 
# Original paper:
# Risso, D., Romualdi, C. (2021). learn2count: Graphical Models for Count Data
# 
# The functions included here are:
# - PCzinb_internal: Main function for structure learning with ZINB models
# - zinb1_noT: ZINB1 algorithm implementation  
# - Helper functions for ZINB likelihood calculations
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

#' Internal PCzinb Function
#' 
#' This function implements PC algorithm for Zero-Inflated Negative 
#' Binomial data
#' Adapted from learn2count package
#' 
#' @param X Count matrix (samples x genes)
#' @param method Algorithm method ("zinb1", "zinb0", "poi", "nb")
#' @param maxcard Maximum cardinality of conditional sets
#' @param alpha Significance level for independence tests
#' @param extend Whether to use union/intersection of tests
#' @param max_iter Maximum iterations for optimization
#' @param tol Tolerance for convergence
#' 
#' @return Adjacency matrix
#' @keywords internal
#' @noRd
PCzinb_internal <- function(X, method = "zinb1", maxcard = 2, alpha = 0.05, 
                            extend = TRUE, max_iter = 100, tol = 1e-6) {
    
    # Input validation
    if (!is.matrix(X)) {
        X <- as.matrix(X)
    }
    
    n <- nrow(X)
    p <- ncol(X)
    
    # Calculate alpha if not provided (following original implementation)
    if (missing(alpha)) {
        alpha <- 2 * pnorm(n^0.2, lower.tail = FALSE)
    }
    
    # Call appropriate method
    if (method == "zinb1") {
        adj <- zinb1_noT(X, maxcard, alpha, extend)
    } else {
        stop("Only method='zinb1' is currently implemented")
    }
    
    return(adj)
}

#' Structure learning with zero-inflated negative binomial model
#'
#' This function estimates the adjacency matrix of a ZINB model given a matrix
#' of counts, using the optim function. Original implementation from 
#' learn2count.
#'
#' This approach assumes that the structure of the graph depends on both the
#' mean parameter and the zero inflation parameter. We call this model `zinb1`.
#'
#' @param X the matrix of counts (n times p).
#' @param maxcard the upper bound of the cardinality of the conditional sets K
#' @param alpha the significant level of the tests
#' @param extend if TRUE it considers the union of the tests, otherwise it
#'   considers the intersection.
#' @return the estimated adjacency matrix of the graph.
#' @keywords internal
#' @noRd
zinb1_noT <- function(X, maxcard, alpha, extend) {
    p <- ncol(X)
    n <- nrow(X)
    
    # Estimate dispersion parameter
    iter.theta <- 2
    stop.epsilon <- 0.0001
    
    # Sequential computation (replacing foreach due to BiocCheck requirements)
    zeta <- matrix(0, nrow = n, ncol = p)
    
    for (i in seq_len(p)) {
        iter <- 1
        local.lik <- rep(NA, iter.theta)
        
        # 1. Estimate zeta
        mean_xi <- mean(X[, i])
        var_xi <- var(X[, i])
        if (var_xi > mean_xi) {
            zeta_i <- rep(mean_xi^2 / (var_xi - mean_xi), n)
        } else {
            zeta_i <- rep(1, n)
        }
        
        # 2. Estimate parameters of ZINB model with zeta_i given by first step
        fitadd <- tryCatch({
            optim_funnoT(
                beta_mu = rep(1, p), 
                gamma_pi = rep(1, p), 
                Y = X[, i],
                X_mu = X[, -i, drop = FALSE], 
                zeta_i, 
                n
            )
        }, error = function(e) {
            fit <- glm(X[, i] ~ X[, -i, drop = FALSE], family = "poisson")
            optim_funnoT(
                beta_mu = fit$coefficients, 
                gamma_pi = rep(1, p), 
                Y = X[, i],
                X_mu = X[, -i, drop = FALSE], 
                zeta_i, 
                n
            )
        })
        
        # Calculate loglikelihood at the first iteration
        local.lik[1] <- zinb.loglik.regression(
            alpha = fitadd, 
            Y = X[, i],
            A.mu = cbind(rep(1, n), X[, -i, drop = FALSE]),
            A.pi = cbind(rep(1, n), X[, -i, drop = FALSE]),
            C.theta = zeta_i
        )
        
        for (iter in 2:iter.theta) {
            # 1. Estimate zeta with initial value alpha=fitadd 
            r <- zinb.regression.parseModel(
                alpha = fitadd, 
                A.mu = cbind(rep(1, n), X[, -i, drop = FALSE]),
                A.pi = cbind(rep(1, n), X[, -i, drop = FALSE])
            )
            
            zeta_temp <- zinbOptimizeDispersion(
                mu = r$logMu, 
                logitPi = r$logitPi,
                Y = X[, i], 
                n = n
            )
            
            # 2. Estimate parameters of ZINB model with zeta given by first step
            fitadd_temp <- optim_funnoT(
                beta_mu = fitadd[seq_len(p)],
                gamma_pi = fitadd[(p + 1):(2 * p)],
                Y = X[, i],
                X_mu = X[, -i, drop = FALSE], 
                zeta_temp, 
                n
            )
            
            local.lik[iter] <- zinb.loglik.regression(
                alpha = fitadd_temp, 
                Y = X[, i],
                A.mu = cbind(rep(1, n), X[, -i, drop = FALSE]),
                A.pi = cbind(rep(1, n), X[, -i, drop = FALSE]),
                C.theta = zeta_temp
            )
            
            if (local.lik[iter] > local.lik[iter - 1]) {
                fitadd <- fitadd_temp
                zeta_i <- zeta_temp
            } else {
                break
            }
            
            if (abs((local.lik[iter] - local.lik[iter - 1]) / 
                    local.lik[iter - 1]) < stop.epsilon) {
                break
            }
        }
        
        zeta[, i] <- zeta_i
    }
    
    # Fallback if zeta estimation fails
    if (any(is.na(zeta)) || any(zeta <= 0)) {
        for (i in seq_len(p)) {
            r <- zinb.regression.parseModel(
                alpha = rep(1, 2 * p), 
                A.mu = cbind(rep(1, n), scale(X[, -i, drop = FALSE])),
                A.pi = cbind(rep(1, n), scale(X[, -i, drop = FALSE]))
            )
            zeta[, i] <- zinbOptimizeDispersion(
                mu = r$logMu, 
                logitPi = r$logitPi,
                Y = X[, i], 
                n = n
            )
        }
    }
    
    # Estimate adjacency matrix
    adj <- matrix(1, p, p)
    diag(adj) <- 0
    ncard <- 0
    
    while (ncard <= maxcard) {
        V <- matrix(0, p, p)
        
        # Sequential computation (replacing foreach)
        for (i in seq_len(p)) {
            neighbor <- which(adj[, i] == 1)
            
            if (length(neighbor) >= ncard) {
                if (ncard == 0) {
                    condset <- list(integer(0))
                } else {
                    condset <- combn(neighbor, ncard, FUN = list, 
                                        simplify = FALSE)
                }
                
                for (j in seq_along(neighbor)) {
                    indcond <- FALSE
                    k <- 1
                    
                    while (!indcond && k <= length(condset)) {
                        if (!(neighbor[j] %in% condset[[k]])) {
                            
                            # Initial value
                            if (length(condset[[k]]) > 0) {
                                cond_vars <- c(neighbor[j], condset[[k]])
                            } else {
                                cond_vars <- neighbor[j]
                            }
                            
                            beta_mu <- tryCatch({
                                glm(X[, i] ~ scale(X[, cond_vars, 
                                                    drop = FALSE]), 
                                    family = "poisson")$coefficients
                            }, error = function(e) {
                                rep(1, length(cond_vars) + 1)
                            })
                            
                            gamma_pi <- rep(0.5, length(cond_vars) + 1)
                            
                            # Fit model with new edges
                            fitadd <- tryCatch({
                                optim_funnoT(
                                    beta_mu = beta_mu, 
                                    gamma_pi = gamma_pi, 
                                    Y = X[, i],
                                    X_mu = scale(X[, cond_vars, drop = FALSE]), 
                                    zeta[, i], 
                                    n
                                )
                            }, error = function(e) {
                                c(beta_mu, gamma_pi)
                            })
                            
                            # Calculate loglikelihood of new model
                            zinb.loglik.add <- tryCatch({
                                zinb.loglik.regression(
                                    alpha = fitadd, 
                                    Y = X[, i],
                                    A.mu = cbind(rep(1, n), 
                                                scale(X[, cond_vars, 
                                                    drop = FALSE])),
                                    A.pi = cbind(rep(1, n), 
                                                scale(X[, cond_vars, 
                                                    drop = FALSE])),
                                    C.theta = zeta[, i]
                                )
                            }, error = function(e) {
                                -Inf
                            })
                            
                            # Fit model without adding new edges
                            if (length(condset[[k]]) > 0) {
                                fitnoadd <- tryCatch({
                                    optim_funnoT(
                                        beta_mu = beta_mu[-2], 
                                        gamma_pi = gamma_pi[-1], 
                                        Y = X[, i],
                                        X_mu = scale(X[, condset[[k]], 
                                                    drop = FALSE]), 
                                        zeta[, i], 
                                        n
                                    )
                                }, error = function(e) {
                                    c(beta_mu[-2], gamma_pi[-1])
                                })
                                
                                # Calculate loglikelihood without new edges
                                zinb.loglik.noadd <- tryCatch({
                                    zinb.loglik.regression(
                                        alpha = fitnoadd, 
                                        Y = X[, i],
                                        A.mu = cbind(rep(1, n), 
                                                    scale(X[, condset[[k]], 
                                                        drop = FALSE])),
                                        A.pi = cbind(rep(1, n), 
                                                    scale(X[, condset[[k]], 
                                                        drop = FALSE])),
                                        C.theta = zeta[, i]
                                    )
                                }, error = function(e) {
                                    -Inf
                                })
                            } else {
                                fitnoadd <- tryCatch({
                                    optim_funnoT(
                                        beta_mu = beta_mu[c(1, 2)], 
                                        gamma_pi = gamma_pi, 
                                        Y = X[, i],
                                        X_mu = rep(0, n), 
                                        zeta[, i], 
                                        n
                                    )
                                }, error = function(e) {
                                    c(beta_mu[c(1, 2)], gamma_pi)
                                })
                                
                                # Calculate loglikelihood without new edges
                                zinb.loglik.noadd <- tryCatch({
                                    zinb.loglik.regression(
                                        alpha = fitadd, 
                                        Y = X[, i],
                                        A.mu = cbind(rep(1, n), rep(0, n)),
                                        A.pi = cbind(rep(1, n), rep(0, n)),
                                        C.theta = zeta[, i]
                                    )
                                }, error = function(e) {
                                    -Inf
                                })
                            }
                            
                            # Deviance tests
                            goodfit.Deviance <- 2 * (zinb.loglik.add - 
                                                        zinb.loglik.noadd)
                            
                            # Handle edge cases
                            if (is.na(goodfit.Deviance) || 
                                is.infinite(goodfit.Deviance)) {
                                p_val <- 1.0
                            } else {
                                p_val <- pchisq(goodfit.Deviance, 2, 
                                                lower.tail = FALSE)
                            }
                            
                            if (p_val > alpha) {
                                V[neighbor[j], i] <- 0
                                indcond <- TRUE
                            } else {
                                V[neighbor[j], i] <- 1
                            }
                        }
                        k <- k + 1
                    }
                    
                    # Set default value if no test performed
                    if (!indcond) {
                        V[neighbor[j], i] <- 1
                    }
                }
            }
        }
        
        if (extend) {
            adj <- V + t(V)
            adj[adj != 0] <- 1
        } else {
            adj <- V * t(V)
        }
        
        ncard <- ncard + 1
    }
    
    return(adj)
}

# Helper functions from original learn2count package

#' Optimize dispersion parameter for ZINB model
#' 
#' Find a single dispersion parameter for a count by 1-dimensional optimization
#' 
#' @param mu the vector mean of the negative binomial
#' @param logitPi the vector of logit of the probabilities of the zero component
#' @param Y the vector of counts
#' @param n length of the returned vector
#' @return vector of dispersion parameters
#' @keywords internal
#' @noRd
zinbOptimizeDispersion <- function(mu, logitPi, Y, n) {
    result <- tryCatch({
        g <- optimize(
            f = zinb.loglik.dispersion, 
            Y = Y, 
            mu = mu,
            logitPi = logitPi, 
            maximum = TRUE,
            interval = c(-100, 100)
        )
        
        zeta.op <- g$maximum
        
        zeta.ot <- tryCatch({
            optim(
                par = zeta.op, 
                fn = zinb.loglik.dispersion,
                gr = zinb.loglik.dispersion.gradient,
                mu = mu,
                logitPi = logitPi,
                Y = Y,
                control = list(fnscale = -1, trace = 0),
                method = "BFGS"
            )$par
        }, error = function(e) {
            zeta.op
        })
        
        if (!inherits(zeta.ot, "try-error")) {
            zeta <- zeta.ot
        } else {
            zeta <- zeta.op
        }
        
        rep(zeta, n)
    }, error = function(e) {
        # Fallback to method of moments
        mu_mean <- mean(mu)
        var_y <- var(Y)
        if (var_y > mu_mean) {
            rep(mu_mean^2 / (var_y - mu_mean), n)
        } else {
            rep(1, n)
        }
    })
    
    return(result)
}

#' Log(1 + exp(x)) function
#' 
#' Copied from copula package to avoid dependencies
#' 
#' @param x input vector
#' @param c0 parameter for numerical stability
#' @param c1 parameter for numerical stability  
#' @param c2 parameter for numerical stability
#' @return log(1 + exp(x))
#' @keywords internal
#' @noRd
log1pexp <- function(x, c0 = -37, c1 = 18, c2 = 33.3) {
    if (has.na <- any(ina <- is.na(x))) {
        y <- x
        x <- x[ok <- !ina]
    }
    r <- exp(x)
    if (any(i <- c0 < x & (i1 <- x <= c1)))
        r[i] <- log1p(r[i])
    if (any(i <- !i1 & (i2 <- x <= c2)))
        r[i] <- x[i] + 1/r[i]
    if (any(i3 <- !i2))
        r[i3] <- x[i3]
    if (has.na) {
        y[ok] <- r
        y
    } else {
        r
    }
}

#' Parse ZINB regression model
#'
#' Given the parameters of a ZINB regression model, this function parses the
#' model and computes the vector of log(mu), logit(pi), and the dimensions of
#' the different components of the vector of parameters.
#'
#' @param alpha the vectors of parameters c(a.mu, a.pi) concatenated
#' @param A.mu matrix of the model (default=empty)
#' @param A.pi matrix of the model (default=empty)
#' @return A list with slots logMu, logitPi, dim.alpha, start.alpha
#' @keywords internal
#' @noRd
zinb.regression.parseModel <- function(alpha, A.mu, A.pi) {
    n <- nrow(A.mu)
    logMu <- 0
    logitPi <- 0
    dim.alpha <- rep(0, 2)
    start.alpha <- rep(NA, 2)
    i <- 0
    
    j <- ncol(A.mu)
    if (j > 0) {
        logMu <- logMu + A.mu %*% alpha[(i + 1):(i + j)]
        dim.alpha[1] <- j
        start.alpha[1] <- i + 1
        i <- i + j
    }
    
    j <- ncol(A.pi)
    if (j > 0) {
        logitPi <- logitPi - A.pi %*% alpha[(i + 1):(i + j)]
        dim.alpha[2] <- j
        start.alpha[2] <- i + 1
        i <- i + j
    }
    
    return(list(
        logMu = logMu, 
        logitPi = logitPi, 
        dim.alpha = dim.alpha,
        start.alpha = start.alpha
    ))
}

#' ZINB optimization function
#' 
#' Structures both in mu and in pi
#' 
#' @param beta_mu parameters for mu
#' @param gamma_pi parameters for pi
#' @param Y response vector
#' @param X_mu design matrix for mu
#' @param zeta dispersion parameters
#' @param n sample size
#' @return optimized parameters
#' @keywords internal
#' @noRd
optim_funnoT <- function(beta_mu, gamma_pi, Y, X_mu, zeta, n) {
    result <- tryCatch({
        optim(
            fn = zinb.loglik.regression,
            gr = zinb.loglik.regression.gradient,
            par = c(beta_mu, gamma_pi),
            Y = Y, 
            A.mu = cbind(rep(1, n), X_mu),
            A.pi = cbind(rep(1, n), X_mu),
            C.theta = matrix(zeta, nrow = n, ncol = 1),
            control = list(fnscale = -1, trace = 0),
            method = "BFGS"
        )$par
    }, error = function(e) {
        c(beta_mu, gamma_pi)
    })
    
    return(result)
}

#' ZINB log-likelihood function
#' 
#' @param Y response vector
#' @param mu mean parameter
#' @param theta dispersion parameter
#' @param logitPi zero-inflation parameter
#' @return log-likelihood value
#' @keywords internal
#' @noRd
zinb.loglik <- function(Y, mu, theta, logitPi) {
    # log-probabilities of counts under the NB model
    # Handle warnings by explicit checks instead of suppressWarnings
    logPnb <- dnbinom(Y, size = theta, mu = mu, log = TRUE)
    
    # Handle potential NaN/Inf values
    if (any(is.na(logPnb))) {
        logPnb[is.na(logPnb)] <- -Inf
    }
    
    # contribution of zero inflation
    lognorm <- -log1pexp(logitPi)
    
    # log-likelihood
    result <- sum(logPnb[Y > 0]) + 
        sum(logPnb[Y == 0] + log1pexp(logitPi[Y == 0] - logPnb[Y == 0])) + 
        sum(lognorm)
    
    return(result)
}

#' ZINB log-likelihood for dispersion optimization
#' 
#' @param zeta log-dispersion parameter
#' @param Y response vector
#' @param mu mean parameter
#' @param logitPi zero-inflation parameter
#' @return log-likelihood value
#' @keywords internal
#' @noRd
zinb.loglik.dispersion <- function(zeta, Y, mu, logitPi) {
    zinb.loglik(Y, mu, exp(zeta), logitPi)
}

#' ZINB log-likelihood gradient for dispersion optimization
#' 
#' @param zeta log-dispersion parameter
#' @param Y response vector
#' @param mu mean parameter
#' @param logitPi zero-inflation parameter
#' @return gradient value
#' @keywords internal
#' @noRd
zinb.loglik.dispersion.gradient <- function(zeta, Y, mu, logitPi) {
    theta <- exp(zeta)
    
    # Check zeros in the count vector
    Y0 <- Y <= 0
    Y1 <- Y > 0
    has0 <- !is.na(match(TRUE, Y0))
    has1 <- !is.na(match(TRUE, Y1))
    
    grad <- 0
    if (has1) {
        grad <- grad + sum(theta * (digamma(Y[Y1] + theta) - 
                                    digamma(theta) + zeta - 
                                    log(mu[Y1] + theta) + 1 -
                                    (Y[Y1] + theta) / (mu[Y1] + theta)))
    }
    
    if (has0) {
        logPnb <- dnbinom(0, size = theta, mu = mu[Y0], log = TRUE)
        # Handle potential NaN/Inf values
        if (any(is.na(logPnb))) {
            logPnb[is.na(logPnb)] <- -Inf
        }
        
        grad <- grad + sum(theta * (zeta - log(mu[Y0] + theta) + 1 -
                                    theta / (mu[Y0] + theta)) / 
                            (1 + exp(logitPi[Y0] - logPnb)))
    }
    
    return(grad)
}

#' ZINB log-likelihood for regression
#' 
#' @param alpha parameter vector
#' @param Y response vector
#' @param A.mu design matrix for mu
#' @param A.pi design matrix for pi
#' @param C.theta dispersion matrix
#' @return log-likelihood value
#' @keywords internal
#' @noRd
zinb.loglik.regression <- function(alpha, Y,
    A.mu = matrix(nrow = length(Y), ncol = 0),
    A.pi = matrix(nrow = length(Y), ncol = 0),
    C.theta = matrix(0, nrow = length(Y), ncol = 1)) {
    
    # Parse the model
    r <- zinb.regression.parseModel(
        alpha = alpha,
        A.mu = A.mu,
        A.pi = A.pi
    )
    
    # Call the log likelihood function
    z <- zinb.loglik(Y, exp(r$logMu), exp(C.theta), r$logitPi)
    
    return(z)
}

#' ZINB log-likelihood gradient for regression
#' 
#' @param alpha parameter vector
#' @param Y response vector
#' @param A.mu design matrix for mu
#' @param A.pi design matrix for pi
#' @param C.theta dispersion matrix
#' @return gradient vector
#' @keywords internal
#' @noRd
zinb.loglik.regression.gradient <- function(alpha, Y,
    A.mu = matrix(nrow = length(Y), ncol = 0),
    A.pi = matrix(nrow = length(Y), ncol = 0),
    C.theta = matrix(0, nrow = length(Y), ncol = 1)) {
    
    # Parse the model
    r <- zinb.regression.parseModel(
        alpha = alpha,
        A.mu = A.mu,
        A.pi = A.pi
    )
    
    theta <- exp(C.theta)
    mu <- exp(r$logMu)
    n <- length(Y)
    
    # Check zeros in the count matrix
    Y0 <- Y <= 0
    Y1 <- Y > 0
    has0 <- !is.na(match(TRUE, Y0))
    has1 <- !is.na(match(TRUE, Y1))
    
    # Check what we need to compute
    need.wres.mu <- r$dim.alpha[1] > 0
    need.wres.pi <- r$dim.alpha[2] > 0
    
    # Compute some useful quantities
    muz <- 1 / (1 + exp(-r$logitPi))
    clogdens0 <- dnbinom(0, size = theta[Y0], mu = mu[Y0], log = TRUE)
    
    # Handle potential NaN/Inf values
    if (any(is.na(clogdens0))) {
        clogdens0[is.na(clogdens0)] <- -Inf
    }
    
    lognorm <- -r$logitPi - log1pexp(-r$logitPi)
    
    dens0 <- muz[Y0] + exp(lognorm[Y0] + clogdens0)
    
    # Compute the partial derivatives we need
    # w.r.t. mu
    if (need.wres.mu) {
        wres_mu <- numeric(length = n)
        if (has1) {
            wres_mu[Y1] <- Y[Y1] - mu[Y1] *
                (Y[Y1] + theta[Y1]) / (mu[Y1] + theta[Y1])
        }
        if (has0) {
            wres_mu[Y0] <- -exp(-log(dens0) + lognorm[Y0] + clogdens0 +
                                C.theta[Y0] - log(mu[Y0] + theta[Y0]) +
                                log(mu[Y0]))
        }
    }
    
    # w.r.t. pi
    if (need.wres.pi) {
        wres_pi <- numeric(length = n)
        if (has1) {
            wres_pi[Y1] <- muz[Y1]
        }
        if (has0) {
            wres_pi[Y0] <- -(1 - exp(clogdens0)) * muz[Y0] * 
                            (1 - muz[Y0]) / dens0
        }
    }
    
    # Make gradient
    grad <- numeric(0)
    
    # w.r.t. a_mu
    if (r$dim.alpha[1] > 0) {
        grad <- c(grad, colSums(wres_mu * A.mu))
    }
    
    # w.r.t. a_pi
    if (r$dim.alpha[2] > 0) {
        grad <- c(grad, colSums(wres_pi * A.pi))
    }
    
    return(grad)
}