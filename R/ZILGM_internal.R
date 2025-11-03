# require(MASS)
# require(glmnet)
# require(mpath)
# require(pscl)

#' Internal ZILGM Network Inference Function with BiocParallel
#'
#' This function provides the core ZILGM (Zero-Inflated Latent Gaussian Models)
#' network inference functionality with BiocParallel for parallelization.
#' It implements zero inflated graphical models for sparse count.
#'
#' @param X A matrix of expression data (samples × genes).
#' @param lambda Regularization parameter(s).
#' @param nlambda Number of lambda values to test. Default: 50.
#' @param family Distribution family: "Poisson", "NBI", or "NBII" .
#' @param update_type Algorithm type: "IRLS" or "MM". Default: "IRLS".
#' @param sym Symmetrization method: "AND" or "OR". Default: "AND".
#' @param theta Dispersion parameter for negative binomial.
#' @param thresh Threshold for coefficient sparsity. Default: 1e-6.
#' @param weights_mat Optional weight matrix for observations.
#' @param penalty_mat Optional penalty factor matrix.
#' @param do_boot Perform bootstrap stability selection. Default:FALSE.
#' @param boot_num Number of bootstrap samples. Default: 10.
#' @param beta Threshold for bootstrap stability. Default: 0.05.
#' @param lambda_min_ratio Minimum lambda ratio. Default: 1e-4.
#' @param init_select Whether to use initialization selection. Default: FALSE.
#' @param nCores Number of cores for parallelization. Uses BiocParallel backend.
#' @param verbose Verbosity level (0, 1). Default: 0.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A list containing:
#'   \item{network}{List of inferred adjacency matrices for each lambda}
#'   \item{coef_network}{Array of coefficient matrices}
#'   \item{lambda}{Vector of lambda values used}
#'   \item{call}{Function call}
#'   \item{v}{Bootstrap stability scores (if do_boot=TRUE)}
#'   \item{opt_index}{Optimal lambda index (if do_boot=TRUE)}
#'   \item{opt_lambda}{Optimal lambda value (if do_boot=TRUE)}
#'
#' @details
#' ZILGM performs network inference using zero-inflated graphical models,
#' which are particularly suitable for sparse count data.
#' The method models both the probability of zero inflation
#' and the count distribution.
#'
#' For each gene, a regularized regression is performed against all other genes,
#' with BiocParallel used to parallelize across genes. The final network is
#' constructed by combining coefficients using AND or OR operations.
#'
#' @importFrom BiocParallel bpparam bplapply MulticoreParam SerialParam
#' @importFrom Matrix Matrix
#' @importFrom mpath glmreg
#' @importFrom glmnet glmnet
#' @importFrom stats optimize optim
#' @keywords internal
#' @noRd
zilgm_internal <- function(
    X, lambda = NULL, nlambda = 50,
    family = c("Poisson", "NBI", "NBII"),
    update_type = c("IRLS", "MM"),
    sym = c("AND", "OR"), theta = NULL,
    thresh = 1e-6, weights_mat = NULL,
    penalty_mat = NULL, do_boot = FALSE,
    boot_num = 10, beta = 0.05,
    lambda_min_ratio = 1e-4,
    init_select = FALSE, nCores = 1,
    verbose = 0, ...) {
    result <- zilgm(
        X = X, lambda = lambda, nlambda = nlambda,
        family = family, update_type = update_type,
        sym = sym, theta = theta, thresh = thresh,
        weights_mat = weights_mat,
        penalty_mat = penalty_mat, do_boot = do_boot,
        boot_num = boot_num, beta = beta,
        lambda_min_ratio = lambda_min_ratio,
        init_select = init_select, nCores = nCores,
        verbose = verbose, ...
    )

    return(result)
}

zilgm <- function(
    X, lambda = NULL, nlambda = 50,
    family = c("Poisson", "NBI", "NBII"),
    update_type = c("IRLS", "MM"),
    sym = c("AND", "OR"), theta = NULL,
    thresh = 1e-6, weights_mat = NULL,
    penalty_mat = NULL, do_boot = FALSE,
    boot_num = 10, beta = 0.05,
    lambda_min_ratio = 1e-4, init_select = FALSE,
    nCores = 1, verbose = 0, ...) {
    family <- match.arg(family)
    update_type <- match.arg(update_type)
    sym <- match.arg(sym)
    fun_call <- match.call()


    if (!any(class(X) %in% "matrix")) {
        X <- as.matrix(X)
    }

    if (!any(class(X) %in% "matrix")) {
        stop("X must be a matrix")
    }

    if (any(lambda < 0)) {
        stop("lambda must be non-negative values")
    }

    n <- NROW(X)
    p <- NCOL(X)

    if (p < 2) {
        stop("X must be a matrix with 2 or more columns")
    }

    penalty <- "LASSO"

    if (verbose > 0) {
        message(
            "learning for ", family, " graphical model \n",
            "nlambda : ", nlambda, "\n",
            "penalty function : ", penalty, "\n",
            "update type : ", update_type, "\n"
        )
    }

    if (is.null(lambda)) {
        if (verbose > 0) message("\n Searching lambda \n")

        rho_max <- find_lammax(X)
        rho_min <- lambda_min_ratio * rho_max
        tmp_lams <- c(exp(seq(log(rho_max), log(rho_min),
            length = 15
        )))

        tmp_net <- zigm_network(
            X = X, lambda = tmp_lams,
            family = family,
            update_type = update_type,
            sym = sym, theta = theta,
            thresh = thresh,
            weights_mat = weights_mat,
            penalty_mat = penalty_mat,
            init_select = init_select,
            nCores = nCores, n = n, p = p,
            verbose = verbose, ...
        )

        nOfEdge <- unlist(lapply(
            tmp_net$hat_net,
            function(x) sum(x != 0)
        ))
        s_lam <- tmp_lams[which.max(nOfEdge > 1)]
        e_lam <- tmp_lams[which.max(nOfEdge)]
        lambda <- seq(s_lam, e_lam, length = nlambda)
        rm(tmp_net)
        gc()
        if (verbose > 0) message("Complete \n")
    } else {
        nlambda <- length(lambda)
    }

    out <- list()

    if (do_boot) {
        if (n < 250) {
            m <- round(0.632 * n)
        } else {
            m <- round(10 * sqrt(n))
        }

        boot_tmp <- vector(mode = "list", length = nlambda)
        for (i in seq_len(nlambda)) {
            boot_tmp[[i]] <- Matrix(0, p, p)
        }

        for (b in seq_len(boot_num)) {
            if (verbose > 0) {
                message(paste("Conducting sampling in progress : ",
                    floor(100 * (b / boot_num)), "%",
                    collapse = ""
                ))
            }

            sub_ind <- sample(seq_len(n), m, replace = FALSE)

            boot_net <- zigm_network(
                X = X[sub_ind, , drop = FALSE],
                lambda = lambda,
                family = family,
                update_type = update_type,
                sym = sym, theta = theta,
                thresh = thresh,
                weights_mat = weights_mat,
                penalty_mat = penalty_mat,
                init_select = init_select,
                nCores = nCores, n = m,
                p = p, verbose = verbose, ...
            )

            for (l in seq_len(nlambda)) {
                boot_tmp[[l]] <- boot_tmp[[l]] + boot_net$hat_net[[l]]
            }
        }

        v <- rep(0, nlambda)
        for (l in seq_len(nlambda)) {
            gv <- as.matrix(boot_tmp[[l]] / boot_num)
            gv_tmp <- 2 * gv * (1 - gv)
            v[l] <- mean(gv_tmp[upper.tri(gv_tmp)])
        }
        rm(boot_tmp)
        gc()

        opt_index <- max(which.max(v >= beta)[1] - 1, 1)
        opt_lambda <- lambda[opt_index]

        out$v <- v
        out$opt_index <- opt_index
        out$opt_lambda <- opt_lambda
    }

    net <- zigm_network(
        X = X, lambda = lambda, family = family,
        update_type = update_type, sym = sym,
        theta = theta, thresh = thresh,
        weights_mat = weights_mat,
        penalty_mat = penalty_mat,
        init_select = init_select,
        nCores = nCores, n = n, p = p,
        verbose = verbose, ...
    )

    out$network <- net$hat_net
    out$coef_network <- net$coef_net
    out$lambda <- lambda
    out$call <- fun_call
    return(out)
}



zigm_network <- function(
    X,
    lambda = NULL,
    family = c("Poisson", "NBI", "NBII"),
    update_type = c("IRLS", "MM"),
    sym = c("AND", "OR"),
    theta = NULL,
    thresh = 1e-6,
    weights_mat = NULL,
    penalty_mat = NULL,
    init_select = FALSE,
    nCores = 1,
    n,
    p,
    verbose = 0,
    ...) {
    family <- match.arg(family)
    update_type <- match.arg(update_type)
    sym <- match.arg(sym)

    coord_fun <- switch(family,
        Poisson = zilgm_poisson,
        NBI = zilgm_negbin,
        NBII = zilgm_negbin2
    )

    nlambda <- length(lambda)
    coef_mat <- array(dim = c(p, p, nlambda))

    if (is.null(weights_mat)) {
        weights_mat <- matrix(1, n, p)
    }

    if (any(weights_mat < 0)) {
        "Weights_mat must have non-negative values"
    }
    if ((NROW(weights_mat) != n) | (NCOL(weights_mat) != p)) {
        "The number of elements in weights_mat not equal to the
        number of rows and columns on X"
    }

    # Setup BiocParallel backend
    if (nCores > 1) {
        BPPARAM <- BiocParallel::MulticoreParam(workers = nCores)
    } else {
        BPPARAM <- BiocParallel::SerialParam()
    }

    coef_tmp <- BiocParallel::bplapply(seq_len(p), FUN = function(j) {
        zigm_wrapper(
            jth = j, X = X, lambda = lambda,
            family = family, update_type = update_type,
            theta = theta, thresh = thresh,
            weights = weights_mat[, j],
            penalty.factor = penalty_mat[, j],
            init_select = init_select, fun = coord_fun,
            n = n, p = p, nlambda = nlambda,
            verbose = verbose, ...
        )
    }, BPPARAM = BPPARAM)

    for (j in seq_len(p)) {
        coef_mat[, j, ] <- as.matrix(coef_tmp[[j]]$Bmat)
    }

    ghat <- lapply(seq_len(nlambda), FUN = function(l) {
        hat_net(coef_mat[, , l], thresh = thresh, type = sym)
    })
    gs <- lapply(seq_len(nlambda), FUN = function(l) as.matrix(ghat[[l]]))

    return(list(hat_net = gs, coef_net = coef_mat))
}


zigm_wrapper <- function(
    jth, X, lambda, family, update_type,
    theta, weights, penalty.factor,
    init_select, fun, n, p, nlambda,
    thresh, verbose = 0, ...) {
    seqP <- seq_len(p)
    Bmat <- Matrix(0, p, nlambda, sparse = TRUE)
    b0 <- rep(0, nlambda)

    if (init_select) {
        fit0 <- glmnet::glmnet(
            x = X[, -jth, drop = FALSE],
            y = X[, jth], standardize = FALSE,
            family = "poisson", nlambda = 100,
            dfmax = p
        )
        bic <- (1 - fit0$dev.ratio) * fit0$nulldev + 2 * fit0$df
        p0.b <- which.min(bic[-1])
        lam_ind <- p0.b
        coeff <- drop(glmnet::predict.glmnet(fit0,
            s = fit0$lambda[lam_ind],
            type = "coefficients"
        ))
        nset <- seqP[-jth][which(abs(coeff[-1]) > (thresh / 100))]

        wthres <- thresh / 100
        for (init_iter in seq_len(100)) {
            if (length(nset) == 0) {
                wthres <- wthres / 10
                nset <- seqP[-jth][which(abs(coeff[-1]) > wthres)]
            } else {
                break
            }
        }
    } else {
        nset <- seqP[-jth]
    }

    # elastic net; if alpha = 1, LASSO penalty, if alpha = 0, ridge
    if (length(nset) == 0) {
        Bmat <- Bmat
        b0 <- b0
    } else {
        for (iter in seq_len(nlambda)) {
            if (verbose == 1) {
                message(
                    "lambda = ", lambda[iter], ", ", jth, "/", p,
                    "th node learning \n"
                )
            }
            coef_res <- fun(
                x = X[, nset, drop = FALSE],
                y = X[, jth], lambda = lambda[iter],
                theta = theta, weights = weights,
                update_type = update_type,
                penalty.factor = penalty.factor,
                thresh = thresh, ...
            )

            Bmat[nset, iter] <- coef_res$bvec[-1]
            b0[iter] <- coef_res$bvec[1]
        }
    }
    return(list(b0 = b0, Bmat = Bmat))
}

dNBI <- function(y, mu, theta, log = FALSE) {
    density <- lgamma(y + theta + 1e-10) - lgamma(y + 1) -
        lgamma(theta + 1e-10) + theta *
            (log(theta + 1e-10) - log(theta + mu + 1e-10)) + y *
            (log(mu + 1e-10) - log(theta + mu + 1e-10))
    if (log == FALSE) {
        density <- exp(density)
    }
    return(density)
}

dP <- function(y, mu, log = FALSE) {
    density <- -mu + y * log(mu + 1e-10) - lgamma(y + 1)
    if (log == FALSE) {
        density <- exp(density)
    }
    return(density)
}

dNBII <- function(y, mu, sigma, log = FALSE) {
    density <- dNBI(y, mu = mu, theta = mu / sigma, log = log)
    return(density)
}



#' Find maximum lambda value for ZILGM regularization
#'
#' This function computes the maximum lambda value for regularization
#' in ZILGM network inference. It is copied from the original ZILGM package.
#'
#' @param X A matrix of expression data (samples × genes)
#' @return Maximum lambda value for regularization
#' @keywords internal
#' @noRd
find_lammax <- function(X) {
    tmp <- t(X) %*% X
    lammax <- 1 / nrow(X) * max(abs(tmp[upper.tri(tmp)]))
    return(lammax)
}

cal_ebic_inflation <- function(results, X, gamval = 1) {
    p <- ncol(X) # p by length of lambda
    n <- nrow(X) # n by length of lambda
    nlam <- ncol(results[[1]]$Bmat)
    risk <- rep(0, nlam)
    df <- rep(0, nlam)

    for (k in seq_len(nlam)) {
        ll <- matrix(0, n, p)
        dfj <- 0

        for (j in seq_len(p)) {
            eta <- results[[j]]$b0[k] + drop(X %*% results[[j]]$Bmat[, k])
            mu <- exp(eta)
            prob0 <- results[[j]]$prob0[k]
            dfj <- dfj + sum(results[[j]]$Bmat[, k] != 0)

            flag0 <- X[, j] == 0
            ll[, j] <- (1 - prob0) * exp(-mu) * mu^X[, j] / factorial(X[, j])
            ll[flag0, j] <- ll[flag0, j] + prob0
        }
        risk[k] <- sum(log(ll))
        df[k] <- dfj
    }
    risk <- -2 * risk
    npair <- p * (p - 1) / 2

    # npair = 1
    bic <- risk + log(n) * (df)
    ebic <- risk + (log(n) + 2 * gamval * log(npair)) * df

    return(list(bic = bic, ebic = ebic, risk = risk, df = df))
}

p_bvec_obj <- function(y, weights, bvec, mu, lambda, penalty.factor) {
    penalty <- lambda * sum(abs(penalty.factor * bvec[-1]))
    pnl <- -sum(weights * dP(y = y, mu = mu, log = FALSE) + 1e-10)
    return(pnl + penalty)
}

nb_bvec_obj <- function(
    y, weights, bvec, mu, theta = NULL,
    lambda, penalty.factor) {
    penalty <- lambda * sum(abs(penalty.factor * bvec[-1]))
    pnl <- -sum(weights * log(dNBI(
        y = y, theta = theta,
        mu = mu, log = FALSE
    ) + 1e-10))
    return(pnl + penalty)
}


p_objective <- function(
    y, weights, prob, bvec, mu, lambda,
    penalty.factor, posz) {
    penalty <- lambda * sum(abs(penalty.factor * bvec[-1]))
    pnl <- -sum(weights * log(prob * posz + (1 - prob) *
        dP(y = y, mu = mu, log = FALSE) + 1e-10))
    return(pnl + penalty)
}

nb_objective <- function(
    y, weights, prob, bvec, mu, theta = NULL,
    lambda, penalty.factor, posz) {
    penalty <- lambda * sum(abs(penalty.factor * bvec[-1]))
    pnl <- -sum(weights * log(prob * posz + (1 - prob) *
        dNBI(
            y = y, theta = theta, mu = mu,
            log = FALSE
        ) + 1e-10))
    return(pnl + penalty)
}

nb2_objective <- function(
    y, weights, prob, bvec, mu,
    sigma = NULL, lambda, penalty.factor, posz) {
    penalty <- lambda * sum(abs(penalty.factor * bvec[-1]))
    pnl <- -sum(weights * log(prob * posz + (1 - prob) *
        dNBII(
            y = y, sigma = sigma, mu = mu,
            log = FALSE
        ) + 1e-10))
    return(pnl + penalty)
}

# thresholding matrix
thresholding_mat <- function(mat, thres = 0.1) {
    Bmat <- Matrix(0, nrow(mat), ncol(mat), sparse = TRUE)
    for (i in seq_len(nrow(mat))) {
        flag <- abs(mat[i, ]) > thres
        Bmat[i, flag] <- mat[i, flag]
    }
    return(Bmat)
}

hat_net <- function(coef_mat, thresh = 1e-6, type = c("AND", "OR")) {
    type <- match.arg(type)

    tmp_mat <- abs(coef_mat) > thresh

    if (type == "AND") {
        res_mat <- tmp_mat * t(tmp_mat)
    }

    if (type == "OR") {
        res_mat <- (tmp_mat + t(tmp_mat) > 0) * 1
    }
    return(res_mat)
}


theta_ml <- function(y, mu, weights = NULL) {
    n <- length(y)
    if (is.null(weights)) {
        weights <- rep(1, n)
    }
    nb_theta <- function(theta, mu, y, weights) {
        return(sum(weights * dNBI(y = y, theta = theta, mu = mu, log = TRUE)))
    }
    fit <- optimize(nb_theta,
        y = y, mu = mu, weights = weights,
        interval = c(1e-4, 5e+3), maximum = TRUE
    )
    theta <- ifelse(fit$maximum > 1e+3, 1e+8, fit$maximum)
    return(theta)
}


sigma_ml <- function(y, mu, weights = NULL) {
    n <- length(y)
    if (is.null(weights)) {
        weights <- rep(1 / n, n)
    }
    NB2_theta <- function(sigma, mu, y, weights) {
        return(sum(n * weights * dNBII(
            y = y, sigma = sigma, mu = mu,
            log = TRUE
        )))
    }
    # start = c(0.01)
    fit <- optimize(NB2_theta,
        y = y, mu = mu, weights = weights,
        interval = c(1e-6, 1000), maximum = TRUE
    )
    sigma <- fit$maximum
    # sigma = ifelse(sigma <= 5e-5, 0, sigma)
    return(sigma)
}



wlasso <- function(
    X, y, eta0 = 0, wID = rep(1, nrow(X)),
    weight = rep(1, ncol(X)),
    maxStep = 1e3, eps = 1e-10, stand.scale = FALSE, trace = FALSE) {
    # trace <- TRUE
    n <- length(y)
    p <- ncol(X) # number of components
    my <- mean(y)
    y <- y - my
    mx <- apply(X, 2, mean)

    X <- X - matrix(1, nrow = n, ncol = 1) %*% mx

    if (stand.scale) {
        sdx <- apply(X, 2, sd)
        X <- sweep(X, 2, sdx, "/")
    }
    tX <- t(X)

    seqN <- seq_len(n)
    seqP <- seq_len(p)

    BetaMatr <- matrix(0, maxStep, p)
    ObjValTrace <- rep(0, maxStep)
    LamTrace <- rep(0, maxStep)

    allPredIncluded <- FALSE

    ####### Initial Solution ######
    conv <- FALSE
    Step <- 1
    wmat <- diag(wID)
    weight[abs(weight) < 1e-10] <- 1e-10

    wgcorr <- -drop(tX %*% wmat %*% y) / weight
    Beta <- rep(0, p)
    eta <- max(abs(wgcorr))

    if (eta <= eta0) {
        conv <- TRUE
        BetaMatr[1, ] <- 0
        LamTrace[1] <- eta0
        ObjValueTrace <- 1
    }

    V <- which(eta == abs(wgcorr))
    Sign <- rep(0, p)
    Sign[V] <- -sign(wgcorr[V])

    ########## Start Loop ###########
    while (!conv) {
        fx <- drop(X[, V, drop = FALSE] %*% Beta[V])
        Residual <- y - fx
        wgcorr <- -drop(t(Residual) %*% wmat %*% X) / weight

        ObjValTrace[Step] <- getObjective(Residual, wID, Beta, eta, weight)
        BetaMatr[Step, ] <- Beta
        LamTrace[Step] <- eta
        if (trace) {
            current_status(
                Step, V, Beta, eta, Sign, wgcorr, Residual,
                ObjValTrace
            )
        }
        if (eta <= eta0) {
            conv <- TRUE
            break
        } else if (Step > maxStep) {
            conv <- FALSE
            break
        }
        nV <- length(V)
        d_beta <- rderiv(
            XV = X[, V, drop = FALSE], wmat, nV,
            weight[V] * Sign[V], eps = eps
        )
        if (trace) {
            message(paste(d_beta, collapse = " "))
        }
        if (sum(abs(d_beta)) < eps) {
            break
        }
        d_wgcorr <- 1 / weight * drop(tX%*%wmat%*%X[,V,drop = FALSE]%*%d_beta)
        d_wgcorr[abs(d_wgcorr) < eps] <- 0
        Events <- Find.Event(p, V, Beta, eta, eta0, wgcorr, d_wgcorr, Sign,
            d_beta,
            trace = trace, eps = eps
        )
        EPO <- Events$EPO
        ds <- Events$ds

        d_S1_var <- Events$d_S1_var
        d_S2_var <- Events$d_S2_var

        d_lamT1 <- Events$d_lamT1
        d_lamT2 <- Events$d_lamT2

        # Update the current solutions and step size
        Beta[V] <- Beta[V] - ds * d_beta
        eta <- eta - ds
        if (EPO == 1) {
            # a predictor reduces to 0
            removedVar <- (seqP[V])[ds == d_S1_var]
            if (Sign[removedVar] > 0) {
                if (trace) {
                    message("Event 1: variable", removedVar, "Removed from V+")
                }
                removedSign <- 1
            } else if (Sign[removedVar] < 0) {
                if (trace) {
                    message("Event 1: variable", removedVar, "Removed from V-")
                }
                removedSign <- -1
            } else {
                if (trace) {
                    message("removedVar is not in V")
                }
            }
            V <- setdiff(V, removedVar)
            Sign[removedVar] <- 0
        } else if (EPO == 2) {
            temp <- min(c(seq_len(p - nV))[abs(ds - d_S2_var) < eps])
            newVar <- seqP[-V][temp]
            if (abs(d_lamT1[temp] - ds) < eps) {
                Sign[newVar] <- 1
                if (trace) {
                    message("Event 2:  New variable", newVar, "Added into V+")
                }
            } else if (abs(d_lamT2[temp] - ds) < eps) {
                Sign[newVar] <- -1
                if (trace) {
                    message("Event 2: New variable", newVar, "Added into V-")
                }
            }
            V <- sort(union(V, newVar))
        } else if (EPO == 3) {
            if (trace) {
                message("reach at ", eta0)
            }
            Step <- Step + 1
            fx <- drop(X[, V, drop = FALSE] %*% Beta[V])
            Residual <- y - fx
            ObjValTrace[Step] <- getObjective(Residual, wID, Beta, eta, weight)
            BetaMatr[Step, ] <- Beta
            LamTrace[Step] <- eta
            break
        } else {
            stop("EPO is not correct value")
        }
        Step <- Step + 1
    }

    ##################### Loop Ends  #################################
    BetaMatr <- BetaMatr[seq_len(Step), ]
    if (!stand.scale) {
        sdx <- rep(1, p)
    }
    BetaMatr <- BetaMatr %*% diag(1 / sdx, p, p)

    Beta0 <- rep(0, Step)
    for (j in seq_len(Step)) {
        Beta0[j] <- my - sum(BetaMatr[j, ] * mx)
    }

    ###
    return(list(
        Beta = BetaMatr,
        Beta0 = Beta0,
        coefficients = c(Beta0[Step], BetaMatr[Step, ]),
        fitted.values = drop(X %*% BetaMatr[Step, ]),
        LamTrace = LamTrace[seq_len(Step)],
        ObjValueTrace = ObjValTrace[seq_len(Step)],
        Step = Step,
        conv = conv
    ))
}
current_status <- function(
    Step, V, Beta, eta, Sign, wgcorr, Residual,
    ObjValTrace) {
    message(
        "V >>> ", paste(V, collapse = " "),
        "Sign[", paste(V, collapse = ","), "]",
        paste(Sign[V], collapse = " ")
    )

    message(
        "Beta[", paste(V, collapse = ","), "] >>",
        paste(Beta[V], collapse = " "), "\t L1 norm >>",
        sum(Sign[V] * Beta[V])
    )
    message("eta   >> ", eta, "\t wgcorr >> ", paste(wgcorr, collapse = " "))
    message("resid >>> ", paste(round(Residual, 3), collapse = " "))
    message("ObjValTrace >>> ", round(ObjValTrace[Step], 4))
}
getObjective <- function(res, wID, Beta, eta, weight) {
    loss <- sum(res^2 * wID) / 2 + eta * sum(weight * abs(Beta))
    return(loss)
}
# right derivatives of d_beta, d_eta
rderiv <- function(XV, wmat, nV, wSignV, eps = 1e-8) {
    A <- t(XV) %*% wmat %*% XV
    b <- -wSignV
    tmpqr <- qr(A)

    if (tmpqr$rank <- nV) {
        sol <- rep(0, nV)
        sol[nV] <- 1e5 # return delta which has a big postive value.
        # break
    } else {
        sol <- array(qr.solve(tmpqr, b)) # gamma_Active
    }
    return(sol)
}
Find.Event <- function(
    p, V, Beta, eta, eta0, wgcorr, d_wgcorr, Sign,
    d_beta, trace = FALSE, eps = 1e-10) {
    nV <- length(V)
    # 1: Active variable becomes inactive.
    index1 <- (Sign[V] > 0 & d_beta > eps) | (Sign[V] < 0 & d_beta < -eps)
    d_S1_var <- rep(Inf, nV)
    if (sum(index1) > 0) {
        d_S1_var[index1] <- (Beta[V])[index1] / d_beta[index1]
        d_S1_var[d_S1_var <= 0] <- Inf
    }
    d_S1 <- min(d_S1_var)

    # 2: Inactive variable joins the active set.
    index2 <- rep(FALSE, p)
    index2[-V] <- TRUE
    d_lamT1 <- (eta + wgcorr[index2]) / (1 + d_wgcorr[index2])
    d_lamT2 <- (eta - wgcorr[index2]) / (1 - d_wgcorr[index2])

    d_lamT1[d_lamT1 <= eps] <- Inf
    d_lamT2[d_lamT2 <= eps] <- Inf

    d_S2_var <- apply(cbind(d_lamT1, d_lamT2), 1, min)
    d_S2_which <- apply(cbind(d_lamT1, d_lamT2), 1, which.min)

    if (sum(d_S2_var == -Inf) == p - nV | sum(d_S2_var < -eps) == p - nV) {
        d_S2_var <- rep(Inf, p - nV)
    }
    if (sum(d_S2_var == Inf) == p - nV) {
        d_S2 <- Inf
    } else {
        d_S2 <- min(d_S2_var)
    }
    # 3: The generalized correlation of active variables reduces to zero
    d_S3 <- ifelse(eta - eta0 < 0, Inf, eta - eta0)

    ds <- min(c(d_S1, d_S2, d_S3))

    if (ds == Inf) {
        # cat("ds is Inf !!! therefore no further update!!!\n")
    }
    EPO <- which.min(c(d_S1, d_S2, d_S3))

    if (trace) {
        # if(1){
        message(
            "d_S1 >>> ", d_S1, "\t",
            "d_S2 >>> ", d_S2, "\t",
            "d_S3 >>> ", d_S3, "\t",
            "Event ", EPO, "occurs!!!"
        )
    }
    return(list(
        EPO = EPO,
        ds = ds,
        d_S1_var = d_S1_var,
        d_S2_var = d_S2_var,
        d_lamT1 = d_lamT1,
        d_lamT2 = d_lamT2
    ))
}


# Poisson regression with l1 regularization using MM algorithm
# glmreg_fit = mpath:::glmreg_fit
wlasso_p <- function(
    y, x, weights, penalty.factor = NULL, eta0 = NULL,
    mu0 = NULL, lambda, thresh = 1e-6, maxit = 100,
    n = NROW(x), p = NCOL(x)) {
    fun_call <- match.call()
    obj_prev <- 1e+150

    if (is.null(eta0)) {
        eta0 <- rep(0, n)
    }

    if (is.null(mu0)) {
        mu0 <- rep(1, n)
    }

    for (i in seq_len(maxit)) {
        sig <- max(n * weights * mu0)

        bobj <- wlasso(
            X = x, y = eta0 + n * weights * (y - mu0) / sig,
            eta0 = lambda / sig, wID = rep(1, n), weight = penalty.factor,
            maxStep = maxit, eps = thresh, stand.scale = FALSE,
            trace = FALSE
        )
        bvec <- bobj$coefficients
        eta <- drop(bvec[1] + x %*% bvec[-1])
        mu <- exp(eta)

        obj <- p_bvec_obj(
            y = y, weights = weights, bvec = bvec, mu = mu,
            lambda = lambda, penalty.factor = penalty.factor
        )
        if (is.nan(obj) | is.infinite(obj)) {
            obj <- obj_prev
        }
        if (abs((obj_prev - obj) / obj_prev) < thresh) {
            bvec <- bvec
            mu <- mu
            eta <- eta
            break
        } else if (obj > obj_prev + 1e-10) {
            bvec <- bvec0
            mu <- mu0
            eta <- eta0
            break
        } else {
            obj_prev <- obj
            bvec0 <- bvec
            mu0 <- mu
            eta0 <- eta
        }
    }
    return(list(bvec = bvec, mu = mu, eta = eta, iter = i))
}


# Poisson regression with l1 regularization for x with 1 columns using IRLS
glm_p <- function(
    y, x, weights, penalty.factor = NULL, eta0 = NULL, mu0 = NULL,
    lambda, thresh = 1e-6, maxit = 100, n = NROW(x), p = NCOL(x)) {
    bobj <- glm.fit(
        x = cbind(1, x), y = y, family = "poisson", intercept = TRUE,
        weights = n * weights,
        control = list(epsilon = thresh, maxit = maxit),
        etastart = eta0, mustart = mu0
    )
    bvec <- bobj$coefficients
    mu <- bobj$fitted.values
    eta <- log(mu)
    return(list(bvec = bvec, mu = mu, eta = eta, iter = 0))
}


irls_p <- function(
    y, x, weights, penalty.factor = NULL, eta0 = NULL,
    mu0 = NULL, lambda, thresh = 1e-7, maxit = 100,
    n = NROW(x), p = NCOL(x)) {
    fun_call <- match.call()
    obj_prev <- 1e+150

    if (is.null(eta0)) {
        eta0 <- rep(0, n)
    }

    if (is.null(mu0)) {
        mu0 <- rep(1, n)
    }

    for (i in seq_len(maxit)) {
        mu0[mu0 < 1e-6] <- 1e-6
        w <- mu0
        z <- eta0 + (y - mu0) / mu0

        bobj <- glmnet::glmnet(
            x = x, y = z, family = "gaussian",
            weights = w * weights,
            lambda = lambda / sum(w * weights),
            standardize = FALSE, alpha = 1,
            thresh = thresh,
            maxit = 10 * maxit, nlambda = 1
        )

        bvec <- drop(coefficients(bobj))
        eta <- drop(bvec[1] + x %*% bvec[-1])
        eta <- ifelse(eta > log(1e+4), log(1e+4), eta)
        mu <- exp(eta)

        obj <- p_bvec_obj(
            y = y, weights = weights, bvec = bvec, mu = mu,
            lambda = lambda, penalty.factor = penalty.factor
        )
        if (is.nan(obj) | is.infinite(obj)) {
            obj <- obj_prev
        }
        if (abs((obj_prev - obj) / obj_prev) < thresh) {
            bvec <- bvec
            mu <- mu
            eta <- eta
            break
        } else {
            obj_prev <- obj
            bvec0 <- bvec
            mu0 <- mu
            eta0 <- eta
        }
    }
    return(list(bvec = bvec, mu = mu, eta = eta, iter = i))
}


pglm_p_mm <- function(
    y, x, weights, penalty.factor = NULL, bvec0 = NULL,
    eta0 = NULL, mu0 = NULL, lambda, thresh = 1e-6,
    maxit = 100, n = NROW(x), p = NCOL(x)) {
    fun_call <- match.call()
    obj_prev <- 1e+150

    if (is.null(eta0)) {
        eta0 <- rep(0, n)
    }

    if (is.null(mu0)) {
        mu0 <- rep(1, n)
    }

    for (i in seq_len(maxit)) {
        sig <- max(n * weights * mu0)
        bobj <- glmnet::glmnet(
            x = x, y = eta0 + n * weights * (y - mu0) / sig,
            family = "gaussian", alpha = 1,
            lambda = lambda / sig,
            penalty.factor = penalty.factor,
            maxit = 10 * maxit,
            thresh = thresh,
            standardize = FALSE
        )
        bvec <- drop(coefficients(bobj, s = lambda / sig))
        eta <- drop(bvec[1] + x %*% bvec[-1])
        mu <- exp(eta)

        obj <- p_bvec_obj(
            y = y, weights = weights, bvec = bvec, mu = mu,
            lambda = lambda, penalty.factor = penalty.factor
        )
        if (is.nan(obj) | is.infinite(obj)) {
            obj <- obj_prev
        }
        if (abs((obj_prev - obj) / obj_prev) < thresh) {
            bvec <- bvec
            mu <- mu
            eta <- eta
            break
        } else if (obj > obj_prev + 1e-10) {
            bvec <- bvec0
            mu <- mu0
            eta <- eta0
            break
        } else {
            obj_prev <- obj
            bvec0 <- bvec
            mu0 <- mu
            eta0 <- eta
        }
    }
    return(list(bvec = bvec, mu = mu, eta = eta, iter = i))
}

pglm_p_irls <- function(
    y, x, weights, bvec0 = NULL, eta0 = NULL, mu0 = NULL,
    lambda, penalty.factor = rep(1, NCOL(x)), thresh = 1e-6,
    maxit = 1e+3, n = NROW(x), p = NCOL(x)) {
    fun_call <- match.call()
    poisson_fit <- try(
        (mpath::glmreg(
            y = y, x = x, weights = weights,
            lambda = lambda, alpha = 1, family = "poisson",
            thresh = thresh, maxit = maxit,
            penalty.factor = penalty.factor,
            start = bvec0, mustart = mu0, etastart = eta0,
            standardize = FALSE, penalty = "enet",
            x.keep = FALSE, y.keep = FALSE, trace = FALSE
        )),
        silent = FALSE
    )
    if (inherits(poisson_fit, "try-error")) {
        poisson_fit <- irls_p(
            y = y, x = x, weights = weights, lambda = lambda,
            thresh = thresh, maxit = maxit,
            penalty.factor = penalty.factor,
            eta0 = eta0, mu0 = mu0
        )
        bvec <- poisson_fit$bvec
        eta <- poisson_fit$eta
        mu <- poisson_fit$mu
    } else {
        bvec <- drop(c(poisson_fit$b0, poisson_fit$beta))
        mu <- poisson_fit$fitted.values
        eta <- log(mu)
    }

    return(list(bvec = bvec, mu = mu, eta = eta))
}

zilgm_poisson <- function(
    y, x, lambda, weights = NULL,
    update_type = c("IRLS", "MM"), penalty.factor = NULL,
    thresh = 1e-6, EM_tol = 1e-5, EM_iter = 3e+2,
    tol = 1e-6, maxit = 3e+2, theta = NULL) {
    update_type <- match.arg(update_type)
    fun_call <- match.call()
    out <- list()

    n <- NROW(x)
    p <- NCOL(x)

    if ((p == 1) & (update_type == "MM")) {
        update_type <- "onecol_MM"
    }
    if ((p == 1) & (update_type == "IRLS")) {
        update_type <- "onecol_IRLS"
    }

    update_fun <- switch(update_type,
        onecol_MM = wlasso_p,
        onecol_irls = glm_p,
        MM = pglm_p_mm,
        IRLS = pglm_p_irls
    )

    pos_zero <- (y == 0)
    pos_nzero <- !pos_zero
    z0 <- z <- rep(1e-6, n)

    if (is.null(penalty.factor)) {
        penalty.factor <- rep(1, p)
    }

    if (is.null(weights)) {
        weights <- rep(1, n)
    }

    if (length(unique(y)) == 1) {
        param <- list(
            bvec = rep(0, p + 1), prob = 0,
            pos_zero = which(pos_zero), iter = 0
        )
        return(param)
    }

    weights <- weights / sum(weights)

    # mu0 = rep(mean(y[y > 0]), n)
    mu0 <- rep(mean(y), n)
    eta0 <- log(mu0)
    bvec0 <- c(eta0[1], rep(0, p))

    prob0 <- (sum(pos_zero) - sum(exp(-mu0))) / n
    prob0 <- ifelse(prob0 < 1e-10, 1e-10, ifelse(prob0 > 1, 1, prob0))

    erisk_prev <- 1e+150

    if (sum(pos_zero) == 0) {
        sol_bvec <- update_fun(
            y = y, x = x, weights = weights, bvec0 = bvec0,
            eta0 = eta0, mu0 = mu0, lambda = lambda,
            penalty.factor = penalty.factor, thresh = tol,
            maxit = maxit, n = n, p = p
        )
        bvec <- sol_bvec$bvec
        eta <- sol_bvec$eta
        mu <- sol_bvec$mu

        prob <- prob0
        iter <- 0
        erisk <- 1e+150
    } else {
        for (iter in seq_len(EM_iter)) {
            # E-step
            tmp_z <- prob0 / (prob0 + (1 - prob0) * dP(0, mu0, log = FALSE))
            tmp_z[is.nan(tmp_z)] <- 1
            tmp_z <- ifelse(tmp_z >= (1 - 1e-6), 1 - 1e-6, tmp_z)
            z[pos_zero] <- tmp_z[pos_zero]

            prob <- sum(z) / n
            prob <- ifelse(prob < 1e-10, 1e-10, ifelse(prob > 1, 1, prob))

            # M-step
            sol_bvec <- update_fun(
                y = y, x = x, weights = weights * (1 - z),
                bvec0 = bvec0, eta0 = eta0, mu0 = mu0,
                lambda = lambda, penalty.factor = penalty.factor,
                thresh = tol, maxit = maxit, n = n, p = p
            )
            bvec <- sol_bvec$bvec
            eta <- sol_bvec$eta
            mu <- sol_bvec$mu

            erisk <- p_objective(
                y = y, weights = weights, prob = prob, bvec = bvec,
                mu = mu, lambda = lambda,
                penalty.factor = penalty.factor, posz = pos_zero
            )
            if (is.infinite(erisk) | is.nan(erisk)) {
                erisk <- 1e+8
            }

            if ((abs((erisk_prev - erisk) / (erisk_prev + 1)) < EM_tol)) {
                bvec <- bvec
                erisk <- erisk
                prob <- prob
                z <- z
                break
                # } else if (erisk > erisk_prev + 1e-10) {
                #   bvec = bvec0
                #   erisk = erisk_prev
                #   prob = prob0
                #   z = z0
                #   break
            } else {
                erisk_prev <- erisk
                bvec0 <- bvec
                eta0 <- eta
                mu0 <- mu
                prob0 <- prob
                z0 <- z
            }
        }
    }
    flag <- abs(bvec) < thresh
    bvec[flag] <- 0

    out$bvec <- bvec
    out$prob <- prob
    out$pos_zero <- which(pos_zero)
    out$iterations <- iter
    out$loglik <- erisk
    out$z <- z
    out$call <- fun_call
    class(out) <- "zilgm"
    return(out)
}


# NB regression with l1 regularization for x with 1 columns
# glmreg_fit = mpath:::glmreg_fit
wlasso_nb <- function(
    y, x, weights, penalty.factor = NULL, eta0 = NULL,
    mu0 = NULL, theta0 = NULL, lambda, thresh = 1e-6,
    maxit = 100, n = NROW(x), p = NCOL(x)) {
    fun_call <- match.call()
    obj_prev <- 1e+150

    if (is.null(eta0)) {
        eta0 <- rep(0, n)
    }

    if (is.null(mu0)) {
        mu0 <- rep(1, n)
    }

    for (i in seq_len(maxit)) {
        sig <- max(n * weights * ((1 + 1 / theta0 * y) * mu0) /
            (1 + 1 / theta0 * mu0)^2)

        bobj <- wlasso(
            X = x,
            y = eta0 + n * weights * ((y - mu0) /
                (1 + 1 / theta0 * mu0)) / sig,
            eta0 = lambda / sig, wID = rep(1, n), weight = penalty.factor,
            maxStep = maxit, eps = thresh, stand.scale = FALSE,
            trace = FALSE
        )
        bvec <- bobj$coefficients
        eta <- drop(bvec[1] + x %*% bvec[-1])
        mu <- exp(eta)

        obj <- nb_bvec_obj(
            y = y, weights = weights, bvec = bvec, mu = mu,
            lambda = lambda, penalty.factor = penalty.factor,
            theta = theta0
        )
        if (is.nan(obj) | is.infinite(obj)) {
            obj <- obj_prev
        }
        if (abs((obj_prev - obj) / obj_prev) < thresh) {
            bvec <- bvec
            mu <- mu
            eta <- eta
            break
        } else if (obj > obj_prev + 1e-10) {
            bvec <- bvec0
            mu <- mu0
            eta <- eta0
            break
        } else {
            obj_prev <- obj
            bvec0 <- bvec
            mu0 <- mu
            eta0 <- eta
        }
    }
    return(list(bvec = bvec, mu = mu, eta = eta, theta = theta0, iter = i))
}

glm_nb <- function(
    y, x, weights, penalty.factor = NULL, eta0 = NULL,
    mu0 = NULL, theta0 = NULL, lambda, thresh = 1e-6,
    maxit = 100, n = NROW(x), p = NCOL(x)) {
    bobj <- glm.fit(
        x = cbind(1, x), y = y,
        family = negative.binomial(theta = theta0),
        intercept = TRUE, weights = n * weights,
        control = list(epsilon = thresh, maxit = maxit),
        etastart = eta0, mustart = mu0
    )
    bvec <- bobj$coefficients
    mu <- bobj$fitted.values
    eta <- log(mu)
    return(list(bvec = bvec, mu = mu, eta = eta, theta = theta0, iter = 0))
}


irls_nb <- function(
    y, x, weights, penalty.factor = NULL, eta0 = NULL,
    mu0 = NULL, theta0 = NULL, lambda, thresh = 1e-6,
    maxit = 100, n = NROW(x), p = NCOL(x)) {
    fun_call <- match.call()
    obj_prev <- 1e+150

    if (is.null(eta0)) {
        eta0 <- rep(0, n)
    }

    if (is.null(mu0)) {
        mu0 <- rep(1, n)
    }

    for (i in seq_len(maxit)) {
        mu0[mu0 < 1e-6] <- 1e-6
        w <- mu0 / (1 + mu0 / theta0)
        z <- eta0 + (y - mu0) / mu0

        bobj <- try((glmnet::glmnet(
            x = x, y = z, family = "gaussian",
            weights = w * weights, lambda = lambda / sum(w * weights),
            standardize = FALSE, alpha = 1, thresh = thresh,
            maxit = 10 * maxit, nlambda = 1,
            penalty.factor = penalty.factor
        )), silent = TRUE)

        if (inherits(bobj, "try-error")) {
            bvec <- rep(0, ncol(x) + 1)
            mu <- rep(1e-8, length(y))
            eta <- log(mu)
        } else {
            bvec <- drop(coefficients(bobj))
            eta <- drop(bvec[1] + x %*% bvec[-1])
            eta <- ifelse(eta > log(1e+4), log(1e+4), eta)
            mu <- exp(eta)
        }

        obj <- nb_bvec_obj(
            y = y, weights = weights, bvec = bvec, mu = mu,
            lambda = lambda, penalty.factor = penalty.factor,
            theta = theta0
        )
        if (is.nan(obj) | is.infinite(obj)) {
            obj <- obj_prev
        }
        if (abs((obj_prev - obj) / obj_prev) < thresh) {
            bvec <- bvec
            mu <- mu
            eta <- eta
            break
        } else {
            obj_prev <- obj
            bvec0 <- bvec
            mu0 <- mu
            eta0 <- eta
        }
    }
    return(list(bvec = bvec, mu = mu, eta = eta, theta = theta0, iter = i))
}

# NB regression with l1 regularization using MM algorithm
pglm_nb_mm <- function(
    y, x, weights, penalty.factor = NULL, bvec0 = NULL,
    eta0 = NULL, mu0 = NULL, theta0 = NULL, lambda,
    thresh = 1e-6, maxit = 100, n = NROW(x), p = NCOL(x)) {
    fun_call <- match.call()
    obj_prev <- 1e+150

    if (is.null(eta0)) {
        eta0 <- rep(0, n)
    }

    if (is.null(mu0)) {
        mu0 <- rep(1, n)
    }

    for (i in seq_len(maxit)) {
        sig <- max(n * weights * ((1 + 1 / theta0 * y) * mu0) /
            (1 + 1 / theta0 * mu0)^2)
        bobj <- glmnet::glmnet(
            x = x,
            y = eta0 + n * weights * ((y - mu0) /
                (1 + 1 / theta0 * mu0)) / sig,
            family = "gaussian", alpha = 1, lambda = lambda / sig,
            penalty.factor = penalty.factor, maxit = 10 * maxit,
            thresh = thresh, standardize = FALSE
        )
        bvec <- drop(coefficients(bobj, s = lambda / sig))
        eta <- drop(bvec[1] + x %*% bvec[-1])
        mu <- exp(eta)

        obj <- nb_bvec_obj(
            y = y, weights = weights, bvec = bvec, mu = mu,
            lambda = lambda, penalty.factor = penalty.factor,
            theta = theta0
        )
        if (is.nan(obj) | is.infinite(obj)) {
            obj <- obj_prev
        }

        if (abs((obj_prev - obj) / obj_prev) < thresh) {
            bvec <- bvec
            mu <- mu
            eta <- eta
            break
        } else if (obj > obj_prev + 1e-10) {
            bvec <- bvec0
            mu <- mu0
            eta <- eta0
            break
        } else {
            obj_prev <- obj
            bvec0 <- bvec
            mu0 <- mu
            eta0 <- eta
        }
    }
    return(list(bvec = bvec, mu = mu, eta = eta, theta = theta0, iter = i))
}


# NB regression with l1 regularization using IRLS algorithm
pglm_nb_irls <- function(
    y, x, weights, theta0 = NULL, bvec0 = NULL,
    eta0 = NULL, mu0 = NULL, lambda,
    penalty.factor = rep(1, NCOL(x)), thresh = 1e-6,
    maxit = 1e+3, n = NROW(x), p = NCOL(x)) {
    fun_call <- match.call()
    negbin_fit <- try((mpath::glmreg(
        y = y, x = x, weights = weights,
        lambda = lambda, alpha = 1, theta = theta0,
        family = "negbin", thresh = thresh, maxit = maxit,
        penalty.factor = penalty.factor, start = bvec0,
        mustart = mu0, etastart = eta0, standardize = FALSE,
        penalty = "enet", x.keep = FALSE, y.keep = FALSE,
        trace = FALSE
    )), silent = TRUE)
    if (inherits(negbin_fit, "try-error")) {
        negbin_fit <- try((irls_nb(
            y = y, x = x, weights = weights,
            lambda = lambda, theta0 = theta0, thresh = thresh,
            maxit = maxit, penalty.factor = penalty.factor,
            eta0 = eta0, mu0 = mu0
        )), silent = TRUE)
        if (inherits(negbin_fit, "try-error")) {
            bvec <- rep(0, ncol(x) + 1)
            mu <- rep(1e-8, length(y))
            eta <- log(mu)
        } else {
            bvec <- negbin_fit$bvec
            eta <- negbin_fit$eta
            mu <- negbin_fit$mu
        }
    } else {
        bvec <- drop(c(negbin_fit$b0, negbin_fit$beta))
        mu <- negbin_fit$fitted.values
        eta <- log(mu)
    }

    return(list(bvec = bvec, mu = mu, eta = eta, theta = theta0))
}



zilgm_negbin <- function(
    y, x, lambda, weights = NULL,
    update_type = c("IRLS", "MM"), penalty.factor = NULL,
    thresh = 1e-6, EM_tol = 1e-5, EM_iter = 3e+2,
    tol = 1e-6, maxit = 3e+2, theta = NULL) {
    update_type <- match.arg(update_type)
    fun_call <- match.call()
    out <- list()

    n <- NROW(x)
    p <- NCOL(x)

    if ((p == 1) & (update_type == "MM")) {
        update_type <- "onecol_MM"
    }
    if ((p == 1) & (update_type == "IRLS")) {
        update_type <- "onecol_IRLS"
    }

    if (!is.null(theta)) {
        fixed_theta <- TRUE
        init_theta <- theta
    } else {
        fixed_theta <- FALSE
    }

    update_fun <- switch(update_type,
        onecol_MM = wlasso_nb,
        onecol_irls = glm_nb,
        MM = pglm_nb_mm,
        IRLS = pglm_nb_irls
    )

    pos_zero <- (y == 0)
    pos_nzero <- !pos_zero
    z <- rep(1e-6, n)

    if (is.null(penalty.factor)) {
        penalty.factor <- rep(1, p)
    }

    if (is.null(weights)) {
        weights <- rep(1, n)
    }

    if (length(unique(y)) == 1) {
        param <- list(
            bvec = rep(0, p + 1), theta = 1e+8, prob = 0,
            pos_zero = which(pos_zero), iter = 0
        )
        return(param)
    }

    weights <- weights / sum(weights)

    mu0 <- rep(mean(y[y > 0]), n)
    eta0 <- log(mu0)
    bvec0 <- c(eta0[1], rep(0, p))

    theta0 <- 1e+8
    prob0 <- (sum(pos_zero) - sum(dNBI(0,
        mu = mu0, theta = theta0,
        log = FALSE
    )))
    prob0 <- ifelse(prob0 < 1e-10, 1e-10, ifelse(prob0 > 1, 1, prob0))

    erisk_prev <- 1e+150

    if (sum(pos_zero) == 0) {
        for (iter in seq_len(EM_iter)) {
            sol_bvec <- update_fun(
                y = y, x = x, weights = weights,
                penalty.factor = penalty.factor, bvec0 = bvec0,
                eta0 = eta0, mu0 = mu0, lambda = lambda,
                theta0 = theta0, thresh = tol, maxit = maxit,
                n = n, p = p
            )
            bvec <- sol_bvec$bvec
            eta <- sol_bvec$eta
            mu <- sol_bvec$mu
            theta <- sol_bvec$theta

            if (fixed_theta) {
                theta <- init_theta
            } else {
                theta <- theta_ml(y = y, mu = mu, weights = weights)
            }

            erisk <- nb_objective(
                y = y, prob = prob0, bvec = bvec, mu = mu,
                lambda = lambda, weights = weights,
                penalty.factor = penalty.factor,
                theta = theta, posz = pos_zero
            )
            if (is.infinite(erisk) | is.nan(erisk)) {
                erisk <- erisk_prev
            }
            if ((abs((erisk_prev - erisk) / (erisk_prev + 1)) < EM_tol)) {
                bvec <- bvec
                theta <- theta
                prob <- 0
                break
            } else if (erisk > erisk_prev + 1e-10) {
                bvec <- bvec0
                theta <- theta
                prob <- 0
                break
            } else {
                erisk_prev <- erisk
                bvec0 <- bvec
                eta0 <- eta
                mu0 <- mu
                theta0 <- theta
                prob <- 0
            }
        }
    } else {
        for (iter in seq_len(EM_iter)) {
            # E-step
            tmp_z <- prob0 / (prob0 + (1 - prob0) * dNBI(0,
                theta = theta0, mu = mu0,
                log = FALSE
            ))
            tmp_z[is.nan(tmp_z)] <- 1
            tmp_z <- ifelse(tmp_z >= (1 - 1e-6), 1 - 1e-6, tmp_z)
            z[pos_zero] <- tmp_z[pos_zero]

            prob <- sum(z) / n
            prob <- ifelse(prob < 1e-10, 1e-10, ifelse(prob > 1, 1, prob))

            # M-step
            sol_bvec <- update_fun(
                y = y, x = x, weights = weights * (1 - z),
                penalty.factor = penalty.factor, bvec0 = bvec0,
                eta0 = eta0, mu0 = mu0, lambda = lambda,
                theta0 = theta0, thresh = tol, maxit = maxit,
                n = n, p = p
            )

            bvec <- sol_bvec$bvec
            eta <- sol_bvec$eta
            mu <- sol_bvec$mu

            if (fixed_theta) {
                theta <- init_theta
            } else {
                theta <- theta_ml(y = y, mu = mu, weights = weights * (1 - z))
            }

            erisk <- nb_objective(
                y = y, prob = prob, bvec = bvec, mu = mu,
                lambda = lambda, weights = weights,
                penalty.factor = penalty.factor,
                theta = theta, posz = pos_zero
            )
            if (is.infinite(erisk) | is.nan(erisk)) {
                erisk <- erisk_prev
            }

            if ((abs((erisk_prev - erisk) / (erisk_prev + 1)) < EM_tol)) {
                bvec <- bvec
                theta <- theta
                prob <- prob
                break
            } else {
                erisk_prev <- erisk
                bvec0 <- bvec
                eta0 <- eta
                mu0 <- mu
                theta0 <- theta
                prob0 <- prob
            }
        }
    }
    flag <- abs(bvec) < thresh
    bvec[flag] <- 0

    out$bvec <- bvec
    out$theta <- theta
    out$prob <- prob
    out$pos_zero <- which(pos_zero)
    out$iterations <- iter
    out$loglik <- erisk
    out$call <- fun_call
    class(out) <- "zilgm"
    return(out)
}

zilgm_negbin2 <- function(
    y, x, lambda, weights = NULL,
    update_type = c("IRLS", "MM"), penalty.factor = NULL,
    tol = 1e-6, EM_tol = 1e-5, EM_iter = 3e+2,
    thresh = 1e-6, maxit = 3e+2, theta = NULL) {
    update_type <- match.arg(update_type)
    fun_call <- match.call()
    out <- list()

    n <- NROW(x)
    p <- NCOL(x)

    if ((p == 1) & (update_type == "MM")) {
        update_type <- "onecol_MM"
    }
    if ((p == 1) & (update_type == "IRLS")) {
        update_type <- "onecol_IRLS"
    }

    if (!is.null(theta)) {
        fixed_theta <- TRUE
        init_theta <- theta
    } else {
        fixed_theta <- FALSE
    }

    update_fun <- switch(update_type,
        onecol_MM = wlasso_p,
        onecol_irls = glm_p,
        MM = pglm_p_mm,
        IRLS = pglm_p_irls
    )

    pos_zero <- (y == 0)
    pos_nzero <- !pos_zero
    z <- rep(1e-6, n)

    if (is.null(penalty.factor)) {
        penalty.factor <- rep(1, p)
    }

    if (is.null(weights)) {
        weights <- rep(1, n)
    }

    if (length(unique(y)) == 1) {
        param <- list(
            bvec = rep(0, p + 1), sigma = 0, prob = 0,
            pos_zero = which(pos_zero), iter = 0
        )
        return(param)
    }

    weights <- weights / sum(weights)

    mu0 <- rep(mean(y[y > 0]), n)
    eta0 <- log(mu0)
    bvec0 <- c(eta0[1], rep(0, p))

    # theta0 = sigma_ml(y = y, mu = mu0)
    theta0 <- 1e-4
    prob0 <- (sum(pos_zero) - sum(dNBII(0,
        mu = mu0, sigma = theta0,
        log = FALSE
    )))
    prob0 <- ifelse(prob0 < 1e-10, 1e-10, ifelse(prob0 > 1, 1, prob0))

    erisk_prev <- 1e+150

    if (sum(pos_zero) == 0) {
        sol_bvec <- update_fun(
            y = y, x = x, weights = weights,
            penalty.factor = penalty.factor,
            bvec0 = bvec0, eta0 = eta0, mu0 = mu0,
            lambda = lambda, thresh = tol, maxit = maxit,
            n = n, p = p
        )
        bvec <- sol_bvec$bvec
        eta <- sol_bvec$eta
        mu <- sol_bvec$mu

        prob <- prob0
        iter <- 0
        erisk <- erisk_prev
        theta <- theta0
    } else {
        for (iter in seq_len(EM_iter)) {
            # E-step
            tmp_z <- prob0 / (prob0 + (1 - prob0) * dNBII(0,
                sigma = theta0,
                mu = mu0, log = FALSE
            ))
            tmp_z[is.nan(tmp_z)] <- 1
            tmp_z <- ifelse(tmp_z >= (1 - 1e-6), 1 - 1e-6, tmp_z)
            z[pos_zero] <- tmp_z[pos_zero]

            prob <- sum(z) / n
            prob <- ifelse(prob < 1e-10, 1e-10, ifelse(prob > 1, 1, prob))

            # M-step
            sol_bvec <- update_fun(
                y = y, x = x, weights = weights * (1 - z),
                penalty.factor = penalty.factor,
                bvec0 = bvec0, eta0 = eta0, mu0 = mu0,
                lambda = lambda, thresh = tol, maxit = maxit,
                n = n, p = p
            )

            bvec <- sol_bvec$bvec
            eta <- sol_bvec$eta
            mu <- sol_bvec$mu

            if (fixed_theta) {
                theta <- init_theta
            } else {
                theta <- sigma_ml(y, mu = mu, weights = weights * (1 - z))
            }

            erisk <- nb2_objective(
                y = y, prob = prob, bvec = bvec, mu = mu,
                lambda = lambda, weights = weights,
                penalty.factor = penalty.factor, sigma = theta,
                posz = pos_zero
            )
            if (is.infinite(erisk) | is.nan(erisk)) {
                erisk <- erisk_prev
            }

            if ((abs((erisk_prev - erisk) / (erisk_prev + 1)) < EM_tol)) {
                bvec <- bvec
                theta <- theta
                prob <- prob
                z <- z
                break
            } else {
                erisk_prev <- erisk
                bvec0 <- bvec
                eta0 <- eta
                mu0 <- mu
                theta0 <- theta
                prob0 <- prob
            }
        }
    }
    flag <- abs(bvec) < thresh
    bvec[flag] <- 0

    out$bvec <- bvec
    out$theta <- theta
    out$prob <- prob
    out$pos_zero <- which(pos_zero)
    out$iterations <- iter
    out$loglik <- erisk
    out$call <- fun_call
    class(out) <- "zilgm"
    return(out)
}
