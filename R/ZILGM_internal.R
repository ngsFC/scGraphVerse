#' Zero-Inflated Latent Gaussian Model (ZILGM) - Internal Implementation
#'
#' This is an internal implementation of the Zero-Inflated Latent Gaussian Model
#' for gene regulatory network inference from count data.
#'
#' @description
#' Zero-Inflated Latent Gaussian Model (ZILGM) is a method for learning sparse
#' graphical models from zero-inflated count data. It handles the excess zeros
#' commonly found in single-cell RNA-seq data by modeling both the zero-inflation
#' and the underlying count distribution.
#'
#' @param X A matrix of count data with samples as rows and genes as columns
#' @param lambda A vector of regularization parameters. If NULL, automatically selected
#' @param nlambda Number of lambda values to use (default: 50)
#' @param family Distribution family: "Poisson", "NBI", or "NBII"
#' @param update_type Update method: "IRLS" or "MM"
#' @param sym Symmetrization method: "AND" or "OR"
#' @param theta Overdispersion parameter for negative binomial (if NULL, estimated)
#' @param thresh Convergence threshold (default: 1e-6)
#' @param weights_mat Weight matrix for observations
#' @param penalty_mat Penalty matrix for regularization
#' @param do_boot Whether to perform bootstrap selection (default: FALSE)
#' @param boot_num Number of bootstrap samples (default: 10)
#' @param beta Significance level for bootstrap selection (default: 0.05)
#' @param lambda_min_ratio Minimum lambda ratio (default: 1e-4)
#' @param init_select Whether to use initial variable selection (default: FALSE)
#' @param nCores Number of cores for parallel processing (default: 1)
#' @param verbose Verbosity level (default: 0)
#' @param ... Additional arguments
#'
#' @return A list containing:
#' \item{network}{List of adjacency matrices for each lambda}
#' \item{coef_network}{Array of coefficient matrices}
#' \item{lambda}{Vector of lambda values used}
#' \item{call}{Function call}
#' \item{v}{Bootstrap variability (if do_boot=TRUE)}
#' \item{opt_index}{Optimal lambda index (if do_boot=TRUE)}
#' \item{opt_lambda}{Optimal lambda value (if do_boot=TRUE)}
#'
#' @references
#' This implementation is based on the ZILGM method described in:
#' 
#' Beomjin Park and Sungkyu Jung (2020).
#' "Zero-Inflated Latent Gaussian Models for High-Dimensional Count Data."
#' 
#' Original implementation available at:
#' https://github.com/bbeomjin/ZILGM/
#'
#' @keywords internal
#' @noRd
zilgm_internal <- function(X, lambda = NULL, nlambda = 50, family = c("Poisson", "NBI", "NBII"), 
                          update_type = c("IRLS", "MM"), sym = c("AND", "OR"), theta = NULL, 
                          thresh = 1e-6, weights_mat = NULL, penalty_mat = NULL,
                          do_boot = FALSE, boot_num = 10, beta = 0.05, lambda_min_ratio = 1e-4,
                          init_select = FALSE, nCores = 1, verbose = 0, ...) {
    
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
    
    if (!is.null(lambda) && any(lambda < 0)) {
        stop("lambda must be non-negative values")
    }
    
    n <- NROW(X)
    p <- NCOL(X)
    
    if (p < 2) {
        stop("X must be a matrix with 2 or more columns")
    }
    
    penalty <- "LASSO"
    
    if (verbose > 0) {
        cat("learning for ", family, " graphical model \n",
            "nlambda : ", nlambda, "\n",
            "penalty function : ", penalty, "\n",
            "update type : ", update_type, "\n", sep = "")
    }
    
    if (is.null(lambda)) {
        if (verbose > 0) {
            cat("\n Searching lambda \n")
        }
        
        rho_max <- find_lammax(X)
        rho_min <- lambda_min_ratio * rho_max
        tmp_lams <- c(exp(seq(log(rho_max), log(rho_min), length = 15)))
        
        tmp_net <- zigm_network(X = X, lambda = tmp_lams, family = family, update_type = update_type, 
                               sym = sym, theta = theta, thresh = thresh, weights_mat = weights_mat, 
                               penalty_mat = penalty_mat, init_select = init_select, nCores = nCores, 
                               n = n, p = p, verbose = verbose, ...)
        
        nOfEdge <- unlist(lapply(tmp_net$hat_net, function(x) sum(x != 0)))
        s_lam <- tmp_lams[which.max(nOfEdge > 1)]
        e_lam <- tmp_lams[which.max(nOfEdge)]
        lambda <- seq(s_lam, e_lam, length = nlambda)
        rm(tmp_net)
        gc()
        
        if (verbose > 0) {
            cat("Complete \n")
        }
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
        for (i in 1:nlambda) {
            boot_tmp[[i]] <- Matrix::Matrix(0, p, p)
        }
        
        for (b in 1:boot_num) {
            if (verbose > 0) {
                cat(paste("Conducting sampling in progress : ", floor(100 * (b/boot_num)), "%", collapse = ""), "\r")
                flush.console()
            }
            
            sub_ind <- sample(1:n, m, replace = FALSE)
            
            boot_net <- zigm_network(X = X[sub_ind, , drop = FALSE], lambda = lambda, family = family, 
                                   update_type = update_type, sym = sym, theta = theta, thresh = thresh, 
                                   weights_mat = weights_mat, penalty_mat = penalty_mat, 
                                   init_select = init_select, nCores = nCores, n = m, p = p, verbose = verbose, ...)
            
            for (l in 1:nlambda) {
                boot_tmp[[l]] <- boot_tmp[[l]] + boot_net$hat_net[[l]]
            }
        }
        
        v <- rep(0, nlambda)
        for (l in 1:nlambda) {
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
    
    net <- zigm_network(X = X, lambda = lambda, family = family, update_type = update_type,
                       sym = sym, theta = theta, thresh = thresh, weights_mat = weights_mat, 
                       penalty_mat = penalty_mat, init_select = init_select, nCores = nCores, 
                       n = n, p = p, verbose = verbose, ...)
    
    out$network <- net$hat_net
    out$coef_network <- net$coef_net
    out$lambda <- lambda
    out$call <- fun_call
    return(out)
}