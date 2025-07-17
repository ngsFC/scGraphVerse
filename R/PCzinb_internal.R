#' PC Algorithm for Zero-Inflated Count Data - Internal Implementation
#'
#' This is an internal implementation of the PC algorithm for learning the structure
#' of graphical models from zero-inflated count data.
#'
#' @description
#' PCzinb implements the PC (Peter-Clark) algorithm for structure learning from
#' zero-inflated count data. It supports multiple distributional assumptions
#' including Poisson, negative binomial, and zero-inflated negative binomial
#' models. The algorithm performs conditional independence testing to learn
#' the structure of the underlying graphical model.
#'
#' @param X A matrix of counts (n times p) where n is the number of observations
#'   and p is the number of variables
#' @param method The algorithm used to estimate the graph: "poi" (Poisson),
#'   "nb" (negative binomial), "zinb0" (zero-inflated negative binomial focusing
#'   on mean), or "zinb1" (zero-inflated negative binomial with zero-inflation)
#' @param alpha The significance level of the tests (default: 2*pnorm(nrow(X)^0.2, lower.tail=FALSE))
#' @param maxcard The upper bound of the cardinality of the conditional sets K (default: 2)
#' @param extend If TRUE, considers the union of the tests; if FALSE, considers
#'   the intersection (default: TRUE)
#' @param max_iter Maximum number of iterations for optimization (default: 100)
#' @param tol Tolerance for convergence (default: 1e-6)
#'
#' @return The estimated adjacency matrix of the graph
#'
#' @references
#' This implementation is based on the PC algorithm for zero-inflated count data
#' described in:
#' 
#' Davide Risso, Koen Van den Berge, Ruth Heller, Sandrine Dudoit (2018).
#' "A General Framework for Flexible Microbiome Data Analysis using Zero-Inflated
#' Graphical Models"
#' 
#' Original implementation available at:
#' https://github.com/drisso/learn2count/
#'
#' @keywords internal
#' @noRd
PCzinb_internal <- function(X,
                           method = c("poi", "nb", "zinb0", "zinb1"),
                           alpha = 2 * pnorm(nrow(X)^0.2, lower.tail = FALSE),
                           maxcard = 2,
                           extend = TRUE,
                           max_iter = 100,
                           tol = 1e-6,
                           nCores = 1) {
    
    method <- match.arg(method)
    
    # Validate inputs
    if (!is.matrix(X)) {
        X <- as.matrix(X)
    }
    
    if (any(X < 0)) {
        stop("X must contain non-negative values")
    }
    
    if (alpha <= 0 || alpha >= 1) {
        stop("alpha must be between 0 and 1")
    }
    
    if (maxcard < 0) {
        stop("maxcard must be non-negative")
    }
    
    # Setup BiocParallel backend
    if (nCores == 1) {
        BPPARAM <- BiocParallel::SerialParam()
    } else {
        BPPARAM <- BiocParallel::MulticoreParam(workers = nCores)
    }
    
    # Call appropriate method
    switch(method,
           poi = pois.wald(X, maxcard, alpha, extend, BPPARAM),
           nb = nb.wald(X, maxcard, alpha, extend, BPPARAM),
           zinb0 = zinb0.noT(X, maxcard, alpha, extend, max_iter, tol, BPPARAM),
           zinb1 = zinb1.noT(X, maxcard, alpha, extend, max_iter, tol, BPPARAM))
}