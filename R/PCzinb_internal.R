#' Internal PCzinb Network Inference Function with BiocParallel
#'
#' This function provides the core PCzinb (PC algorithm for Zero-Inflated models) 
#' network inference functionality with BiocParallel support for parallelization.
#' It implements PC algorithms for count data using Poisson, Negative Binomial,
#' and Zero-Inflated Negative Binomial models.
#'
#' @param X A matrix of expression data (samples × genes).
#' @param method The algorithm to use: "poi" (Poisson), "nb" (Negative Binomial), 
#'   "zinb0" (ZINB with structure only in mu), or "zinb1" (ZINB with structure in both mu and pi).
#' @param alpha Significance level for conditional independence tests. Default: 2*pnorm(n^0.2, lower.tail=FALSE).
#' @param maxcard Maximum cardinality of conditioning sets. Default: 2.
#' @param extend If TRUE, considers union of tests; if FALSE, considers intersection. Default: TRUE.
#' @param nCores Number of cores for parallelization. Uses BiocParallel backend.
#' @param verbose Logical. If TRUE, print progress messages. Default: FALSE.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return An adjacency matrix representing the inferred network structure.
#'
#' @details
#' PCzinb performs structure learning using PC algorithms adapted for count data.
#' Different methods handle different distributional assumptions:
#' \itemize{
#'   \item "poi": Poisson distribution
#'   \item "nb": Negative Binomial distribution  
#'   \item "zinb0": Zero-Inflated NB with structure only in mean parameter
#'   \item "zinb1": Zero-Inflated NB with structure in both mean and zero-inflation parameters
#' }
#'
#' The algorithm uses conditional independence tests appropriate for each distribution
#' to determine network structure, with BiocParallel used for parallelization.
#'
#' @importFrom BiocParallel bpparam bplapply MulticoreParam SerialParam
#' @importFrom stats pnorm
#' @importFrom methods setGeneric setMethod signature standardGeneric
#' @keywords internal
#' @noRd
PCzinb_internal <- function(X, method = c("poi", "nb", "zinb0", "zinb1"),
                           alpha = NULL, maxcard = 2, extend = TRUE, 
                           nCores = 1, verbose = FALSE, ...) {
  
  method <- match.arg(method)
  
  # Set default alpha if not provided
  if (is.null(alpha)) {
    alpha <- 2 * pnorm(nrow(X)^0.2, lower.tail = FALSE)
  }
  
  # Call the appropriate method
  result <- switch(method,
                   poi = pois.wald(X, maxcard, alpha, extend, nCores),
                   nb = nb.wald(X, maxcard, alpha, extend, nCores),
                   zinb0 = zinb0.noT(X, maxcard, alpha, extend, nCores),
                   zinb1 = zinb1.noT(X, maxcard, alpha, extend, nCores))
  
  return(result)
}

#' Structure learning for count data using PC algorithms
#'
#' This function performs structure learning for count data using various PC algorithms
#' adapted for different distributional assumptions including Poisson, Negative Binomial,
#' and Zero-Inflated Negative Binomial models.
#'
#' @param x A matrix of count data or SummarizedExperiment object
#' @param ... Additional arguments passed to specific methods
#' @return An adjacency matrix or SummarizedExperiment with adjacency matrix
#' @export
#' @rdname PCzinb
setGeneric("PCzinb", function(x, ...) standardGeneric("PCzinb"))

#' @rdname PCzinb
#' @importFrom SummarizedExperiment assay assayNames `assay<-`
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors metadata `metadata<-`
#' @param ... Arguments to pass to the matrix method.
#' @param whichAssay The assay to use as input for the matrix method.
#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @examples
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(matrix(rpois(50, 5), ncol=10))
#' PCzinb(se, method="poi")
setMethod(
    f = "PCzinb",
    signature = signature(x = "SummarizedExperiment"),
    definition = function(x, whichAssay = "processed", ...){
        if(whichAssay == "processed" & !whichAssay %in% assayNames(x)) {
            warning("We recommend to use QPtransform() before learning the graph.")
            whichAssay <- 1
        }
        adj <- PCzinb(t(assay(x, whichAssay)), ...)
        metadata(x)$adj_mat <- as(adj, "dgCMatrix")
        rownames(metadata(x)$adj_mat) <- rownames(x)
        colnames(metadata(x)$adj_mat) <- rownames(x)

        return(x)
})

#' @param x the matrix of counts (n times p) or a SummarizedExperiment containing such matrix (transposed).
#' @param alpha the significant level of the tests
#' @param maxcard the uper bound of the cardinality of the conditional sets K
#' @param extend if TRUE it considers the union of the tests, otherwise it
#'   considers the intersection.
#' @param method the algorithm used to estimate the graph: `poi`, `nb`, `zinb0`,
#'   or `zinb1`. See details below.
#' @return if x is a matrix, the estimated adjacency matrix of the graph; if x
#'   is a SummarizedExperiment, a SummarizedExperiment object with the adjacency
#'   matrix as metadata.
#' @rdname PCzinb
#' @importFrom stats pnorm
#' @export
#' @examples
#' mat <- matrix(rpois(50, 5), nrow=10)
#' PCzinb(mat, method="poi")
setMethod(
    f = "PCzinb",
    signature = signature(x ="matrix"),
    definition = function(x,
                          method=c("poi", "nb", "zinb0", "zinb1"),
                          alpha=2*pnorm(nrow(x)^.2,lower.tail=FALSE),
                          maxcard=2,
                          extend=TRUE,
                          nCores=1) {

        method <- match.arg(method)

        switch(method,
               poi = pois.wald(x, maxcard, alpha, extend, nCores),
               nb = nb.wald(x, maxcard, alpha, extend, nCores),
               zinb0 = zinb0.noT(x, maxcard, alpha, extend, nCores),
               zinb1 = zinb1.noT(x, maxcard, alpha, extend, nCores))

})

#' @rdname PCzinb
#' @importFrom SingleCellExperiment `rowPair<-`
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @param ... Arguments to pass to the matrix method.
#' @param whichAssay The assay to use as input for the matrix method.
#' @export
#' @examples
#' library(SingleCellExperiment)
#' se <- SingleCellExperiment(matrix(rpois(50, 5), ncol=10))
#' se <- PCzinb(se, method="poi")
#' rowPair(se)
setMethod(
    f = "PCzinb",
    signature = signature(x = "SingleCellExperiment"),
    definition = function(x, whichAssay = "processed", ...){
        if(whichAssay == "processed" & !whichAssay %in% assayNames(x)) {
            warning("We recommend to use QPtransform() before learning the graph.")
            whichAssay <- 1
        }
        adj <- PCzinb(t(assay(x, whichAssay)), ...)
        rowPair(x) <- as(adj, "dgCMatrix")

        return(x)
})

#' Structure learning with negative binomial model using glm
#'
#' This function estimates the adjacency matrix of a NB model given a matrix of
#' counts, using the glm.nb function.
#'
#' @param X the matrix of counts (n times p).
#' @param alpha the significant level of the tests
#' @param maxcard the uper bound of the cardinality of the conditional sets K
#' @param extend if TRUE it considers the union of the tests, otherwise it
#'   considers the intersection.
#' @return the estimated adjacency matrix of the graph.
#' @export
#' @importFrom MASS glm.nb negative.binomial
#' @importFrom utils combn
#' @importFrom stats as.formula glm
nb.wald <- function(X, maxcard, alpha, extend, nCores = 1){
  p <- ncol(X)
  n <- nrow(X)
  adj <- matrix(1,p,p)
  diag(adj) <- 0
  ncard <- 0
  
  # Setup BiocParallel backend
  if (nCores > 1) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = nCores)
  } else {
    BPPARAM <- BiocParallel::SerialParam()
  }
  
  while (ncard <= maxcard) {
    adj.est <- BiocParallel::bplapply(1:p, function(i) {
      neighbor <- which (adj[, i] == 1)
      if (length(neighbor) >= ncard){
        condset <- combn(neighbor, ncard, FUN = list)
        for (j in 1: length(neighbor)){
          condset.temp <- condset
          indcond <- FALSE
          k <- 1
          while (!indcond & k <= length(condset.temp)){
            if (!(neighbor[j] %in% condset.temp[[k]])){
              # fit model with new edges c(neighbor[j]
              X_new <- scale(as.matrix(cbind(X[,c(neighbor[j], condset.temp[[k]])]),
                                 nrow=n, ncol=ncard+1))
              data <- data.frame(cbind(X[,i],X_new))
              colnames(data) <- paste("V", 1:(ncard+2), sep="")
              fmla <- as.formula(paste("V1 ~ ", paste(colnames(data)[-1], collapse= "+")))
              fitadd <- try(glm.nb(fmla, data = data,link = "log"),silent = TRUE)
              if (is(fitadd, "try-error")){
                fitadd <- glm(X[,i]~scale(X_new),family = negative.binomial(theta = 1))
              }

              ########## wald type tests

              if (summary(fitadd)$coefficients[2,4] > alpha){
                adj[neighbor[j], i] <- 0
                indcond <- TRUE
              }
            }
            k <- k+1
          }
        }
      }
      return(adj[,i])
    }, BPPARAM = BPPARAM)
    
    adj.est <- do.call(cbind, adj.est)

    if (extend == TRUE){
      adj <- adj.est + t(adj.est)
      adj[which(adj != 0)] <-1
    }else
    adj <- adj.est * t(adj.est)
    ncard <- ncard + 1
  }
  return(adj)
}

#' Structure learning with Poisson models
#'
#' This function estimates the adjacency matrix of a Poisson model given a
#' matrix of counts, using the glm function.
#'
#' @param X the matrix of counts (n times p).
#' @param alpha the significant level of the tests
#' @param maxcard the uper bound of the cardinality of the conditional sets K
#' @param extend if TRUE it considers the union of the tests, otherwise it
#'   considers the intersection.
#' @return the estimated adjacency matrix of the graph.
#' @export
#' @importFrom stats coefficients
pois.wald <- function(X, maxcard, alpha, extend, nCores = 1){
  p <- ncol(X)
  n <- nrow(X)
  adj <- matrix(1,p,p)
  diag(adj) <- 0
  ncard <- 0
  
  # Setup BiocParallel backend
  if (nCores > 1) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = nCores)
  } else {
    BPPARAM <- BiocParallel::SerialParam()
  }

  while (ncard <= maxcard) {
    V <- BiocParallel::bplapply(1:p, function(i) {
      neighbor <- which (adj[, i] == 1)
      if (length(neighbor) >= ncard){
        condset <- combn(neighbor, ncard, FUN = list)
        for (j in 1: length(neighbor)){
          condset.temp <- condset
          indcond <- FALSE
          k <- 1
          while (!indcond & k <= length(condset.temp)){
            if (!(neighbor[j] %in% condset.temp[[k]])){
              fit <- glm(X[,i] ~ scale(X[,c(neighbor[j], condset.temp[[k]])]),
                         family="poisson")
              if (coefficients(summary(fit))[2,4] > alpha){
                adj[neighbor[j], i] <- 0
                indcond <- TRUE
              }
            }
            k <- k+1
          }
        }
      }
      return(adj[,i])
    }, BPPARAM = BPPARAM)
    
    V <- do.call(cbind, V)

    if (extend == TRUE){
      adj <- V + t(V)
      adj[which(adj != 0)] <-1
    }else

      adj <- V * t(V)

    ncard <- ncard + 1

  }
  return(adj)
}

#' Structure learning with Poisson models using the Poisson K2 (PK2) algorithm
#'
#' This function finds the best fitting structure of a Poisson model given a
#' matrix of counts and topological ordering, using a given criterion ("AIC",
#' "BIC"). The PK2 algorithm is a modification of the K2 algorithm of Cooper and
#' Herskovits (1992) able to deal with Poisson data. See Nguyen et al. (2022)
#' for details.
#'
#' @references
#' Cooper, G. F. and Herskovits, E. (1992). A Bayesian method for the induction
#' of probabilistic networks from data. Machine learning, 9(4), 309–347.
#' @references
#' Nguyen, Chiogna, Risso, Banzato (2022). Guided structure learning of
#' DAGs for count data. arXiv:2206.09754.
#' @param X the matrix of counts (n times p).
#' @param maxcard the uper bound of the cardinality of the parent sets.
#' @param order the topological ordering of variables (names of nodes).
#' @param criterion the score function that measure the fitting of structures,
#'   could be "AIC" or "BIC".
#' @return a list containing the estimated adjacency matrix of the graph and a
#'   graphNEL object of the same graph.
#' @export
#' @importFrom stats coefficients glm AIC BIC
#' @import methods
#' @importClassesFrom graph graphNEL
Poisk2 <- function(X, order, criterion="BIC", maxcard) {

    if (length(colnames(X)) >0){
        nodes <- colnames(X)
    } else {
        nodes <- seq(1,dim(X)[2],1)
    }
    p <- length(nodes)
    pa_list <- list()

    ###### first auxiliary function to calculate the score
    f <- function(i, pa){
        Y <- X[,pa]
        q <- nrow(as.matrix(unique(Y)))
        if (q==0){
            fit <- glm(X[,i]~1,family="poisson")
        }else{
            fit <- glm(X[,i]~Y,family="poisson")
        }
        if (criterion =="BIC"){
            res <- -BIC(fit)
        }
        if(criterion =="AIC"){
            res <- -AIC(fit)
        }
        return (res)
    }


    ###### second auxiliary function to find new candidates for parent sets
    findmax <- function(i, parents) {
        gmax <- f(i, parents)
        z <- integer()
        if (pos==1) {return(z)} else{
            candidates <- setdiff(order[1:(pos-1)], parents)
            for (j in 1:length(candidates)) {
                pa <- c(parents,candidates[j])
                gnew <- f(i,pa)
                if (gnew > gmax) {gmax <-gnew
                z <- candidates[j]
                }}
            return(z)}}


    ### estimate the adjacency matrix
    Adj <- matrix(0,nrow = p, ncol = p)
    colnames( Adj) <- rownames( Adj)<- nodes
    # Use BiocParallel instead of foreach (assuming nCores parameter is available)
    BPPARAM <- BiocParallel::SerialParam()  # Default to serial for compatibility
    Adj.est <- BiocParallel::bplapply(seq_len(p), function(i) {
        pa_list[[i]] <- integer()
        pos <- which(order == nodes[i])
        gold <- f(i, c())
        OK <- TRUE
        counter <- 0

        while ((OK) & (length(pa_list[[i]]) < min(maxcard,pos-1)) ){
            counter <- counter + 1
            z <- findmax(i, pa_list[[i]])
            if (length(z) == 0) {
                OK <- FALSE
            }
            else {
                pa_list[[i]] <- c(pa_list[[i]], z)
            }
        }

        Adj[pa_list[[i]],i] <- 1
        return(Adj[,i])
    }, BPPARAM = BPPARAM)
    
    Adj.est <- do.call(cbind, Adj.est)
    colnames(Adj.est) <- rownames(Adj.est) <- nodes
    list(adjm = Adj.est, graphN = as(Adj.est, "graphNEL"))
}

# Find a single dispersion parameter for a count by 1-dimensional optimization of the likelihood
# Given a vector of count, this function computes a single dispersion parameter
# (log(theta)) the counts under a zero-inflated negative binomial
#' (ZINB) model. The ZINB distribution is parametrized by three
#' parameters: the mean value and the dispersion of the negative binomial
#' distribution, and the probability of the zero component.
#'
#' @param Y the vector of counts
#' @param mu the vector mean  of the negative binomial
#' @param logitPi the vector of logit of the probabilities of the zero component
#' @param n length of the returned vector
#' @examples
#' n <- 10
#' mu <- seq(10,50,length.out=n)
#' logitPi <- rnorm(1)
#' zeta <- rnorm(1)
#' Y <- rnbinom(n=n, size=exp(zeta), mu=mu)
#' learn2count:::zinbOptimizeDispersion ( mu, logitPi,Y,n)
zinbOptimizeDispersion <- function( mu, logitPi,Y,n) {

  g <- optimize(f=zinb.loglik.dispersion, Y=Y, mu=mu,
                logitPi=logitPi, maximum=TRUE,interval=c(-100,100))

  zeta.op <- g$maximum

  zeta.ot <- try(optim( par=zeta.op, fn=zinb.loglik.dispersion ,
                        gr=zinb.loglik.dispersion.gradient,mu=mu,
                        logitPi=logitPi,Y=Y,control=list(fnscale=-1,trace=0),
                        method="BFGS")$par,silent = TRUE)
  if (!is(zeta.ot, "try-error")){
    zeta <- zeta.ot
  }else{
    zeta <- zeta.op
  }


  zeta <- rep((zeta),n)
  zeta
}

# Copied on 5/14/2019 from log1pexp of the copula package (v. 0.999.19) by
# Marius Hofert, Ivan Kojadinovic, Martin Maechler, Jun Yan,
# Johanna G. Neslehova
# Copied here to avoid dependence on gsl which causes troubles.
log1pexp <- function (x, c0 = -37, c1 = 18, c2 = 33.3)
{
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
  }
  else r
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
#' @return A list with slots \code{logMu}, \code{logitPi}, \code{dim.alpha} (a
#'   vector of length 2 with the dimension of each of the vectors \code{a.mu},
#'   \code{a.pi}  in \code{alpha}), and \code{start.alpha} (a vector
#'   of length 2 with the starting indices of the 2 vectors in \code{alpha})
zinb.regression.parseModel <- function(alpha, A.mu,A.pi) {

  n <- nrow(A.mu)
  logMu <-0
  logitPi <- 0
  dim.alpha <- rep(0,2)
  start.alpha <- rep(NA,2)
  i <- 0

  j <- ncol(A.mu)
  if (j>0) {
    logMu <- logMu + A.mu %*% alpha[(i+1):(i+j)]
    dim.alpha[1] <- j
    start.alpha[1] <- i+1
    i <- i+j
  }

  j <- ncol(A.pi)
  if (j>0) {
    logitPi <- logitPi - A.pi %*% alpha[(i+1):(i+j)]
    dim.alpha[2] <- j
    start.alpha[2] <- i+1
    i <- i+j
  }


  return(list(logMu=logMu, logitPi=logitPi, dim.alpha=dim.alpha,
              start.alpha=start.alpha))
}

#####################
#########optimal functions

###zinb1: there are structures both in $mu$ and in $pi$
optim_funnoT <- function(beta_mu, gamma_pi, Y, X_mu, zeta, n) {
    optim( fn=zinb.loglik.regression,
           gr=zinb.loglik.regression.gradient,
           par=c(beta_mu, gamma_pi),
           Y=Y, A.mu=cbind(rep(1,n),X_mu),
           A.pi=cbind(rep(1,n),X_mu),
           C.theta=matrix(zeta, nrow = n, ncol = 1),
           control=list(fnscale=-1,trace=0),
           method="BFGS")$par
}
###zinb0: there are only structures in $mu$
optim_fun0noT <- function(beta_mu, gamma_pi, Y, X_mu, zeta, n) {
    optim( fn=zinb.loglik.regression,
           gr=zinb.loglik.regression.gradient,
           par=c(beta_mu, gamma_pi),
           Y=Y, A.mu=cbind(rep(1,n),X_mu),
           A.pi= matrix(rep(1,n),n,1),
           C.theta=matrix(zeta, nrow = n, ncol = 1),
           control=list(fnscale=-1,trace=0),
           method="BFGS")$par
}

zinb.loglik <- function(Y, mu, theta, logitPi) {

  # log-probabilities of counts under the NB model
  logPnb <- suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))

  # contribution of zero inflation
  lognorm <- - log1pexp(logitPi)

  # log-likelihood
  sum(logPnb[Y>0]) + sum(logPnb[Y==0] + log1pexp(logitPi[Y==0] -
                                                           logPnb[Y==0])) + sum(lognorm)
}

zinb.loglik.dispersion <- function(zeta, Y, mu, logitPi) {
  zinb.loglik(Y, mu, exp(zeta), logitPi)
}

zinb.loglik.dispersion.gradient <- function(zeta, Y, mu, logitPi) {
  theta <- exp(zeta)

  # Check zeros in the count vector
  Y0 <- Y <= 0
  Y1 <- Y > 0
  has0 <- !is.na(match(TRUE,Y0))
  has1 <- !is.na(match(TRUE,Y1))

  grad <- 0
  if (has1) {
    grad <- grad + sum( theta * (digamma(Y[Y1] + theta) - digamma(theta) +
                                   zeta - log(mu[Y1] + theta) + 1 -
                                   (Y[Y1] + theta)/(mu[Y1] + theta) ) )
  }

  if (has0) {
    logPnb <- suppressWarnings(dnbinom(0, size = theta, mu = mu[Y0],
                                       log = TRUE))
    grad <- grad + sum( theta * (zeta - log(mu[Y0] + theta) + 1 -
                                   theta/(mu[Y0] + theta)) / (1+exp(logitPi[Y0] - logPnb)))
    # *exp(- copula::log1pexp( -logPnb + logitPi[Y0])))

  }

  grad
}

zinb.loglik.regression <- function(alpha, Y,
                                   A.mu = matrix(nrow=length(Y), ncol=0),
                                   A.pi = matrix(nrow=length(Y), ncol=0),
                                   C.theta = matrix(0, nrow=length(Y), ncol=1)) {

  # Parse the model
  r <- zinb.regression.parseModel(alpha=alpha,
                                  A.mu = A.mu,
                                  A.pi = A.pi)

  # Call the log likelihood function
  z <- zinb.loglik(Y, exp(r$logMu), exp(C.theta), r$logitPi)
  #return z
  z
}

zinb.loglik.regression.gradient <- function(alpha, Y,
                                            A.mu = matrix(nrow=length(Y), ncol=0),
                                            A.pi = matrix(nrow=length(Y), ncol=0),
                                            C.theta = matrix(0, nrow=length(Y), ncol=1)) {

  # Parse the model
  r <- zinb.regression.parseModel(alpha=alpha,
                                  A.mu = A.mu,
                                  A.pi = A.pi)

  theta <- exp(C.theta)
  mu <- exp(r$logMu)
  n <- length(Y)

  # Check zeros in the count matrix
  Y0 <- Y <= 0
  Y1 <- Y > 0
  has0 <- !is.na(match(TRUE,Y0))
  has1 <- !is.na(match(TRUE,Y1))

  # Check what we need to compute,
  # depending on the variables over which we optimize
  need.wres.mu <- r$dim.alpha[1] >0
  need.wres.pi <- r$dim.alpha[2] >0

  # Compute some useful quantities
  muz <- 1/(1+exp(-r$logitPi))
  clogdens0 <- dnbinom(0, size = theta[Y0], mu = mu[Y0], log = TRUE)

  lognorm <- -r$logitPi - log1pexp(-r$logitPi)

  dens0 <- muz[Y0] + exp(lognorm[Y0] + clogdens0)

  # Compute the partial derivatives we need
  ## w.r.t. mu
  if (need.wres.mu) {
    wres_mu <- numeric(length = n)
    if (has1) {
      wres_mu[Y1] <- Y[Y1] - mu[Y1] *
        (Y[Y1] + theta[Y1])/(mu[Y1] + theta[Y1])
    }
    if (has0) {

      wres_mu[Y0] <- -exp(-log(dens0) + lognorm[Y0] + clogdens0 +
                            C.theta[Y0] - log(mu[Y0] + theta[Y0]) +
                            log(mu[Y0]))
    }
  }

  ## w.r.t. pi
  if (need.wres.pi) {
    wres_pi <- numeric(length = n)
    if (has1) {

      wres_pi[Y1] <-  muz[Y1]

    }
    if (has0) {

      wres_pi[Y0] <- -(1 - exp(clogdens0)) * muz[Y0] * (1-muz[Y0]) / dens0
    }
  }

  # Make gradient
  grad <- numeric(0)

  ## w.r.t. a_mu
  if (r$dim.alpha[1] >0) {
    istart <- r$start.alpha[1]
    iend <- r$start.alpha[1]+r$dim.alpha[1]-1
    grad <- c(grad , colSums(wres_mu * A.mu))
  }

  ## w.r.t. a_pi
  if (r$dim.alpha[2] >0) {
    istart <- r$start.alpha[2]
    iend <- r$start.alpha[2]+r$dim.alpha[2]-1
    grad <- c(grad , colSums(wres_pi * A.pi) )
  }



  grad
}

# Find a single dispersion parameter for a count by 1-dimensional optimization of the likelihood
# Given a vector of count, this function computes a single dispersion parameter
# (log(theta)) the counts under a  negative binomial
#' (NB) model. The NB distribution is parametrized by two
#' parameters: the mean value and the dispersion of the negative binomial distribution
#'
#' @param Y the vector of counts
#' @param mu the vector mean  of the negative binomial
#' @param n the length of the vector to return
#'   Note that theta is sometimes called inverse dispersion parameter (and
#'   phi=1/theta is then called the dispersion parameter). We follow the
#'   convention that the variance of the NB variable with mean mu and dispersion
#'   theta is mu + mu^2/theta.
nb.OptimizeDispersion <- function( mu,Y,n) {

  g <- optimize(f=nb.loglik.dispersion, Y=Y, mu=mu,
                maximum=TRUE,interval=c(-100,100))

  zeta.op <- g$maximum

  zeta.ot <- try(optim( par=zeta.op, fn=nb.loglik.dispersion ,
                        gr=nb.loglik.dispersion.gradient,mu=mu,
                        Y=Y,control=list(fnscale=-1,trace=0),
                        method="BFGS")$par,silent = TRUE)
  if (!is(zeta.ot, "try-error")){
    zeta <- zeta.ot
  }else{
    zeta <- zeta.op
  }


  zeta <- rep((zeta),n)
  zeta
}



#' Log-likelihood of the  negative binomial model
#' Given a vector of counts, this function computes the sum of the
#' log-probabilities of the counts under a  negative binomial
#' (NB) model. The NB distribution is parametrized by two
#' parameters: the mean value and the dispersion of the negative binomial distribution
#' @param Y the vector of counts
#' @param mu the vector of mean parameters of the negative binomial
#' @param theta the vector of dispersion parameters of the negative binomial, or
#'   a single scalar is also possible if the dispersion parameter is constant.
#'   Note that theta is sometimes called inverse dispersion parameter (and
#'   phi=1/theta is then called the dispersion parameter). We follow the
#'   convention that the variance of the NB variable with mean mu and dispersion
#'   theta is mu + mu^2/theta.
#'
#' @return the log-likelihood of the model.
#' @importFrom stats dnbinom optim optimize rbinom rnbinom runif var
nb.loglik <- function(Y, mu, theta) {

  # log-probabilities of counts under the NB model
  logPnb <- suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))

  sum(logPnb)

}


#' Log-likelihood of  negative binomial model, for a fixed
#' dispersion parameter
#'
#' Given a unique dispersion parameter and a set of counts, together with a
#' corresponding set of mean parameters,
#' this function computes the sum of the log-probabilities of the counts under
#' the NB model. The dispersion parameter is provided to the function through
#' zeta = log(theta), where theta is sometimes called the inverse dispersion
#' parameter.
#'
#' @param zeta a vector, the log of the inverse dispersion parameters of the
#'   negative binomial model
#' @param Y a vector of counts
#' @param mu a vector of mean parameters of the negative binomial
#' @return the log-likelihood of the model.
nb.loglik.dispersion <- function(zeta, Y, mu){

  nb.loglik(Y, mu, exp(zeta))

}

#' Parse ZINB regression model
#'
#' Given the parameters of a NB regression model, this function parses the
#' model and computes the vector of log(mu), and the dimensions of
#' the different components of the vector of parameters. See
#' \code{\link{nb.loglik.regression}} for details of the NB regression model
#' and its parameters.
#'
#' @param alpha the vectors of parameters c(a.mu) concatenated
#' @param A.mu matrix of the model (default=empty)
#' @return A list with slot \code{logMu},
#' @seealso \code{\link{nb.loglik.regression}}
nb.regression.parseModel <- function(alpha, A.mu) {

  n <- nrow(A.mu)
  logMu <-0
  dim.alpha <- rep(0,1)
  i <- 0

  j <- ncol(A.mu)
  if (j>0) {
    logMu <- logMu + A.mu %*% alpha[(i+1):(i+j)]
    dim.alpha[1] <- j
  }


  return(list(logMu=logMu, dim.alpha=dim.alpha))

}



#'log-likelihood of the NB regression model
#'
#' This function computes the log-likelihood of a NB regression
#' model given a vector of counts.
#'
#' @param alpha the vectors of parameters a.mu concatenated
#' @param Y the vector of counts
#' @param A.mu matrix of the model (see Details, default=empty)
#' @param C.theta matrix of the model (\eqn{log(\theta)}, default=zero)
#' @details The regression model is parametrized as follows: \deqn{log(\mu) =
#'   A_\mu * a_\mu}  \deqn{log(\theta) = C_\theta}
#'   where \eqn{\mu, \theta} are
#'   respectively the vector of mean parameters of the NB distribution,
#'    and the vector of inverse   dispersion parameters.  The
#'   log-likelihood of a vector of parameters \eqn{\alpha = a_\mu}
#' @return the log-likelihood.
nb.loglik.regression <- function(alpha, Y,
                                   A.mu = matrix(nrow=length(Y), ncol=0),
                                   C.theta = matrix(0, nrow=length(Y), ncol=1)) {

  # Parse the model
  r <- nb.regression.parseModel(alpha=alpha,
                                  A.mu = A.mu)

  # Call the log likelihood function
  z <- nb.loglik(Y, exp(r$logMu), exp(C.theta))
  #return z
  z
}

#' Gradient of the  log-likelihood of the NB regression model
#'
#' This function computes the gradient of the log-likelihood of a NB
#' regression model given a vector of counts.
#'
#' @param alpha the vectors of parameters a.mu concatenated
#' @param Y the vector of counts
#' @param A.mu matrix of the model (see Details, default=empty)
#' @param C.theta matrix of the model (see Details, default=zero)
#' @details The regression model is described in
#'   \code{\link{nb.loglik.regression}}.
#' @seealso \code{\link{nb.loglik.regression}}
#' @return The gradient of the log-likelihood.
nb.loglik.regression.gradient <- function(alpha, Y,
                                            A.mu = matrix(nrow=length(Y), ncol=0),
                                            C.theta = matrix(0, nrow=length(Y), ncol=1)) {

  # Parse the model
  r <- nb.regression.parseModel(alpha=alpha,
                                  A.mu = A.mu)

  theta <- exp(C.theta)
  mu <- exp(r$logMu)
  n <- length(Y)

  # Check what we need to compute,
  # depending on the variables over which we optimize
  need.wres.mu <- r$dim.alpha[1] >0

  # Compute the partial derivatives we need
  ## w.r.t. mu
  if (need.wres.mu) {
    wres_mu <- numeric(length = n)
    wres_mu <- Y - mu *
      (Y + theta)/(mu + theta)
    wres_mu <- as.vector(wres_mu)
  }


  # Make gradient
  grad <- numeric(0)

  ## w.r.t. a_mu
  if (r$dim.alpha[1] >0) {
    grad <- c(grad , colSums(wres_mu * A.mu) )
  }



  grad
}

nb.loglik.dispersion.gradient <- function(zeta, Y, mu) {
  theta <- exp(zeta)

  grad <- 0
  grad <- grad + sum( theta * (digamma(Y + theta) - digamma(theta) +
                                 zeta - log(mu + theta) + 1 -
                                 (Y + theta)/(mu + theta) ) )

  grad
}


nb.optim_funnoT <- function(beta_mu, Y, X_mu, zeta, n) {
  optim( fn=nb.loglik.regression,
         gr=nb.loglik.regression.gradient,
         par=beta_mu,
         Y=Y, A.mu=cbind(rep(1,n),X_mu),
         C.theta=matrix(zeta, nrow = n, ncol = 1),
         control=list(fnscale=-1,trace=0),
         method="BFGS")$par
}

#' Prediction scores for estimated graph
#'
#' This function returns a set of metrics useful to evaluate the goodness of the
#' estimation of the graph on simulated data, or whenever a true graph is
#' available.
#'
#' @return A vector with the number of true positives (TP), the number of false
#'   positives (FP), the number of false negatives (FN), the positive predictive
#'   value (PPV), the sensitivity (Se), and the F1 score (F1).
#'
#' @param trueG the adjacency matrix of the true graph.
#' @param estimatedG the adjacency matrix of the estimated graph.
#' @param type the type of the true graph, could be "UG" undirected graph, or "DAG" directed acyclic graph.
#'
#' @export
#' @examples
#' set.seed(123)
#' adj <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow=3)
#' mat <- simdata(n=100, p=3, B=adj, mu=5, mu_noise=1)
#' res <- PCzinb(mat, method="poi", alpha=0.05)
#' prediction_scores(adj, res, type= "UG")
prediction_scores <- function(trueG ,estimatedG, type){
  if (type == "UG"){
    TP <-  sum(trueG *estimatedG)/2
    FP <-  sum((estimatedG-trueG)==1)/2
    FN <- sum((estimatedG-trueG)==-1)/2
    PPV  <- TP/(TP+FP)
    Se <- TP/(TP+FN)
    F1 <- 2*(PPV*Se)/(PPV+Se)
    return (c(TP=TP,FP=FP,FN=FN,PPV=PPV,Se=Se,F1=F1))
  }else{
    TP <-  sum(trueG *estimatedG)
    FP <-  sum((estimatedG-trueG)==1)
    FN <- sum((estimatedG-trueG)==-1)
    PPV  <- TP/(TP+FP)
    Se <- TP/(TP+FN)
    F1 <- 2*(PPV*Se)/(PPV+Se)
    return (c(TP=TP,FP=FP,FN=FN,PPV=PPV,Se=Se,F1=F1))
  }
}

#' Structure learning with zero-inflated negative binomial model (mean only)
#'
#' This function estimates the adjacency matrix of a ZINB model given a matrix
#' of counts, using the optim function. Uses BiocParallel for parallelization.
#'
#' This approach assumes that the structure of the graph only depends on the
#' mean parameter, treating zero inflation as a technical noise effect. We call
#' this model `zinb0`.
#'
#' @param X the matrix of counts (n times p).
#' @param alpha the significant level of the tests
#' @param maxcard the uper bound of the cardinality of the conditional sets K
#' @param extend if TRUE it considers the union of the tests, otherwise it
#'   considers the intersection.
#' @param nCores number of cores for parallelization
#' @return the estimated adjacency matrix of the graph.
zinb0.noT <- function(X, maxcard, alpha, extend, nCores = 1){
  p <- ncol(X)
  n <- nrow(X)
  
  # Setup BiocParallel backend
  if (nCores > 1) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = nCores)
  } else {
    BPPARAM <- BiocParallel::SerialParam()
  }
  
  ######estimate dispersion parameter
  iter.theta <- 2
  stop.epsilon <- .0001
  zeta <- try(BiocParallel::bplapply(1:p, function(i) {
          iter <- 1
          local.lik <- rep(NA,iter.theta)
          zeta.i <- rep(mean(X[,i])^2/(var(X[,i])-mean(X[,i])),n)
          #2. Estimate parameters of ZINB model with zeta.i given by the first step
          fitadd <- try(fitadd <- optim_fun0noT (beta_mu= rep(1,p), gamma_pi=1, Y=X[,i],
                                                X_mu=X[,-i], zeta.i, n),silent = TRUE)
          if(is(fitadd, "try-error")){
             fit <- glm(X[,i]~X[,-i],family = "poisson")
             fitadd <- optim_fun0noT (beta_mu= fit$coefficients, gamma_pi=1, Y=X[,i],
                                    X_mu=X[,-i], zeta.i, n)
            }
          ####### Calculate loglikelihood at the first iteration with
          #    alpha=fitadd, and C.theta = zeta.i obtained from the above step

            local.lik[1] <- zinb.loglik.regression (alpha=fitadd, Y=X[,i],
                                                   A.mu = cbind(rep(1,n),X_mu=X[, -i]),
                                                    A.pi= matrix(rep(1,n),n,1),
                                                    C.theta = zeta.i)
            for (iter in 2:iter.theta) {
                #1. Estimate zeta with initial value alpha=fitadd given by previuous iteration
                r <- zinb.regression.parseModel (alpha=fitadd, A.mu=cbind(rep(1,n),X[,-i]),
                                                 A.pi= matrix(rep(1,n),n,1))
                zeta.temp <- zinbOptimizeDispersion ( mu=r$logMu, logitPi=r$logitPi,Y=X[,i],n)
               #2. Estimate parameters of ZINB model with zeta given by the first step
                fitadd.temp <- optim_fun0noT (beta_mu=fitadd[1:p],
                                            gamma_pi=fitadd[(p+1)],Y=X[,i],
                                            X_mu=X[,-i], zeta.temp, n)
                local.lik[iter] <- zinb.loglik.regression (alpha=fitadd.temp, Y=X[,i],
                                                           A.mu = cbind(rep(1,n),X_mu=X[,-i]),
                                                           A.pi= matrix(rep(1,n),n,1),
                                                           C.theta = zeta.temp)

                if (local.lik[iter] > local.lik[iter-1]){
                   fitadd <- fitadd.temp
                   zeta.i <- zeta.temp
                  }else break
                if(abs((local.lik[iter]-local.lik[iter-1]) /local.lik[iter-1])< stop.epsilon)
                   break

                iter <- iter+1
                 }

             return(zeta.i)
        }, BPPARAM = BPPARAM), silent = TRUE)

  if (!is.matrix(zeta)){
    zeta <- BiocParallel::bplapply(1:p, function(i) {
      r <- zinb.regression.parseModel (alpha=rep(1,2*p), A.mu=cbind(rep(1,n),scale(X[,-i])),
                                       A.pi= matrix(rep(1,n),n,1))
      zeta.i <- zinbOptimizeDispersion ( mu=r$logMu, logitPi=r$logitPi,Y=X[,i],n)
      return(zeta.i)
    }, BPPARAM = BPPARAM)
    zeta <- do.call(cbind, zeta)
  }
  ############### Estimate adjacency matrix

  adj <- matrix(1,p,p)
  diag(adj) <- 0

  ncard <- 0
  while (ncard <= maxcard) {
    V <- BiocParallel::bplapply(1:p, function(i) {

    neighbor <- which (adj[, i] == 1)
    if (length(neighbor) >= ncard){
      condset <- combn(neighbor, ncard, FUN = list)
      for (j in 1: length(neighbor)){
        condset.temp <- condset
        indcond <- FALSE
        k <- 1
        while (!indcond & k <= length(condset.temp)){
           if (!(neighbor[j] %in% condset.temp[[k]])){

             # initial value
             beta_mu <- c(glm(X[,i]~scale(X[,c(neighbor[j], condset.temp[[k]])])
                              ,family = "poisson")$coefficients)
             gamma_pi <- 0.5

             # fit model with new edges
               fitadd <- optim_fun0noT (beta_mu=beta_mu, gamma_pi, Y=X[,i],
                                    X_mu=scale(X[,c(neighbor[j], condset.temp[[k]])]), zeta[,i], n)
                # calculate loglikelihood of new model
                zinb.loglik.add <- zinb.loglik.regression (alpha=fitadd, Y=X[,i],
                                                          A.mu = cbind(rep(1,n),X_mu=scale(X[,c(neighbor[j], condset.temp[[k]])])),
                                                          A.pi= matrix(rep(1,n),n,1),
                                                          C.theta = zeta[,i])
                # fit model without adding new edges
                if (length(condset.temp[[k]])>0){
                    fitnoadd <- optim_fun0noT (beta_mu=beta_mu[-2], gamma_pi, Y=X[,i],
                                           X_mu=scale(X[, condset.temp[[k]]]), zeta[,i], n)
                    # calculate loglikelihood of model without adding new edges
                    zinb.loglik.noadd <- zinb.loglik.regression (alpha=fitnoadd, Y=X[,i],
                                                                A.mu = cbind(rep(1,n),X_mu=scale(X[, condset.temp[[k]]])),
                                                                A.pi= matrix(rep(1,n),n,1),
                                                                C.theta = zeta[,i])
                     }else{
                          fitnoadd <- optim_fun0noT (beta_mu=beta_mu[c(1,2)], gamma_pi, Y=X[,i],
                                                 X_mu=rep(0,n), zeta[,i], n)
                          # calculate loglikelihood of model without adding new edges
                          zinb.loglik.noadd <- zinb.loglik.regression (alpha=fitadd, Y=X[,i],
                                                                       A.mu = cbind(rep(1,n),rep(0,n)),
                                                                       A.pi= matrix(rep(1,n),n,1),
                                                                       C.theta = zeta[,i])
                        }

                      goodfit.Deviance <- 2*abs(zinb.loglik.noadd -zinb.loglik.add )

                      if (1-pchisq(goodfit.Deviance,1) > alpha){
                          adj[neighbor[j], i] <- 0
                          indcond <- TRUE
                        }
                      }
                      k <- k+1
                    }
                  }
                }
              return(adj[,i])
              }, BPPARAM = BPPARAM)
    
    V <- do.call(cbind, V)
    if (extend){
      adj <- V + t(V)
      adj[which(adj != 0)] <-1
    }else
      adj <- V * t(V)

    ncard <- ncard + 1
  }
  return(adj)
}

#' Structure learning with zero-inflated negative binomial model
#'
#' This function estimates the adjacency matrix of a ZINB model given a matrix
#' of counts, using the optim function. Uses BiocParallel for parallelization.
#'
#' This approach assumes that the structure of the graph depends on both the
#' mean parameter and the zero inflation parameter. We call this model `zinb1`.
#'
#' @param X the matrix of counts (n times p).
#' @param alpha the significant level of the tests
#' @param maxcard the uper bound of the cardinality of the conditional sets K
#' @param extend if TRUE it considers the union of the tests, otherwise it
#'   considers the intersection.
#' @param nCores number of cores for parallelization
#' @return the estimated adjacency matrix of the graph.
zinb1.noT <- function(X, maxcard, alpha, extend, nCores = 1){
  p <- ncol(X)
  n <- nrow(X)
  
  # Setup BiocParallel backend
  if (nCores > 1) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = nCores)
  } else {
    BPPARAM <- BiocParallel::SerialParam()
  }
  
  ###### Estimate dispersion parameter
  iter.theta <- 2
  stop.epsilon <- .0001
  zeta <- try(BiocParallel::bplapply(1:p, function(i) {
          iter <- 1
          local.lik <- rep(NA,iter.theta)
          #1. Estimate zeta
          zeta.i <- rep(mean(X[,i])^2/(var(X[,i])-mean(X[,i])),n)
          #2. Estimate parameters of ZINB model with zeta.i given by the first step
          fitadd <- try(fitadd <- optim_funnoT (beta_mu= rep(1,p), gamma_pi=rep(1,p), Y=X[,i],
                                                X_mu=X[,-i], zeta.i, n),silent = TRUE)
          if(is(fitadd, "try-error")){
             fit <- glm(X[,i]~X[,-i],family = "poisson")
             fitadd <- optim_funnoT (beta_mu= fit$coefficients, gamma_pi=rep(1,p), Y=X[,i],
                                    X_mu=X[,-i], zeta.i, n)
            }
          ####### Calculate loglikelihood at the first iteration with
          #    alpha=fitadd, and C.theta = zeta.i obtained from the above step

            local.lik[1] <- zinb.loglik.regression (alpha=fitadd, Y=X[,i],
                                                   A.mu = cbind(rep(1,n),X_mu=X[, -i]),
                                                    A.pi = cbind(rep(1,n),X_mu=X[, -i]),
                                                    C.theta = zeta.i)
            for (iter in 2:iter.theta) {
                #1. Estimate zeta with initial value alpha=fitadd given by previuous iteration
                r <- zinb.regression.parseModel (alpha=fitadd, A.mu=cbind(rep(1,n),X[,-i]),
                                                 A.pi=cbind(rep(1,n),X_mu=X[, -i]))
                zeta.temp <- zinbOptimizeDispersion ( mu=r$logMu, logitPi=r$logitPi,Y=X[,i],n)
               #2. Estimate parameters of ZINB model with zeta given by the first step
                fitadd.temp <- optim_funnoT (beta_mu=fitadd[1:p],
                                            gamma_pi=fitadd[(p+1):(2*p)],Y=X[,i],
                                            X_mu=X[,-i], zeta.temp, n)
                local.lik[iter] <- zinb.loglik.regression (alpha=fitadd.temp, Y=X[,i],
                                                           A.mu = cbind(rep(1,n),X_mu=X[,-i]),
                                                           A.pi = cbind(rep(1,n),X_mu=X[, -i]),
                                                           C.theta = zeta.temp)

                if (local.lik[iter] > local.lik[iter-1]){
                   fitadd <- fitadd.temp
                   zeta.i <- zeta.temp
                  }else break
                if(abs((local.lik[iter]-local.lik[iter-1]) /local.lik[iter-1])< stop.epsilon)
                   break

                iter <- iter+1
                 }

             return(zeta.i)
        }, BPPARAM = BPPARAM), silent = TRUE)
  
  if (!is.matrix(zeta)){
    zeta <- BiocParallel::bplapply(1:p, function(i) {
      r <- zinb.regression.parseModel (alpha=rep(1,2*p), A.mu=cbind(rep(1,n),scale(X[,-i])),
                                          A.pi=cbind(rep(1,n),scale(X[,-i])))
      zeta.i <- zinbOptimizeDispersion ( mu=r$logMu, logitPi=r$logitPi,Y=X[,i],n)
      return(zeta.i)
    }, BPPARAM = BPPARAM)
    zeta <- do.call(cbind, zeta)
  }

  ############### Estimate adjacency matrix

  adj <- matrix(1,p,p)
  diag(adj) <- 0
  ncard <- 0
  while (ncard <= maxcard) {
    V <- BiocParallel::bplapply(1:p, function(i) {

    neighbor <- which (adj[, i] == 1)
    if (length(neighbor) >= ncard){
      condset <- combn(neighbor, ncard, FUN = list)
      for (j in 1: length(neighbor)){
        condset.temp <- condset
        indcond <- FALSE
        k <- 1
        while (!indcond & k <= length(condset.temp)){
           if (!(neighbor[j] %in% condset.temp[[k]])){

             # initial value
             beta_mu <- c(glm(X[,i]~scale(X[,c(neighbor[j], condset.temp[[k]])])
                              ,family = "poisson")$coefficients)
             gamma_pi <- rep(0.5,2+ncard)

             # fit model with new edges
               fitadd <- optim_funnoT (beta_mu=beta_mu, gamma_pi, Y=X[,i],
                                    X_mu=scale(X[,c(neighbor[j], condset.temp[[k]])]), zeta[,i], n)
                # calculate loglikelihood of new model
                zinb.loglik.add <- zinb.loglik.regression (alpha=fitadd, Y=X[,i],
                                                          A.mu = cbind(rep(1,n),X_mu=scale(X[,c(neighbor[j], condset.temp[[k]])])),
                                                          A.pi = cbind(rep(1,n),X_mu=scale(X[,c(neighbor[j], condset.temp[[k]])])),
                                                          C.theta = zeta[,i])
                # fit model without adding new edges
                if (length(condset.temp[[k]])>0){
                    fitnoadd <- optim_funnoT (beta_mu=beta_mu[-2], gamma_pi[-1], Y=X[,i],
                                           X_mu=scale(X[, condset.temp[[k]]]), zeta[,i], n)
                    # calculate loglikelihood of model without adding new edges
                    zinb.loglik.noadd <- zinb.loglik.regression (alpha=fitnoadd, Y=X[,i],
                                                                A.mu = cbind(rep(1,n),X_mu=scale(X[, condset.temp[[k]]])),
                                                                A.pi = cbind(rep(1,n),X_mu=scale(X[, condset.temp[[k]]])),
                                                                C.theta = zeta[,i])
                     }else{
                          fitnoadd <- optim_funnoT (beta_mu=beta_mu[c(1,2)], gamma_pi, Y=X[,i],
                                                 X_mu=rep(0,n), zeta[,i], n)
                          # calculate loglikelihood of model without adding new edges
                          zinb.loglik.noadd <- zinb.loglik.regression (alpha=fitadd, Y=X[,i],
                                                                       A.mu = cbind(rep(1,n),rep(0,n)),
                                                                       A.pi = cbind(rep(1,n),rep(0,n)),
                                                                       C.theta = zeta[,i])
                        }
                        ### Deviance tests
                      goodfit.Deviance <- 2*(zinb.loglik.add-zinb.loglik.noadd  )

                      if (pchisq(goodfit.Deviance,2,lower.tail=FALSE) > alpha){
                          adj[neighbor[j], i] <- 0
                          indcond <- TRUE
                        }
                      }
                      k <- k+1
                    }
                  }
                }
              return(adj[,i])
              }, BPPARAM = BPPARAM)
    
    V <- do.call(cbind, V)
    if (extend){
      adj <- V + t(V)
      adj[which(adj != 0)] <-1
    }else
      adj <- V * t(V)

    ncard <- ncard + 1
  }
  return(adj)
}

utils::globalVariables("i")
