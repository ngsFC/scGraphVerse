#' Internal JRF implementation using native C functions
#' Adapted from JRF 0.1-4 (CRAN)
#' @keywords internal
NULL

# Extract variable importance (from original JRF)
.jrf_importance <- function(x, scale = TRUE) {
  type <- NULL
  class <- NULL
  if (!inherits(x, "randomForest"))
    stop("x is not of class randomForest")

  imp <- x$importance
  if (scale) {
    SD <- x$importanceSD
    imp[, -ncol(imp)] <- imp[, -ncol(imp), drop = FALSE] /
      ifelse(SD < .Machine$double.eps, 1, SD)
  }
  imp <- imp[, 2]
  imp
}

# JRF_onetarget (C backend)
.jrf_onetarget <- function(x, y, ntree, mtry, sampsize, nclasses, importance = TRUE) {
  # Total number of samples across all classes
  totsize <- sum(sampsize)
  nrnodes <- 2 * trunc(totsize / max(1, 5 - 4)) + 1
  ncat <- rep(1L, nrow(x))   # all predictors numeric
  maxcat <- 1L
  # This calls the native compiled code via .C like original JRF
  rfout <- .C("regRF",
    # Arguments in correct order:
    x,                                        # 1
    y,                                        # 2
    as.double(1 / sampsize),                  # 3 ww (weights)
    as.integer(c(totsize, nrow(x))),          # 4 xdim: (samples, predictors)
    as.integer(totsize),                      # 5 nsample
    as.integer(5),                            # 6 nodesize (fixed to 5 like original)
    as.integer(nrnodes),                      # 7 nrnodes
    as.integer(ntree),                        # 8 ntree
    as.integer(mtry),                         # 9 mtry
    as.integer(c(importance, FALSE, 1)),      #10 imp: (importance, localImp, nPerm)
    as.integer(ncat),                         #11 ncat (vector)
    as.integer(maxcat),                       #12 maxcat
    as.integer(0),                            #13 do.trace
    as.integer(0),                            #14 proximity
    as.integer(0),                            #15 oob.prox
    as.integer(0),                            #16 corr.bias
    ypred = double(totsize * nclasses),       #17 ypred
    impout = double(nrow(x) * 2),             #18 impout
    impmat = double(1),                       #19 impmat
    impSD = double(nrow(x)),                  #20 impSD
    prox = double(1),                         #21 prox
    ndbigtree = integer(ntree),               #22 ndbigtree
    nodestatus = integer(nrnodes * ntree * nclasses), #23 nodestatus
    leftDaughter = integer(1),                #24 leftDaughter
    rightDaughter = integer(1),               #25 rightDaughter
    nodepred = double(1),                     #26 nodepred
    bestvar = integer(1),                     #27 bestvar
    xbestsplit = double(1),                   #28 xbestsplit
    mse = double(ntree * nclasses),           #29 mse
    keep = as.integer(c(1, 0)),               #30 keep (keep.forest=1, keep.inbag=0)
    replace = as.integer(1),                  #31 replace (bootstrap)
    testdat = as.integer(0),                  #32 testdat
    xts = double(1),                          #33 xts
    ntest = as.integer(1),                    #34 ntest
    yts = as.double(0),                       #35 yts
    labelts = as.integer(0),                  #36 labelts
    ytestpred = double(1),                    #37 ytestpred
    proxts = double(1),                       #38 proxts
    msets = double(1),                        #39 msets
    coef = double(2),                         #40 coef
    oob.times = integer(totsize),             #41 oob.times
    inbag = integer(1),                       #42 inbag
    as.integer(nclasses)                      #43 nclasses
  )

  # Build an object similar to randomForest
  out <- list(
    importance = matrix(rfout$impout, ncol = 2)
  )
  class(out) <- "randomForest"
  return(out)
}

# Main JRF function equivalent

JRF_internal <- function(X, ntree = 1000, mtry = NULL, genes.name = NULL) {
  nclasses <- length(X)
  if (is.null(genes.name)) genes.name <- rownames(X[[1]])
  if (is.null(mtry)) mtry <- round(sqrt(length(genes.name) - 1))

  p <- length(genes.name)
  sampsize <- vapply(X, ncol, numeric(1))
  tot <- max(sampsize)

  imp <- array(0, c(p, p, nclasses))
  vec1 <- matrix(rep(genes.name, p), p, p)
  vec2 <- t(vec1)
  vec1 <- vec1[lower.tri(vec1, diag = FALSE)]
  vec2 <- vec2[lower.tri(vec2, diag = FALSE)]

  for (j in seq_len(p)) {
    covar <- matrix(0, (p - 1) * nclasses, tot)
    y <- matrix(0, nclasses, tot)

    for (c in seq_len(nclasses)) {
      y[c, seq_len(sampsize[c])] <- as.matrix(X[[c]][j, ])
      covar[((c - 1) * (p - 1) + 1):(c * (p - 1)), seq_len(sampsize[c])] <- X[[c]][-j, ]
    }

    jrf.out <- .jrf_onetarget(
      x = covar, y = y,
      ntree = ntree, mtry = mtry,
      sampsize = sampsize, nclasses = nclasses,
      importance = TRUE
    )

    for (s in seq_len(nclasses)) {
      imp[-j, j, s] <- .jrf_importance(jrf.out)[seq((p - 1) * (s - 1) + 1, (p - 1) * (s - 1) + p - 1)]
    }
  }

  imp.final <- matrix(0, p * (p - 1) / 2, nclasses)
  for (s in seq_len(nclasses)) {
    imp.s <- imp[, , s]
    t.imp <- t(imp.s)
    imp.final[, s] <- (imp.s[lower.tri(imp.s, diag = FALSE)] + t.imp[lower.tri(t.imp, diag = FALSE)]) / 2
  }

  out <- cbind(vec1, vec2, as.data.frame(imp.final), stringsAsFactors = FALSE)
  colnames(out) <- c("gene1", "gene2", paste0("importance", seq_len(nclasses)))
  return(out)
}

