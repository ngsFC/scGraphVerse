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
  # This calls the native compiled code via .C like original JRF
  rfout <- .C("regRF",
    x,
    y,
    as.double(1 / sampsize),           # weights
    as.integer(c(sampsize[1], nrow(x))),
    as.integer(sum(sampsize)),
    as.integer(5),                     # nodesize
    as.integer(2 * trunc(sum(sampsize) / max(1, 5 - 4)) + 1), # nrnodes
    as.integer(ntree),
    as.integer(mtry),
    as.integer(c(importance, FALSE, 1)),
    as.integer(rep(1, nrow(x))),
    as.integer(1),
    as.integer(0), # do.trace
    as.integer(0), # proximity
    as.integer(0), # oob.prox
    as.integer(0), # corr.bias
    ypred = double(sum(sampsize) * nclasses),
    impout = double(nrow(x) * 2),
    impmat = double(1),
    impSD = double(nrow(x)),
    prox = double(1),
    ndbigtree = integer(ntree),
    nodestatus = integer((2 * trunc(sum(sampsize) / max(1, 5 - 4)) + 1) * ntree * nclasses),
    leftDaughter = integer(1),
    rightDaughter = integer(1),
    nodepred = double(1),
    bestvar = integer(1),
    xbestsplit = double(1),
    mse = double(ntree * nclasses),
    keep = as.integer(c(1, 0)),
    replace = as.integer(1),
    testdat = as.integer(0),
    xts = double(1),
    ntest = as.integer(1),
    yts = as.double(0),
    labelts = as.integer(0),
    ytestpred = double(1),
    proxts = double(1),
    msets = double(1),
    coef = double(2),
    oob.times = integer(sum(sampsize)),
    inbag = integer(1),
    as.integer(nclasses)
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

