# =========================
# Internal helper for JRF
# =========================

# 1) JRF on a single target (C backend)
.jrf_onetarget <- function(x, y, ntree, mtry, sampsize, nclasses, importance = TRUE) {
  totsize <- sum(sampsize)  # total observations
  nrnodes <- 2 * trunc(ncol(x) / max(1, 5 - 4)) + 1  # based on columns (samples)
  # For numeric-only data
  ncat <- rep(1L, nrow(x))       # integer vector
  maxcat <- 1L                   # single integer
  keep <- as.integer(c(1, 0))    # keep flags


  rfout <- .C("regRF",
    # ---- Input arguments ----
    x,                                        # 1 double matrix (flattened)
    as.vector(y),                             # 2 double vector (flattened)
    as.double(1 / sampsize),                  # 3 weights
    as.integer(c(nrow(x), totsize)),          # 4 xdim: (variables, samples)
    as.integer(sampsize),                     # 5 sampsize per class
    as.integer(totsize),                      # 6 total sample size
    as.integer(5),                            # 7 nodesize (fixed to 5)
    as.integer(nrnodes),                      # 8 nrnodes
    as.integer(ntree),                        # 9 ntree
    as.integer(mtry),                         #10 mtry
    as.integer(c(importance, FALSE, 1)),      #11 importance flags
    as.integer(ncat),                         #12 ncat
    as.integer(maxcat),                       #13 maxcat
    as.integer(0),                            #14 do.trace
    as.integer(0),                            #15 proximity
    as.integer(0),                            #16 oob.prox
    as.integer(0),                            #17 corr.bias
    # ---- Output placeholders ----
    ypred = double(totsize * nclasses),       #18 predictions
    impout = double(nrow(x) * 2),             #19 variable importance
    impmat = double(1),                       #20
    impSD = double(nrow(x)),                  #21
    prox = double(1),                         #22
    ndbigtree = integer(ntree),               #23
    nodestatus = integer(nrnodes * ntree * nclasses), #24
    leftDaughter = integer(1),                #25
    rightDaughter = integer(1),               #26
    nodepred = double(1),                     #27
    bestvar = integer(1),                     #28
    xbestsplit = double(1),                   #29
    mse = double(ntree * nclasses),           #30
    keep = as.integer(c(1, 0)),               #31 keep forest, keep inbag
    replace = as.integer(1),                  #32 bootstrap
    testdat = as.integer(0),                  #33
    xts = double(1),                          #34
    ntest = as.integer(1),                    #35
    yts = as.double(0),                       #36
    labelts = as.integer(0),                  #37
    ytestpred = double(1),                    #38
    proxts = double(1),                       #39
    msets = double(1),                        #40
    coef = double(2),                         #41
    oob.times = integer(totsize),             #42
    inbag = integer(1),                       #43
    as.integer(nclasses)                      #44
  )

  # Return minimal object (importance matrix)
  out <- list(
    importance = matrix(rfout$impout, ncol = 2)
  )
  class(out) <- "randomForest"
  return(out)
}

# 2) Main JRF network inference
.jrf_network <- function(data_list, ntree = 1000, mtry = NULL) {
  nclasses <- length(data_list)
  sampsize <- vapply(data_list, ncol, FUN.VALUE = integer(1))
  p <- nrow(data_list[[1]])
  if (is.null(mtry)) mtry <- max(floor(sqrt(p - 1)), 1)

  genes.name <- rownames(data_list[[1]])
  if (is.null(genes.name)) genes.name <- paste0("G", seq_len(p))

  # Storage for importance
  imp <- array(0, dim = c(p, length(genes.name), nclasses))

  # For each target gene
  for (j in seq_len(p)) {
    # Build covariate matrix: (p-1) * nclasses rows, totsize cols
    totsize <- sum(sampsize)
    covar <- matrix(0, (p - 1) * nclasses, totsize)
    y <- matrix(0, nclasses, totsize)

    for (c in seq_len(nclasses)) {
      idx_start <- if (c == 1) 1 else sum(sampsize[1:(c - 1)]) + 1
      idx_end <- sum(sampsize[1:c])
      y[c, idx_start:idx_end] <- data_list[[c]][j, ]
      covar[((c - 1) * (p - 1) + 1):(c * (p - 1)), idx_start:idx_end] <- data_list[[c]][-j, ]
    }

    # Run regRF on combined data
    jrf.out <- .jrf_onetarget(x = covar, y = y, ntree = ntree, mtry = mtry,
                               sampsize = sampsize, nclasses = nclasses)

    # Extract importance scores per class
    for (s in seq_len(nclasses)) {
      imp[-j, j, s] <- jrf.out$importance[seq((p - 1) * (s - 1) + 1, (p - 1) * s), 2]
    }
  }

  # Average importance for symmetric matrix
  imp.final <- matrix(0, nrow = p * (p - 1) / 2, ncol = nclasses)
  vec1 <- matrix(rep(genes.name, p), p, p)
  vec2 <- t(vec1)
  vec1 <- vec1[lower.tri(vec1, diag = FALSE)]
  vec2 <- vec2[lower.tri(vec2, diag = FALSE)]

  for (s in seq_len(nclasses)) {
    imp.s <- imp[, , s]
    t.imp <- t(imp.s)
    imp.final[, s] <- (imp.s[lower.tri(imp.s, diag = FALSE)] +
                         t.imp[lower.tri(t.imp, diag = FALSE)]) / 2
  }

  out <- cbind(as.character(vec1), as.character(vec2), as.data.frame(imp.final), stringsAsFactors = FALSE)
  colnames(out) <- c("gene1", "gene2", paste0("importance", seq_len(nclasses)))
  return(out)
}
