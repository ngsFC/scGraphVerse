# =========================
# JRF Implementation - EXACT copy of working JRF package functions
# =========================

.jrf_onetarget <- function(
    x, y = NULL, xtest = NULL, ytest = NULL, ntree, sampsize,
    totsize = if (replace) ncol(x) else ceiling(0.632 * ncol(x)),
    mtry = if (!is.null(y) && !is.factor(y)) {
        max(
            floor(nrow(x) / 3),
            1
        )
    } else {
        floor(sqrt(nrow(x)))
    }, replace = TRUE, classwt = NULL,
    cutoff, strata, nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1,
    maxnodes = NULL, importance = FALSE, localImp = FALSE, nPerm = 1,
    proximity, oob.prox = proximity, norm.votes = TRUE, do.trace = FALSE,
    keep.forest = !is.null(y) && is.null(xtest), corr.bias = FALSE,
    keep.inbag = FALSE, nclasses, ...) {
    ww <- 1 / sampsize # EXACT calculation from working version
    nclass <- mylevels <- ipi <- sw <- NULL
    addclass <- is.null(y)
    classRF <- addclass || is.factor(y)
    if (!classRF && length(unique(y)) <= 5) {
        warning(
            "The response has five or fewer unique values. ",
            "Are you sure you want to do regression?"
        )
    }
    if (classRF && !addclass && length(unique(y)) < 2) {
        stop("Need at least two classes to do classification.")
    }
    n <- ncol(y) 
    p <- nrow(x) / nclasses
    if (n == 0) {
        stop("data (x) has 0 rows")
    }
    x.row.names <- rownames(x)
    x.col.names <- if (is.null(colnames(x))) {
        seq_len(ncol(x))
    } else {
        colnames(x)
    }
    keep.forest <- !is.null(y)
    xtest <- NULL
    ytest <- NULL
    testdat <- !is.null(xtest)
    if (testdat) {
        if (ncol(x) != ncol(xtest)) {
            stop("x and xtest must have same number of columns")
        }
        ntest <- nrow(xtest)
        xts.row.names <- rownames(xtest)
    }
    prox <- proxts <- double(1)
    if (any(is.na(x))) {
        stop("NA not permitted in predictors")
    }
    if (testdat && any(is.na(xtest))) {
        stop("NA not permitted in xtest")
    }
    if (any(is.na(y))) {
        stop("NA not permitted in response")
    }
    if (!is.null(ytest) && any(is.na(ytest))) {
        stop("NA not permitted in ytest")
    }
    if (is.data.frame(x)) {
        xlevels <- lapply(x, function(z) if (is.factor(z)) levels(z) else NULL)
        ncat <- vapply(xlevels, length, integer(1))
        ncat <- ifelse(vapply(x, is.ordered, logical(1)), 1, ncat)
        x <- data.matrix(x)
        if (testdat) {
            if (!is.data.frame(xtest)) {
                stop("xtest must be data frame if x is")
            }
            xfactor <- which(vapply(xtest, is.factor, logical(1)))
            if (length(xfactor) > 0) {
                for (i in xfactor) {
                    if (any(!levels(xtest[[i]]) %in% xlevels[[i]])) {
                        stop("New factor levels in xtest not present in x")
                    }
                    xtest[[i]] <- factor(xlevels[[i]][match(
                        xtest[[i]],
                        xlevels[[i]]
                    )], levels = xlevels[[i]])
                }
            }
            xtest <- data.matrix(xtest)
        }
    } else {
        ncat <- rep(1, p)
        xlevels <- as.list(rep(0, p))
    }
    maxcat <- max(ncat)
    if (maxcat > 32) {
        stop(
            "Can not handle categorical predictors with more ",
            "than 32 categories."
        )
    }
    addclass <- FALSE
    proximity <- addclass
    impout <- matrix(0, p * nclasses, 2)
    impSD <- matrix(0, p * nclasses, 1)
    nsample <- if (addclass) {
        2 * n
    } else {
        n
    }
    Stratify <- length(n) > 1
    nodesize <- 5
    nrnodes <- 2 * trunc(n / max(1, nodesize - 4)) + 1
    maxnodes <- NULL
    if (!is.null(maxnodes)) {
        maxnodes <- 2 * maxnodes - 1
        if (maxnodes > nrnodes) {
            warning("maxnodes exceeds its max value.")
        }
        nrnodes <- min(c(nrnodes, max(c(maxnodes, 1))))
    }
    storage.mode(x) <- "double"
    xtest <- double(1)
    ytest <- double(1)
    ntest <- 1
    labelts <- FALSE
    nt <- if (keep.forest) {
        ntree
    } else {
        1
    }
    nPerm <- 1
    do.trace <- FALSE
    oob.prox <- FALSE
    corr.bias <- FALSE
    keep.inbag <- FALSE
    impmat <- double(1)
    replace <- TRUE

    rfout <- .C("regRF", x, y, ww, as.integer(c(
        totsize,
        p
    )),
    sampsize = as.integer(sampsize), as.integer(totsize),
    as.integer(nodesize), as.integer(nrnodes), as.integer(ntree),
    as.integer(mtry), as.integer(c(
        importance, localImp,
        nPerm
    )), as.integer(ncat), as.integer(maxcat),
    as.integer(do.trace), as.integer(proximity), as.integer(oob.prox),
    as.integer(corr.bias), ypred = double(n * nclasses),
    impout = impout, impmat = impmat, impSD = impSD,
    prox = prox, ndbigtree = integer(ntree),nodestatus=matrix(integer(nrnodes *
        nt * nclasses), ncol = nt), leftDaughter = matrix(integer(nrnodes *
        nt * nclasses), ncol = nt), rightDaughter = matrix(integer(nrnodes *
        nt * nclasses), ncol = nt), nodepred = matrix(double(nrnodes *
        nt * nclasses), ncol = nt), bestvar = matrix(integer(nrnodes *
        nt * nclasses), ncol = nt), xbestsplit = matrix(double(nrnodes *
        nt * nclasses), ncol = nt), mse = double(ntree *
        nclasses), keep = as.integer(c(keep.forest, keep.inbag)),
    replace = as.integer(replace), testdat = as.integer(testdat),
    xts = xtest, ntest = as.integer(ntest), yts = as.double(ytest),
    labelts = as.integer(labelts), ytestpred = double(ntest),
    proxts = proxts, msets = double(if (labelts) ntree else 1),
    coef = double(2), oob.times = integer(n),
    inbag = if (keep.inbag) matrix(integer(n * ntree), n) else integer(1),
    as.integer(nclasses), PACKAGE = "scGraphVerse"
    )[c(
        16:28,
        36:41
    )]

    if (keep.forest) {
        max.nodes <- max(rfout$ndbigtree)
        rfout$nodestatus <- rfout$nodestatus[seq_len(max.nodes), ,
            drop = FALSE
        ]
        rfout$bestvar <- rfout$bestvar[seq_len(max.nodes), , drop = FALSE]
        rfout$nodepred <- rfout$nodepred[seq_len(max.nodes), , drop = FALSE]
        rfout$xbestsplit <- rfout$xbestsplit[seq_len(max.nodes), ,
            drop = FALSE
        ]
        rfout$leftDaughter <- rfout$leftDaughter[seq_len(max.nodes), ,
            drop = FALSE
        ]
        rfout$rightDaughter <- rfout$rightDaughter[seq_len(max.nodes), ,
            drop = FALSE
        ]
    }
    cl <- match.call()
    cl[[1]] <- as.name("randomForest")
    ypred <- rfout$ypred
    if (any(rfout$oob.times < 1)) {
        ypred[rfout$oob.times == 0] <- NA
    }
    out <- list(
        call = cl, type = "regression", predicted = 0,
        mse = rfout$mse, rsq = 1 - rfout$mse / (var(y[1, ]) *
            (n - 1) / n),oob.times=rfout$oob.times,importance=if (importance) {
            matrix(
                rfout$impout,
                p * nclasses, 2
            )
        } else {
            matrix(rfout$impout, ncol = 1)
        },
        importanceSD = if (importance) rfout$impSD else NULL,
        localImportance = if (localImp) {
            matrix(rfout$impmat,
                p, n,
                dimnames = list(x.col.names, x.row.names)
            )
        } else {
            NULL
        },
        proximity = if (proximity) {
            matrix(rfout$prox, n,
                n,
                dimnames = list(x.row.names, x.row.names)
            )
        } else {
            NULL
        },
        ntree = ntree, mtry = mtry,
        forest = if (keep.forest) {
            c(
                rfout[c(
                    "ndbigtree",
                    "nodestatus", "leftDaughter", "rightDaughter",
                    "nodepred", "bestvar", "xbestsplit"
                )], list(ncat = ncat),
                list(nrnodes = max.nodes), list(ntree = ntree),
                list(xlevels = xlevels)
            )
        } else {
            NULL
        },
        coefs = if (corr.bias) rfout$coef else NULL,
        y = y, test = if (testdat) {
            list(
                predicted = structure(rfout$ytestpred, names = x.row.names),
                mse = if (labelts) rfout$msets else NULL, rsq = if (labelts) {
                    1 -
                        rfout$msets / (var(ytest) * (n - 1) / n)
                } else {
                    NULL
                },
                proximity = if (proximity) {
                    matrix(rfout$proxts / ntree,
                        nrow = ntest, dimnames = list(
                            x.row.names,
                            c(x.row.names, x.row.names)
                        )
                    )
                } else {
                    NULL
                }
            )
        } else {
            NULL
        }, inbag = if (keep.inbag) {
            matrix(rfout$inbag,
                nrow(rfout$inbag),
                dimnames = list(
                    x.row.names,
                    NULL
                )
            )
        } else {
            NULL
        }
    )
    class(out) <- "randomForest"
    return(out)
}

.importance <- function(x, scale = TRUE) {
    type <- NULL
    class <- NULL
    if (!inherits(x, "randomForest")) {
        stop("x is not of class randomForest")
    }
    classRF <- x$type != "regression"
    hasImp <- !is.null(dim(x$importance)) || ncol(x$importance) ==
        1
    hasType <- !is.null(type)
    if (hasType && type == 1 && !hasImp) {
        stop("That measure has not been computed")
    }
    allImp <- is.null(type) && hasImp
    if (hasType) {
        if (!(type %in% seq_len(2))) {
            stop("Wrong type specified")
        }
        if (type == 2 && !is.null(class)) {
            stop("No class-specific measure for that type")
        }
    }
    imp <- x$importance
    if (hasType && type == 2) {
        if (hasImp) {
            imp <- imp[, ncol(imp), drop = FALSE]
        }
    } else {
        if (scale) {
            SD <- x$importanceSD
            imp[, -ncol(imp)] <- imp[, -ncol(imp), drop = FALSE] / ifelse(SD <
                .Machine$double.eps, 1, SD)
        }
        if (!allImp) {
            if (is.null(class)) {
                imp <- imp[, ncol(imp) - 1, drop = FALSE]
            } else {
                whichCol <- if (classRF) {
                    match(class, colnames(imp))
                } else {
                    1
                }
                if (is.na(whichCol)) {
                    stop("Class ", class, " not found.")
                }
                imp <- imp[, whichCol, drop = FALSE]
            }
        }
    }
    imp <- imp[, 2]
    imp
}

.jrf_network <- function(data_list, ntree = 1000, mtry = NULL) {
    X <- data_list
    nclasses <- length(X)
    sampsize <- rep(0, nclasses)

    for (j in seq_len(nclasses)) sampsize[j] <- dim(X[[j]])[2]
    tot <- max(sampsize) 
    p <- dim(X[[1]])[1]

    genes.name <- rownames(X[[1]])
    if (is.null(genes.name)) genes.name <- paste0("G", seq_len(p))

    # Set default mtry if not provided
    if (is.null(mtry)) mtry <- round(sqrt(p - 1))

    imp <- array(0, c(p, length(genes.name), nclasses))
    imp.final <- matrix(0, p * (p - 1) / 2, nclasses)
    vec1 <- matrix(rep(genes.name, p), p, p)
    vec2 <- t(vec1)
    vec1 <- vec1[lower.tri(vec1, diag = FALSE)]
    vec2 <- vec2[lower.tri(vec2, diag = FALSE)]
    index <- seq(1, p)

    for (j in seq_along(genes.name)) {
        covar <- matrix(0, (p - 1) * nclasses, tot)
        y <- matrix(0, nclasses, tot)

        for (c in seq_len(nclasses)) {
            y[c, seq(1, sampsize[c])] <- as.matrix(X[[c]][j, ])
            covar[
                seq((c - 1) * (p - 1) + 1, c * (p - 1)),
                seq(1, sampsize[c])
            ] <- X[[c]][-j, ]
        }

        jrf.out <- .jrf_onetarget(
            x = covar, y = y, mtry = mtry,
            importance = TRUE, sampsize = sampsize,
            nclasses = nclasses, ntree = ntree
        )

        for (s in seq_len(nclasses)) {
            imp[-j, j, s] <- .importance(jrf.out, scale = FALSE)[
                seq((p - 1) * (s - 1) + 1, (p - 1) * (s - 1) + p - 1)
            ]
        }
    }

    for (s in seq_len(nclasses)) {
        imp.s <- imp[, , s]
        t.imp <- t(imp.s)
        imp.final[, s] <- (imp.s[lower.tri(imp.s, diag = FALSE)] +
            t.imp[lower.tri(t.imp, diag = FALSE)]) / 2
    }

    # Convert to scGraphVerse format
    result_list <- vector("list", nclasses)
    for (s in seq_len(nclasses)) {
        result_list[[s]] <- data.frame(
            gene1 = as.character(vec1),
            gene2 = as.character(vec2),
            importance = imp.final[, s],
            stringsAsFactors = FALSE
        )
    }

    return(result_list)
}
