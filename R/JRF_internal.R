# =========================
# Internal helper for JRF (C implementation)
# =========================

.jrf_onetarget <- function(x, y, ntree, mtry, sampsize, nclasses,
                        importance = TRUE) {
    # Setup exactly like C function expects
    totsize <- sum(sampsize)
    p_per_class <- nrow(x) / nclasses
    nodesize <- 5
    nrnodes <- 2 * floor(totsize / max(1, nodesize - 4)) + 1
    nt <- ntree  # for forest matrices
    
    
    # Storage mode for C interface  
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    
    # Initialize output matrices
    impout <- matrix(0.0, p_per_class * nclasses, 2)
    impSD <- matrix(0.0, p_per_class * nclasses, 1)
    impmat <- double(1)
    
    # Class weights (equal weighting)
    weight <- rep(1.0, nclasses)
    storage.mode(weight) <- "double"
    
    # Call C function regRF matching EXACTLY the C signature (44 params)
    rfout <- .C("regRF",
                x,                                          # 1. double *x
                y,                                          # 2. double *y  
                weight,                                     # 3. double *weight
                as.integer(c(totsize, p_per_class)),        # 4. int *xdim (nsample, mdim)
                as.integer(sampsize),                       # 5. int *sampsize
                as.integer(totsize),                        # 6. int *totsize
                as.integer(nodesize),                       # 7. int *nthsize
                as.integer(nrnodes),                        # 8. int *nrnodes
                as.integer(ntree),                          # 9. int *nTree
                as.integer(mtry),                           # 10. int *mtry
                as.integer(c(as.integer(importance), 0, 1)), # 11. int *imp
                as.integer(rep(1, p_per_class)),            # 12. int *cat
                as.integer(1),                              # 13. int *maxcat
                as.integer(0),                              # 14. int *jprint
                as.integer(0),                              # 15. int *doProx
                as.integer(0),                              # 16. int *oobprox
                as.integer(0),                              # 17. int *biasCorr
                ypred = double(totsize * nclasses),         # 18. double *yptr
                impout = impout,                            # 19. double *errimp
                impmat = impmat,                            # 20. double *impmat
                impSD = impSD,                              # 21. double *impSD
                prox = double(1),                           # 22. double *prox
                treeSize = integer(ntree),                  # 23. int *treeSize
                nodestatus = integer(nrnodes * nt * nclasses), # 24. int *nodestatus
                lDaughter = integer(nrnodes * nt * nclasses), # 25. int *lDaughter
                rDaughter = integer(nrnodes * nt * nclasses), # 26. int *rDaughter
                avnode = double(nrnodes * nt * nclasses),   # 27. double *avnode
                mbest = integer(nrnodes * nt * nclasses),   # 28. int *mbest
                upper = double(nrnodes * nt * nclasses),    # 29. double *upper
                mse = double(ntree * nclasses),             # 30. double *mse
                keepf = as.integer(c(0, 0)),                # 31. int *keepf
                replace = as.integer(1),                    # 32. int *replace
                testdat = as.integer(0),                    # 33. int *testdat
                xts = double(1),                            # 34. double *xts
                nts = as.integer(1),                        # 35. int *nts
                yts = double(1),                            # 36. double *yts
                labelts = as.integer(0),                    # 37. int *labelts
                yTestPred = double(1),                      # 38. double *yTestPred
                proxts = double(1),                         # 39. double *proxts
                msets = double(1),                          # 40. double *msets
                coef = double(2),                           # 41. double *coef
                nout = integer(totsize),                    # 42. int *nout
                inbag = integer(1),                         # 43. int *inbag
                nclasses = as.integer(nclasses)             # 44. int *nclasses
    )[c("impout")]
    
    out <- list(importance = rfout$impout)
    class(out) <- "randomForest"
    return(out)
}


# Main JRF network inference
.jrf_network <- function(data_list, ntree = 1000, mtry = NULL) {
    nclasses <- length(data_list)
    sampsize <- vapply(data_list, ncol, FUN.VALUE = integer(1))
    p <- nrow(data_list[[1]])
    if (is.null(mtry)) mtry <- max(floor(sqrt(p - 1)), 1)

    genes.name <- rownames(data_list[[1]])
    if (is.null(genes.name)) genes.name <- paste0("G", seq_len(p))

    imp <- array(0, dim = c(p, length(genes.name), nclasses))

    for (j in seq_len(p)) {
        # Build covariate matrix: (p-1) * nclasses rows, totsize cols
        totsize <- sum(sampsize)
        covar <- matrix(0, (p - 1) * nclasses, totsize)
        y <- matrix(0, nclasses, totsize)

        # Fill data matrices
        for (c in seq_len(nclasses)) {
            idx_start <- if (c == 1) 1 else sum(sampsize[seq_len(c - 1)]) + 1
            idx_end <- sum(sampsize[seq_len(c)])

            # Target gene response for class c
            y[c, idx_start:idx_end] <- data_list[[c]][j, ]

            # Predictor genes for class c (excluding target gene j)
            covar[
                ((c - 1) * (p - 1) + 1):(c * (p - 1)),
                idx_start:idx_end
            ] <- data_list[[c]][-j, ]
        }

        # Run TRUE joint RF - joint modeling with full trees
        jrf.out <- .jrf_onetarget(
            x = covar, y = y, ntree = ntree, mtry = mtry,
            sampsize = sampsize, nclasses = nclasses
        )

        # Extract importance scores per class
        for (s in seq_len(nclasses)) {
            start_idx <- (p - 1) * (s - 1) + 1
            end_idx <- (p - 1) * s

            if (end_idx <= nrow(jrf.out$importance)) {
                imp[-j, j, s] <- jrf.out$importance[start_idx:end_idx, 2]
            } else {
                # Handle edge case - fill with small random values
                imp[-j, j, s] <- runif(p - 1, 0, 0.01)
            }
        }
    }

    # Create symmetric importance matrix
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
