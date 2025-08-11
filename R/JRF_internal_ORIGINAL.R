# =========================
# TRUE JOINT JRF Implementation (Mathematically Equivalent to Original)
# =========================

.jrf_onetarget_joint <- function(x, y, ntree, mtry, sampsize, nclasses,
                                importance = TRUE) {
    # TRUE joint modeling - single C call with full data
    # This matches the original JRF mathematical formulation exactly
    
    totsize <- sum(sampsize)
    p_per_class <- nrow(x) / nclasses  # p-1 predictors per class
    
    # Prepare data structures exactly as in original JRF
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    
    # Original JRF parameter setup
    n <- totsize  # total sample size across all classes
    p <- p_per_class  # predictors per class (p-1 in the context)
    
    # Key parameters matching original JRF_onetarget
    ww <- rep(1.0, nclasses)  # weights per class
    storage.mode(ww) <- "double"
    
    purity <- 0.0  # default purity parameter
    nodesize <- 5
    nrnodes <- 2 * floor(n / max(1, nodesize - 4)) + 1
    
    # Importance matrices - dimensions match original exactly
    impout <- matrix(0.0, p * nclasses, 2)
    impSD <- matrix(0.0, p * nclasses, 1)
    
    # Other parameters
    ncat <- rep(1, p * nclasses)  # all continuous variables
    maxcat <- 1
    
    # Trees and forest storage (keep.forest = TRUE equivalent)
    nt <- ntree  # number of trees to keep
    
    # Call regRF with JOINT data structure - single call for all classes
    rfout <- .C("regRF",
                x,                                            # 1. double *x: joint covariate matrix
                y,                                            # 2. double *y: joint response matrix  
                ww,                                           # 3. double *weight: class weights
                as.double(purity),                            # 4. double *purity
                as.integer(c(totsize, p)),                    # 5. int *xdim: (total_samples, predictors_per_class)
                as.integer(totsize),                          # 6. int *sampsize: total sample size
                as.integer(nodesize),                         # 7. int *nthsize
                as.integer(nrnodes),                          # 8. int *nrnodes
                as.integer(ntree),                            # 9. int *nTree
                as.integer(mtry),                             # 10. int *mtry
                as.integer(c(as.integer(importance), 0, 1)),  # 11. int *imp: (importance, localImp, nPerm)
                as.integer(ncat),                             # 12. int *cat: categorical indicators
                as.integer(maxcat),                           # 13. int *maxcat
                as.integer(0),                                # 14. int *jprint
                as.integer(0),                                # 15. int *doProx
                as.integer(0),                                # 16. int *oobprox
                as.integer(0),                                # 17. int *biasCorr
                ypred = double(n * nclasses),                 # 18. double *yptr
                impout = impout,                              # 19. double *errimp: importance output
                impmat = double(1),                           # 20. double *impmat: importance matrix  
                impSD = impSD,                                # 21. double *impSD: importance SD
                prox = double(1),                             # 22. double *prox: proximity
                treeSize = integer(ntree),                    # 23. int *treeSize
                nodestatus = matrix(integer(nrnodes * nt * nclasses), ncol = nt),  # 24. int *nodestatus
                lDaughter = matrix(integer(nrnodes * nt * nclasses), ncol = nt),   # 25. int *lDaughter  
                rDaughter = matrix(integer(nrnodes * nt * nclasses), ncol = nt),   # 26. int *rDaughter
                avnode = matrix(double(nrnodes * nt * nclasses), ncol = nt),       # 27. double *avnode
                mbest = matrix(integer(nrnodes * nt * nclasses), ncol = nt),       # 28. int *mbest
                upper = matrix(double(nrnodes * nt * nclasses), ncol = nt),        # 29. double *upper
                mse = double(ntree * nclasses),               # 30. double *mse
                keepf = as.integer(c(1, 0)),                  # 31. int *keepf: (keep.forest, keep.inbag)
                replace = as.integer(1),                      # 32. int *replace
                testdat = as.integer(0),                      # 33. int *testdat
                xts = double(1),                              # 34. double *xts
                nts = as.integer(1),                          # 35. int *nts
                yts = double(1),                              # 36. double *yts
                labelts = as.integer(0),                      # 37. int *labelts
                yTestPred = double(1),                        # 38. double *yTestPred
                proxts = double(1),                           # 39. double *proxts
                msets = double(1),                            # 40. double *msets
                coef = double(2),                             # 41. double *coef
                nout = integer(n),                            # 42. int *nout (was oob.times)
                inbag = integer(1),                           # 43. int *inbag
                nclasses = as.integer(nclasses)               # 44. int *nclasses: ACTUAL number of classes
    )[c("impout")]
    
    # Return in randomForest format
    out <- list(importance = rfout$impout)
    class(out) <- "randomForest"
    return(out)
}


# Main JRF network inference with TRUE joint modeling
.jrf_network_joint <- function(data_list, ntree = 1000, mtry = NULL) {
    nclasses <- length(data_list)
    sampsize <- vapply(data_list, ncol, FUN.VALUE = integer(1))
    p <- nrow(data_list[[1]])
    if (is.null(mtry)) mtry <- max(floor(sqrt(p - 1)), 1)

    genes.name <- rownames(data_list[[1]])
    if (is.null(genes.name)) genes.name <- paste0("G", seq_len(p))

    imp <- array(0, dim = c(p, length(genes.name), nclasses))

    for (j in seq_len(p)) {
        # Build covariate matrix: (p-1) * nclasses rows, totsize cols
        # This is IDENTICAL to original JRF data structure
        totsize <- sum(sampsize)
        covar <- matrix(0, (p - 1) * nclasses, totsize)
        y <- matrix(0, nclasses, totsize)

        # Fill data matrices - IDENTICAL to original JRF
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

        # Run TRUE joint RF - SINGLE C call with cross-class modeling
        jrf.out <- .jrf_onetarget_joint(
            x = covar, y = y, ntree = ntree, mtry = mtry,
            sampsize = sampsize, nclasses = nclasses
        )

        # Extract importance scores per class - IDENTICAL to original
        for (s in seq_len(nclasses)) {
            start_idx <- (p - 1) * (s - 1) + 1
            end_idx <- (p - 1) * s

            if (end_idx <= nrow(jrf.out$importance)) {
                # Use column 2 for importance (matches original importance() function)
                imp[-j, j, s] <- jrf.out$importance[start_idx:end_idx, 2]
            } else {
                warning("Importance matrix dimensions mismatch for gene ", j, ", class ", s)
                imp[-j, j, s] <- rep(0, p - 1)
            }
        }
    }

    # Create symmetric importance matrix - IDENTICAL to original
    imp.final <- matrix(0, nrow = p * (p - 1) / 2, ncol = nclasses)

    # Generate gene pair names - IDENTICAL to original
    vec1 <- matrix(rep(genes.name, p), p, p)
    vec2 <- t(vec1)
    vec1 <- vec1[lower.tri(vec1, diag = FALSE)]
    vec2 <- vec2[lower.tri(vec2, diag = FALSE)]

    # Aggregate importance scores (symmetric) - IDENTICAL to original
    for (s in seq_len(nclasses)) {
        imp.s <- imp[, , s]
        t.imp <- t(imp.s)
        imp.final[, s] <- (imp.s[lower.tri(imp.s, diag = FALSE)] +
                          t.imp[lower.tri(t.imp, diag = FALSE)]) / 2
    }

    # Create result list - one network per class
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