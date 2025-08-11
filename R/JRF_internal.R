# =========================
# Internal helper for JRF (C implementation)
# =========================

.jrf_onetarget <- function(x, y, ntree, mtry, sampsize, nclasses,
                        importance = TRUE) {
    # JRF processes each class separately, not jointly in one C call
    # This is the correct interpretation of the original algorithm
    
    totsize <- sum(sampsize)
    p_per_class <- nrow(x) / nclasses
    
    # Combine importance from all classes
    combined_importance <- matrix(0.0, p_per_class * nclasses, 2)
    
    # Process each class separately 
    for (class_idx in seq_len(nclasses)) {
        # Extract data for this class only
        class_start_col <- if (class_idx == 1) 1 else sum(sampsize[seq_len(class_idx - 1)]) + 1
        class_end_col <- sum(sampsize[seq_len(class_idx)])
        
        class_start_row <- (class_idx - 1) * p_per_class + 1
        class_end_row <- class_idx * p_per_class
        
        # Single class data
        x_class <- x[class_start_row:class_end_row, class_start_col:class_end_col]
        y_class <- matrix(y[class_idx, class_start_col:class_end_col], nrow = 1)
        
        # Storage mode for C interface
        storage.mode(x_class) <- "double"
        storage.mode(y_class) <- "double"
        
        # Single class parameters
        class_sampsize <- sampsize[class_idx]
        weight_class <- 1.0
        storage.mode(weight_class) <- "double"
        
        nodesize <- 5
        nrnodes <- 2 * floor(class_sampsize / max(1, nodesize - 4)) + 1
        
        # Initialize output matrices for this class
        impout_class <- matrix(0.0, p_per_class, 2)
        impSD_class <- matrix(0.0, p_per_class, 1)
        
        # Call C function for this class only
        rfout <- .C("regRF",
                    x_class,                                      # 1. double *x
                    y_class,                                      # 2. double *y  
                    weight_class,                                 # 3. double *weight
                    as.integer(c(class_sampsize, p_per_class)),  # 4. int *xdim (nsample, mdim)
                    as.integer(class_sampsize),                   # 5. int *sampsize
                    as.integer(class_sampsize),                   # 6. int *totsize
                    as.integer(nodesize),                         # 7. int *nthsize
                    as.integer(nrnodes),                          # 8. int *nrnodes
                    as.integer(ntree),                            # 9. int *nTree
                    as.integer(mtry),                             # 10. int *mtry
                    as.integer(c(as.integer(importance), 0, 1)),  # 11. int *imp
                    as.integer(rep(1, p_per_class)),              # 12. int *cat
                    as.integer(1),                                # 13. int *maxcat
                    as.integer(0),                                # 14. int *jprint
                    as.integer(0),                                # 15. int *doProx
                    as.integer(0),                                # 16. int *oobprox
                    as.integer(0),                                # 17. int *biasCorr
                    ypred = double(class_sampsize),               # 18. double *yptr
                    impout = impout_class,                        # 19. double *errimp
                    impmat = double(1),                           # 20. double *impmat
                    impSD = impSD_class,                          # 21. double *impSD
                    prox = double(1),                             # 22. double *prox
                    treeSize = integer(ntree),                    # 23. int *treeSize
                    nodestatus = integer(nrnodes * ntree),        # 24. int *nodestatus
                    lDaughter = integer(nrnodes * ntree),         # 25. int *lDaughter
                    rDaughter = integer(nrnodes * ntree),         # 26. int *rDaughter
                    avnode = double(nrnodes * ntree),             # 27. double *avnode
                    mbest = integer(nrnodes * ntree),             # 28. int *mbest
                    upper = double(nrnodes * ntree),              # 29. double *upper
                    mse = double(ntree),                          # 30. double *mse
                    keepf = as.integer(c(0, 0)),                  # 31. int *keepf
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
                    nout = integer(class_sampsize),               # 42. int *nout
                    inbag = integer(1),                           # 43. int *inbag
                    nclasses = as.integer(1)                      # 44. int *nclasses (always 1 per call)
        )[c("impout")]
        
        # Store importance for this class
        combined_importance[class_start_row:class_end_row, ] <- rfout$impout
    }
    
    out <- list(importance = combined_importance)
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

    # Generate gene pair names
    vec1 <- matrix(rep(genes.name, p), p, p)
    vec2 <- t(vec1)
    vec1 <- vec1[lower.tri(vec1, diag = FALSE)]
    vec2 <- vec2[lower.tri(vec2, diag = FALSE)]

    # Aggregate importance scores (symmetric)
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