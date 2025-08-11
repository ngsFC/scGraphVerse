# =========================
# HYBRID JRF Implementation - Mathematically Closer to Original 
# =========================
# This implementation addresses the C function limitations while maximizing
# mathematical similarity to the original JRF

.jrf_onetarget_hybrid <- function(x, y, ntree, mtry, sampsize, nclasses,
                                 importance = TRUE) {
    # Hybrid approach: Simulate joint modeling effects while using stable C calls
    # Key insight: Run each class but with AUGMENTED information from other classes
    
    totsize <- sum(sampsize)
    p_per_class <- nrow(x) / nclasses
    
    # Initialize combined importance matrix
    combined_importance <- matrix(0.0, p_per_class * nclasses, 2)
    
    # Key improvement: Use cross-class information for better approximation
    for (class_idx in seq_len(nclasses)) {
        # Extract primary class data
        class_start_col <- if (class_idx == 1) 1 else sum(sampsize[seq_len(class_idx - 1)]) + 1
        class_end_col <- sum(sampsize[seq_len(class_idx)])
        
        class_start_row <- (class_idx - 1) * p_per_class + 1
        class_end_row <- class_idx * p_per_class
        
        # Primary class data
        x_primary <- x[class_start_row:class_end_row, class_start_col:class_end_col]
        y_primary <- matrix(y[class_idx, class_start_col:class_end_col], nrow = 1)
        
        # HYBRID ENHANCEMENT: Augment with representative samples from other classes
        # This simulates the cross-class information sharing of true joint modeling
        if (nclasses > 1) {
            other_classes <- setdiff(1:nclasses, class_idx)
            
            # Sample representative data from other classes (10% of samples)
            aug_samples <- list()
            aug_responses <- list()
            aug_count <- 0
            
            for (other_class in other_classes) {
                other_start_col <- if (other_class == 1) 1 else sum(sampsize[seq_len(other_class - 1)]) + 1
                other_end_col <- sum(sampsize[seq_len(other_class)])
                other_start_row <- (other_class - 1) * p_per_class + 1
                other_end_row <- other_class * p_per_class
                
                # Sample a subset of other class data
                n_aug <- max(1, floor(0.1 * sampsize[other_class]))
                if (sampsize[other_class] > 1) {
                    sample_idx <- sample(other_start_col:other_end_col, n_aug)
                } else {
                    sample_idx <- other_start_col:other_end_col
                }
                
                aug_samples[[length(aug_samples) + 1]] <- x[other_start_row:other_end_row, sample_idx, drop = FALSE]
                aug_responses[[length(aug_responses) + 1]] <- matrix(y[other_class, sample_idx], nrow = 1)
                aug_count <- aug_count + length(sample_idx)
            }
            
            # Combine primary and augmented data
            if (aug_count > 0) {
                x_combined <- cbind(x_primary, do.call(cbind, aug_samples))
                y_combined <- cbind(y_primary, do.call(cbind, aug_responses))
                combined_sampsize <- ncol(x_combined)
            } else {
                x_combined <- x_primary
                y_combined <- y_primary
                combined_sampsize <- sampsize[class_idx]
            }
        } else {
            x_combined <- x_primary
            y_combined <- y_primary
            combined_sampsize <- sampsize[class_idx]
        }
        
        # Ensure proper storage mode
        storage.mode(x_combined) <- "double"
        storage.mode(y_combined) <- "double"
        
        # Adjusted parameters for combined data
        weight_class <- 1.0
        storage.mode(weight_class) <- "double"
        
        nodesize <- 5
        nrnodes <- 2 * floor(combined_sampsize / max(1, nodesize - 4)) + 1
        
        # Initialize output matrices for this class
        impout_class <- matrix(0.0, p_per_class, 2)
        impSD_class <- matrix(0.0, p_per_class, 1)
        
        # ENHANCED C call with augmented data
        rfout <- .C("regRF",
                    x_combined,                                           # 1. double *x
                    y_combined,                                           # 2. double *y  
                    weight_class,                                         # 3. double *weight
                    as.integer(c(combined_sampsize, p_per_class)),        # 4. int *xdim 
                    as.integer(combined_sampsize),                        # 5. int *sampsize
                    as.integer(combined_sampsize),                        # 6. int *totsize
                    as.integer(nodesize),                                 # 7. int *nthsize
                    as.integer(nrnodes),                                  # 8. int *nrnodes
                    as.integer(ntree),                                    # 9. int *nTree
                    as.integer(mtry),                                     # 10. int *mtry
                    as.integer(c(as.integer(importance), 0, 1)),          # 11. int *imp
                    as.integer(rep(1, p_per_class)),                      # 12. int *cat
                    as.integer(1),                                        # 13. int *maxcat
                    as.integer(0),                                        # 14. int *jprint
                    as.integer(0),                                        # 15. int *doProx
                    as.integer(0),                                        # 16. int *oobprox
                    as.integer(0),                                        # 17. int *biasCorr
                    ypred = double(combined_sampsize),                    # 18. double *yptr
                    impout = impout_class,                                # 19. double *errimp
                    impmat = double(1),                                   # 20. double *impmat
                    impSD = impSD_class,                                  # 21. double *impSD
                    prox = double(1),                                     # 22. double *prox
                    treeSize = integer(ntree),                            # 23. int *treeSize
                    nodestatus = integer(nrnodes * ntree),                # 24. int *nodestatus
                    lDaughter = integer(nrnodes * ntree),                 # 25. int *lDaughter
                    rDaughter = integer(nrnodes * ntree),                 # 26. int *rDaughter
                    avnode = double(nrnodes * ntree),                     # 27. double *avnode
                    mbest = integer(nrnodes * ntree),                     # 28. int *mbest
                    upper = double(nrnodes * ntree),                      # 29. double *upper
                    mse = double(ntree),                                  # 30. double *mse
                    keepf = as.integer(c(0, 0)),                          # 31. int *keepf
                    replace = as.integer(1),                              # 32. int *replace
                    testdat = as.integer(0),                              # 33. int *testdat
                    xts = double(1),                                      # 34. double *xts
                    nts = as.integer(1),                                  # 35. int *nts
                    yts = double(1),                                      # 36. double *yts
                    labelts = as.integer(0),                              # 37. int *labelts
                    yTestPred = double(1),                                # 38. double *yTestPred
                    proxts = double(1),                                   # 39. double *proxts
                    msets = double(1),                                    # 40. double *msets
                    coef = double(2),                                     # 41. double *coef
                    nout = integer(combined_sampsize),                    # 42. int *nout
                    inbag = integer(1),                                   # 43. int *inbag
                    nclasses = as.integer(1)                              # 44. int *nclasses (1 per call)
        )[c("impout")]
        
        # Store importance for this class
        combined_importance[class_start_row:class_end_row, ] <- rfout$impout
    }
    
    out <- list(importance = combined_importance)
    class(out) <- "randomForest"
    return(out)
}


.jrf_network_hybrid <- function(data_list, ntree = 1000, mtry = NULL) {
    nclasses <- length(data_list)
    sampsize <- vapply(data_list, ncol, FUN.VALUE = integer(1))
    p <- nrow(data_list[[1]])
    if (is.null(mtry)) mtry <- max(floor(sqrt(p - 1)), 1)

    genes.name <- rownames(data_list[[1]])
    if (is.null(genes.name)) genes.name <- paste0("G", seq_len(p))

    imp <- array(0, dim = c(p, length(genes.name), nclasses))

    for (j in seq_len(p)) {
        # Build covariate matrix: (p-1) * nclasses rows, totsize cols
        # IDENTICAL to original JRF data structure
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

        # Run HYBRID RF - Enhanced cross-class modeling
        jrf.out <- .jrf_onetarget_hybrid(
            x = covar, y = y, ntree = ntree, mtry = mtry,
            sampsize = sampsize, nclasses = nclasses
        )

        # Extract importance scores per class - IDENTICAL to original
        for (s in seq_len(nclasses)) {
            start_idx <- (p - 1) * (s - 1) + 1
            end_idx <- (p - 1) * s

            if (end_idx <= nrow(jrf.out$importance)) {
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