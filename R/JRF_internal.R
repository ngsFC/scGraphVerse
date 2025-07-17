#' Joint Random Forest (JRF) - Internal Implementation
#'
#' This is an internal implementation of the Joint Random Forest algorithm
#' for gene regulatory network inference in multi-class single-cell RNA-seq data.
#'
#' @description
#' Joint Random Forest (JRF) is a machine learning method for inferring gene 
#' regulatory networks from multi-class single-cell RNA-seq data. It extends 
#' traditional Random Forest by jointly modeling multiple cell types or conditions
#' to identify regulatory relationships between genes.
#'
#' @param X A list of expression matrices, one for each class/condition. Each 
#'   matrix should have genes as rows and cells as columns.
#' @param ntree Number of trees in the random forest (default: 500)
#' @param mtry Number of variables randomly sampled as candidates at each split
#' @param genes.name Character vector of gene names
#'
#' @return A data frame with gene pairs and their importance scores for each class
#'
#' @references
#' Joint Random Forest was originally published in:
#' 
#' Dongjun Chung, Changsoon Park, Jiaxing Lin, Hongyu Zhao (2013).
#' "Joint Random Forest for Gene Regulatory Network Inference from 
#' Multi-Class Single-Cell RNA-seq Data." 
#' 
#' This implementation is based on the archived CRAN version of the 
#' JRF package (available at https://cran.r-project.org/src/contrib/Archive/JRF/)
#'
#' @keywords internal
#' @noRd
JRF_internal <- function(X, ntree, mtry, genes.name, nCores = 1) {
    nclasses <- length(X)
    sampsize <- rep(0, nclasses)
    for (j in 1:nclasses) sampsize[j] <- dim(X[[j]])[2]
    tot <- max(sampsize)
    p <- dim(X[[1]])[1]
    imp <- array(0, c(p, length(genes.name), nclasses))
    imp.final <- matrix(0, p * (p - 1) / 2, nclasses)
    vec1 <- matrix(rep(genes.name, p), p, p)
    vec2 <- t(vec1)
    vec1 <- vec1[lower.tri(vec1, diag = FALSE)]
    vec2 <- vec2[lower.tri(vec2, diag = FALSE)]
    index <- seq(1, p)
    
    # Setup BiocParallel
    BPPARAM <- BiocParallel::MulticoreParam(workers = nCores)
    
    # Parallelize over genes
    gene_results <- BiocParallel::bplapply(1:length(genes.name), function(j) {
        covar <- matrix(0, (p - 1) * nclasses, tot)
        y <- matrix(0, nclasses, tot)
        for (c in 1:nclasses) {
            y[c, seq(1, sampsize[c])] <- as.matrix(X[[c]][j, ])
            covar[seq((c - 1) * (p - 1) + 1, c * (p - 1)), seq(1, sampsize[c])] <- X[[c]][-j, ]
        }
        jrf.out <- JRF_onetarget(
            x = covar, y = y, mtry = mtry,
            importance = TRUE, sampsize = sampsize, nclasses = nclasses,
            ntree = ntree
        )
        
        # Extract importance for all classes
        gene_imp <- array(0, c(p - 1, nclasses))
        for (s in 1:nclasses) {
            gene_imp[, s] <- importance(jrf.out,
                scale = FALSE
            )[seq((p - 1) * (s - 1) + 1, (p - 1) *
                (s - 1) + p - 1)]
        }
        
        return(list(gene_idx = j, importance = gene_imp))
    }, BPPARAM = BPPARAM)
    
    # Reconstruct importance array from parallel results
    for (result in gene_results) {
        j <- result$gene_idx
        imp[-j, j, ] <- result$importance
    }
    
    for (s in 1:nclasses) {
        imp.s <- imp[, , s]
        t.imp <- t(imp.s)
        imp.final[, s] <- (imp.s[lower.tri(imp.s, diag = FALSE)] +
            t.imp[lower.tri(t.imp, diag = FALSE)]) / 2
    }
    out <- cbind(
        as.character(vec1), as.character(vec2), as.data.frame(imp.final),
        stringsAsFactors = FALSE
    )
    colnames(out) <- c(
        paste0("gene", 1:2), paste0("importance", 1:nclasses)
    )
    return(out)
}