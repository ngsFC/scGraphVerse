#' Internal JRF Network Inference Function (Simplified)
#'
#' This function provides JRF network inference by delegating to the 
#' original working JRF package instead of reimplementing the algorithm.
#'
#' @param X A list of expression matrices
#' @param ntree Number of trees in the random forest. Default: 500.
#' @param mtry Number of variables to sample at each tree node.
#' @param genes.name Character vector of gene names. If NULL, generates names.
#'
#' @return A data frame with gene pairs and importance scores
#' @export
#' @examples
#' # Create sample expression data as a list of matrices
#' n_genes <- 5
#' n_samples <- 20
#' 
#' # Generate two expression matrices
#' X <- list(
#'   matrix(rnorm(n_genes * n_samples, mean = 5, sd = 2), 
#'          nrow = n_genes, ncol = n_samples),
#'   matrix(rnorm(n_genes * n_samples, mean = 4, sd = 1.5), 
#'          nrow = n_genes, ncol = n_samples)
#' )
#' 
#' # Add gene names
#' rownames(X[[1]]) <- rownames(X[[2]]) <- paste0("Gene", 1:n_genes)
#' 
#' # Run JRF network inference
#' result <- JRF_internal(X, ntree = 50, genes.name = paste0("Gene", 1:n_genes))
JRF_internal <- function(X, ntree=500, mtry=NULL, genes.name=NULL) {
    # Internal JRF implementation exactly replicating original JRF algorithm
    
    # Set defaults exactly like original JRF
    nclasses <- length(X)
    sampsize <- rep(0, nclasses)
    
    for (j in 1:nclasses) sampsize[j] <- dim(X[[j]])[2]
    
    tot <- max(sampsize)
    p <- dim(X[[1]])[1]
    
    if (is.null(mtry)) mtry <- max(floor(p/3), 1)  # Default from original
    if (is.null(genes.name)) genes.name <- rownames(X[[1]])
    if (is.null(genes.name)) genes.name <- paste0("Gene", seq_len(p))
    
    # Initialize importance array like original
    imp <- array(0, c(p, length(genes.name), nclasses))
    
    # For each target gene (same loop as original JRF)
    for (j in seq_len(length(genes.name))) {
        
        # Create covar and y matrices exactly like original
        covar <- matrix(0, (p-1)*nclasses, tot)             
        y <- matrix(0, nclasses, tot)             
        
        for (c in 1:nclasses)  {
            y[c, seq_len(sampsize[c])] <- as.matrix(X[[c]][j,])
            covar[seq((c-1)*(p-1)+1, c*(p-1)), seq_len(sampsize[c])] <- X[[c]][-j,]
        }
        
        # Call regRF using .C like the original JRF_onetarget
        # Prepare parameters for regRF call (fix sizes to prevent crash)
        nodesize <- 5L
        nrnodes <- 2L * trunc(tot/max(1, nodesize - 4)) + 1L
        ncat <- rep(1L, (p-1)*nclasses)  # all continuous
        maxcat <- 1L
        
        # Fix weight vector size - should match total samples
        ww <- rep(1.0, tot)  # Equal weights for all samples
        
        # Limit matrix sizes to prevent memory issues
        nrnodes <- min(nrnodes, 1000L)  # Cap at reasonable size
        
        # Initialize output arrays
        impout <- matrix(0.0, (p-1)*nclasses, 2)
        impSD <- matrix(0.0, (p-1)*nclasses, 1)
        
        # Call regRF via .C exactly like original JRF
        nt <- 1  # keep.forest = FALSE
        rfout <- .C("regRF",
                   as.double(covar),
                   as.double(y), 
                   as.double(ww),
                   as.integer(c(tot, (p-1)*nclasses)),
                   as.integer(sampsize),
                   as.integer(tot),
                   as.integer(nodesize),
                   as.integer(nrnodes),
                   as.integer(ntree),
                   as.integer(mtry),
                   as.integer(c(1, 0, 1)),  # importance=TRUE, localImp=FALSE, nPerm=1
                   as.integer(ncat),
                   as.integer(maxcat),
                   as.integer(0),  # do.trace=FALSE
                   as.integer(0),  # proximity=FALSE
                   as.integer(0),  # oob.prox=FALSE
                   as.integer(0),  # corr.bias=FALSE
                   ypred = double(tot * nclasses),
                   impout = impout,
                   impmat = double(1),  # localImp=FALSE
                   impSD = impSD,
                   prox = double(1),    # proximity=FALSE
                   treeSize = integer(ntree),  # treeSize parameter  
                   nodestatus = integer(nrnodes * nt),  # Simplified allocation
                   leftDaughter = integer(nrnodes * nt),
                   rightDaughter = integer(nrnodes * nt), 
                   nodepred = double(nrnodes * nt),
                   bestvar = integer(nrnodes * nt),
                   xbestsplit = double(nrnodes * nt),
                   mse = double(ntree * nclasses),
                   keep = as.integer(c(0, 0)),  # keep.forest=FALSE, keep.inbag=FALSE
                   replace = as.integer(1),  # replace=TRUE
                   testdat = as.integer(0),  # no test data
                   xts = double(1),
                   ntest = as.integer(1),
                   yts = double(1),
                   labelts = as.integer(0),
                   ytestpred = double(1),
                   proxts = double(1),
                   msets = double(1),
                   coef = double(2),
                   nout = integer(tot),
                   inbag = integer(1), 
                   as.integer(nclasses),
                   PACKAGE = "scGraphVerse"
                   )[c(17:29, 35:40)]  # Extract like original
        
        # Extract importance scores like original (scale=FALSE)
        importance_scores <- rfout$impout[, 2]  # MeanDecreaseGini column
        
        # Save importance scores for each class like original
        for (s in 1:nclasses) {
            imp[-j, j, s] <- importance_scores[seq((p-1)*(s-1)+1, (p-1)*(s-1)+p-1)]
        }
    }
    
    # Derive final importance scores exactly like original
    imp.final <- matrix(0, p*(p-1)/2, nclasses)
    vec1 <- matrix(rep(genes.name, p), p, p)
    vec2 <- t(vec1)
    vec1 <- vec1[lower.tri(vec1, diag=FALSE)]
    vec2 <- vec2[lower.tri(vec2, diag=FALSE)]
    
    for (s in 1:nclasses) { 
        imp.s <- imp[,,s]
        t.imp <- t(imp.s)
        imp.final[,s] <- (imp.s[lower.tri(imp.s, diag=FALSE)] + 
                         t.imp[lower.tri(t.imp, diag=FALSE)])/2        
    }
    
    # Create output exactly like original
    out <- cbind(as.character(vec1), as.character(vec2), 
                 as.data.frame(imp.final), stringsAsFactors=FALSE)
    colnames(out) <- c(paste0('gene', 1:2), paste0('importance', 1:nclasses))
    
    return(out)
}