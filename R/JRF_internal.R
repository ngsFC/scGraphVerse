#' Internal JRF Network Inference Function (Simplified)
#'
#' This function provides JRF network inference by delegating to the 
#' original working JRF package instead of reimplementing the algorithm.
#'
#' @param X A list of expression matrices
#' @param ntree Number of trees in the random forest. Default: 500.
#' @param mtry Number of variables to sample at each tree node. If NULL, defaults to sqrt(p).
#' @param genes.name Character vector of gene names. If NULL, generates names.
#'
#' @return A data frame with gene pairs and importance scores
#' @export
JRF_internal <- function(X, ntree=500, mtry=NULL, genes.name=NULL) {
    # SOLUTION: Use the working original JRF package directly
    # scGraphVerse's JRF_onetarget implementation has bugs, so delegate to original
    
    # Set defaults exactly like original
    p <- dim(X[[1]])[1]
    if (is.null(mtry)) mtry <- floor(sqrt(p))
    if (is.null(genes.name)) genes.name <- rownames(X[[1]])
    if (is.null(genes.name)) genes.name <- paste0("Gene", 1:p)
    
    # Check if JRF package is available
    if (!requireNamespace("JRF", quietly = TRUE)) {
        stop("JRF package is required but not available. Please install it.")
    }
    
    # Call the original working JRF function
    result <- JRF::JRF(X, ntree=ntree, mtry=mtry, genes.name=genes.name)
    
    return(result)
}