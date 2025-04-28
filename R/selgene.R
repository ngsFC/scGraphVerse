#' Extract Top Highly Expressed Genes
#'
#' Identifies and extracts the top \code{n} most highly expressed genes across all cells
#' from a \linkS4class{Seurat} object, a \linkS4class{SingleCellExperiment} object, or a numeric matrix.
#'
#' @param object A \linkS4class{Seurat} object, a \linkS4class{SingleCellExperiment} object, or a numeric expression matrix (genes x cells).
#' @param top_n Integer. Number of top expressed genes to return.
#'
#' @return
#' A character vector containing the top \code{n} most highly expressed genes.
#'
#' @details
#' The function automatically detects the input type, extracts normalized expression data,
#' computes the mean expression per gene across all cells, and selects the top genes based on average expression.
#'
#' @note
#' Requires the \pkg{Seurat} and \pkg{SingleCellExperiment} packages if S4 objects are provided.
#'
#' @export
#'
#' @examples
#' counts <- matrix(rpois(100, lambda = 5), nrow = 10)
#' rownames(counts) <- paste0("Gene", 1:10)
#' colnames(counts) <- paste0("Cell", 1:10)
#' selected_genes <- selgene(counts, top_n = 5)

selgene <- function(object, top_n) {
  if (missing(top_n) || !is.numeric(top_n)) {
    stop("Please provide a valid 'top_n' value (positive integer).")
  }
  
  # Extract normalized data based on object class
  if (inherits(object, "Seurat")) {
    expr <- object[["RNA"]]@data
  } else if (inherits(object, "SingleCellExperiment")) {
    expr <- SummarizedExperiment::assay(object, "logcounts")
  } else if (is.matrix(object)) {
    expr <- object
  } else {
    stop("Input must be a Seurat object, SingleCellExperiment object, or an expression matrix.")
  }
  
  # Sanitize matrix
  expr <- as.matrix(expr)
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]
  
  # Compute mean expression per gene
  avg_expression <- rowMeans(expr, na.rm = TRUE)
  
  # Select top genes
  sorted_genes <- names(sort(avg_expression, decreasing = TRUE))
  selected_genes <- head(sorted_genes, top_n)
  message("Top ", top_n, " genes selected based on mean expression.")
  
  return(selected_genes)
}
