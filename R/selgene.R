#' Extract Top Highly Expressed Genes
#'
#' Identifies and extracts the top \code{n} most highly expressed genes across all cells
#' from a \linkS4class{Seurat} object, a \linkS4class{SingleCellExperiment} object,
#' or a numeric expression matrix.
#'
#' @param object   A Seurat object, SingleCellExperiment object, or numeric matrix (genes × cells).
#' @param top_n    Integer. Number of top expressed genes to return.
#' @param assay    (SCE only) Name of the assay to use; if NULL, tries "logcounts" then "counts".
#'
#' @return
#' A character vector of the top \code{n} most highly expressed genes.
#'
#' @export
selgene <- function(object, top_n, assay = NULL) {
  # validate top_n
  if (missing(top_n) || length(top_n) != 1 || !is.numeric(top_n) || top_n <= 0) {
    stop("Please provide a valid 'top_n' (a single positive integer).")
  }
  
  ### Seurat path: pick 'counts' slot first
  if (inherits(object, "Seurat")) {
    # determine default assay name
    assay_name <- Seurat::DefaultAssay(object)
    seurat_assay <- object[[assay_name]]
    # list available slots
    slots_avail <- methods::slotNames(seurat_assay)
    if ("counts" %in% slots_avail) {
      expr <- seurat_assay@counts
      message("Using Seurat assay '", assay_name, "' slot 'counts'.")
    } else if ("data" %in% slots_avail) {
      expr <- seurat_assay@data
      message("No 'counts' slot in assay '", assay_name,
              "'; using 'data' (normalized) instead.")
    } else {
      stop("Assay '", assay_name, "' has neither 'counts' nor 'data' slots. ",
           "Available slots: ", paste(slots_avail, collapse = ", "))
    }
    
    ### SingleCellExperiment path
  } else if (inherits(object, "SingleCellExperiment")) {
    available_assays <- SummarizedExperiment::assayNames(object)
    # choose assay name
    if (is.null(assay)) {
      if ("logcounts" %in% available_assays) {
        assay_to_use <- "logcounts"
      } else if ("counts" %in% available_assays) {
        assay_to_use <- "counts"
      } else {
        stop("No 'logcounts' or 'counts' assays found. Available assays: ",
             paste(available_assays, collapse = ", "))
      }
    } else {
      if (! assay %in% available_assays) {
        stop("Requested assay '", assay, "' not found. Available assays: ",
             paste(available_assays, collapse = ", "))
      }
      assay_to_use <- assay
    }
    expr <- SummarizedExperiment::assay(object, assay_to_use)
    message("Using SCE assay '", assay_to_use, "'.")
    
    ### Matrix path
  } else if (is.matrix(object)) {
    expr <- object
    
  } else {
    stop("Input must be a Seurat object, SingleCellExperiment object, or numeric matrix.")
  }
  
  # clean up duplicates
  expr <- as.matrix(expr)
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]
  
  # compute mean expression per gene
  avg_expression <- rowMeans(expr, na.rm = TRUE)
  
  # select top genes
  sorted_genes   <- names(sort(avg_expression, decreasing = TRUE))
  selected_genes <- head(sorted_genes, top_n)
  message("Top ", top_n, " genes selected based on mean expression.")
  
  return(selected_genes)
}
