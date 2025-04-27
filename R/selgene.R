#' Extract Top Highly Expressed Genes for a Specific Cell Type
#'
#' Identifies and extracts the top \code{n} most highly expressed genes for a specified
#' cell type from the RNA assay of a \linkS4class{Seurat} object.
#'
#' @param seurat_object A \linkS4class{Seurat} object containing RNA expression data.
#' @param cell_type Character string specifying the target cell type of interest.
#' @param top_n Integer. Number of top expressed genes to return.
#'
#' @return
#' A character vector containing the top \code{n} most highly expressed genes for the specified cell type.
#'
#' @details
#' The function verifies that the specified \code{cell_type} exists in the metadata of the \code{Seurat} object.
#' It retrieves the normalized expression matrix from the \code{"RNA"} assay, computes mean expression per gene
#' across all cells of the specified type, and returns the genes with the highest average expression.
#'
#' @note
#' Requires the \pkg{Seurat} package.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Select top 20 genes for the "T_cells" cluster
#' selected_genes <- selgene(my_seurat_object, cell_type = "T_cells", top_n = 20)
#' print(selected_genes)
#' }
selgene <- function(seurat_object, cell_type, top_n = NULL) {
  # Validate cell type
  if (!"cell_type" %in% colnames(seurat_object@meta.data)) {
    stop("'cell_type' column not found in the Seurat metadata.")
  }
  if (!cell_type %in% seurat_object@meta.data$cell_type) {
    stop("The selected cell type is not present in the metadata.")
  }

  cells_in_type <- which(seurat_object@meta.data$cell_type == cell_type)
  if (length(cells_in_type) == 0) {
    stop("No cells found for the selected cell type.")
  }

  normalized_data <- as.matrix(seurat_object[["RNA"]]@data)
  normalized_data <- normalized_data[!duplicated(rownames(normalized_data)), , drop = FALSE]
  normalized_data <- normalized_data[, !duplicated(colnames(normalized_data)), drop = FALSE]

  if (is.null(top_n)) {
    stop("The 'top_n' parameter must be provided to select highly expressed genes.")
  }

  avg_expression <- rowMeans(normalized_data[, cells_in_type, drop = FALSE])
  sorted_genes <- names(sort(avg_expression, decreasing = TRUE))
  selected_genes <- head(sorted_genes, top_n)
  message("Top ", top_n, " genes selected based on mean expression.")

  return(selected_genes)
}
