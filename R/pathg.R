#' Select and Extract Highly Expressed Genes in a Seurat Object
#'
#' This function identifies and extracts the top `n` most highly expressed genes 
#' for a given cell type in a Seurat object.
#'
#' @param seurat_object A Seurat object containing RNA expression data.
#' @param cell_type A character string specifying the cell type of interest.
#' @param top_n An integer indicating the number of top expressed genes to select.
#' @return A character vector of the selected genes.
#' @examples
#' \dontrun{
#' selected_genes <- pathg(my_seurat_object, cell_type = "T_cells", top_n = 20)
#' }
#' @export
pathg <- function(seurat_object, cell_type, top_n = NULL) {
  
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

