#' Extract Top Highly Expressed Genes for a Specific Cell Type in a Seurat Object
#'
#' This function identifies and extracts the top `n` most highly expressed genes 
#' for a given cell type from the RNA assay in a Seurat object. It computes the 
#' average expression of each gene across all cells belonging to the specified 
#' cell type and selects the highest-expressed ones.
#'
#' @param seurat_object A Seurat object containing RNA expression data.
#' @param cell_type A character string specifying the cell type of interest.
#' @param top_n An integer indicating the number of top expressed genes to select.
#' @return A character vector of the selected top `n` highly expressed genes.
#' @details 
#' The function first validates the presence of the specified `cell_type` in the 
#' metadata and extracts the corresponding cells. It then retrieves the normalized 
#' gene expression matrix from the `"RNA"` assay and calculates the mean expression 
#' of each gene across the selected cells. The genes are sorted in descending order 
#' based on their average expression, and the top `n` genes are returned.
#'
#' @examples
#' \dontrun{
#' selected_genes <- pathg(my_seurat_object, cell_type = "T_cells", top_n = 20)
#' print(selected_genes)
#' }
#' 
#' @export
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

