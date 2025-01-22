pathg <- function(seurat_object, cell_type, pathway_genes = NULL, top_n = NULL) {
  
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
  
  # Extract normalized RNA expression data
  normalized_data <- as.matrix(seurat_object[["RNA"]]@data)
  
  # Remove duplicate row names
  normalized_data <- normalized_data[!duplicated(rownames(normalized_data)), , drop = FALSE]
  # Remove duplicate column names
  normalized_data <- normalized_data[, !duplicated(colnames(normalized_data)), drop = FALSE]
  
  # Determine the genes to use
  if (!is.null(pathway_genes)) {
    # Use pathway genes if provided
    pathway_filtered_genes <- intersect(rownames(normalized_data), pathway_genes)
    message("Genes found in the pathway: ", length(pathway_filtered_genes))
    print(pathway_filtered_genes)
    
    if (length(pathway_filtered_genes) == 0) {
      warning("No genes matched the pathway. Returning an empty list.")
      return(NULL)
    }
    selected_genes <- pathway_filtered_genes
    
  } else if (!is.null(top_n)) {
    # Select top n expressed genes if specified
    avg_expression <- rowMeans(normalized_data[, cells_in_type, drop = FALSE])
    sorted_genes <- names(sort(avg_expression, decreasing = TRUE))
    selected_genes <- head(sorted_genes, top_n)
    message("Top ", top_n, " genes selected based on mean expression.")
    
  } else {
    stop("Neither pathway genes nor top_n parameter provided. Please specify one.")
  }
  
  # Extract expression data for the selected genes
  selected_expression <- normalized_data[selected_genes, cells_in_type, drop = FALSE]
  
  # Remove duplicated genes (rows) based on maximum expression value
  selected_expression <- selected_expression[!duplicated(rownames(selected_expression)), , drop = FALSE]
  
  # Melt the data for ggplot visualization
  melted_data <- reshape2::melt(selected_expression)
  colnames(melted_data) <- c("Gene", "Cell", "Expression")
  
  # Plot the average expression per gene
  avg_expression <- rowMeans(selected_expression)
  sorted_genes <- names(sort(avg_expression, decreasing = TRUE))
  
  return(selected_genes)
}
