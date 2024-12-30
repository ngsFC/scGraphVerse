pathg <- function(seurat_object, cell_type, mart, pathway_genes) {
  
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
  
  # Filter by pathway or functional genes
  if (!is.null(pathway_genes)) {
    pathway_filtered_genes <- intersect(rownames(normalized_data), pathway_genes)
    message("Genes found in the pathway: ", length(pathway_filtered_genes))
    print(pathway_filtered_genes)
  } else {
    stop("No pathway genes provided. Please supply a list of pathway genes.")
  }
  
  if (length(pathway_filtered_genes) > 0) {
    selected_expression <- normalized_data[pathway_filtered_genes, cells_in_type, drop = FALSE]
    melted_data <- melt(selected_expression)
    colnames(melted_data) <- c("Gene", "Cell", "Expression")
    
    plot <- ggplot(melted_data, aes(x = Gene, y = Expression)) +
      geom_violin(trim = TRUE, fill = "lightblue", color = "darkblue") +
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
      theme_minimal() +
      labs(
        title = paste("Expression Distribution of Pathway Genes in", cell_type),
        x = "Gene",
        y = "Expression Level"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(plot)
  } else {
    warning("No genes matched the pathway for plotting.")
  }
  
  return(pathway_filtered_genes)
}
