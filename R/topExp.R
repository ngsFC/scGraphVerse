topExp <- function(seurat_object, cell_type, mart, pathway_genes, n = 1000) {
  
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
  
  # Compute average expression
  avg_expression <- rowMeans(normalized_data[, cells_in_type, drop = FALSE])
  
  top_ensembl_ids <- names(sort(avg_expression, decreasing = TRUE))[1:n]
  gene_mapping <- tryCatch(
    {
      biomaRt::getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol"),
        filters = "ensembl_gene_id",
        values = top_ensembl_ids,
        mart = mart
      )
    },
    error = function(e) {
      stop("Error querying biomaRt: ", e$message)
    }
  )
  
  # Map Ensembl IDs to gene names
  gene_mapping <- gene_mapping[gene_mapping$ensembl_gene_id %in% top_ensembl_ids, ]
  top_gene_names <- gene_mapping$hgnc_symbol[gene_mapping$hgnc_symbol != ""]
  
  # Filter by pathway or functional genes
  if (!is.null(pathway_genes)) {
    pathway_filtered_genes <- intersect(top_gene_names, pathway_genes)
    message("Top genes filtered by pathway: ", length(pathway_filtered_genes))
    print(pathway_filtered_genes)
  } else {
    message("No pathway filtering applied. Returning top genes.")
    pathway_filtered_genes <- top_gene_names
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
        title = paste("Expression Distribution of Selected Genes in", cell_type),
        x = "Gene",
        y = "Expression Level"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(plot)
  } else {
    warning("No genes matched the pathway filter for plotting.")
  }
  
  return(pathway_filtered_genes)
}
