# Function to calculate top N genes expressed in a cell type, returning gene names
topExp <- function(seurat_object, cell_type, mart, n = 300) {
  # Validate cell type
  if (!"cell_type" %in% colnames(seurat_object@meta.data)) {
    stop("'cell_type' column not found in the Seurat metadata.")
  }
  if (!cell_type %in% seurat_object@meta.data$cell_type) {
    stop("The selected cell type is not present in the metadata.")
  }
  
  # Subset cells of the selected type
  cells_in_type <- which(seurat_object@meta.data$cell_type == cell_type)
  if (length(cells_in_type) == 0) {
    stop("No cells found for the selected cell type.")
  }
  
  # Extract normalized data (dense matrix may be necessary here for rowMeans)
  normalized_data <- as.matrix(seurat_object[["RNA"]]@data)
  
  # Fetch gene lengths using biomaRt
  gene_lengths <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "transcript_length"),
    filters = "ensembl_gene_id",
    values = rownames(normalized_data),
    mart = mart
  )
  
  # Ensure gene lengths align with expression data
  valid_genes <- intersect(rownames(normalized_data), gene_lengths$ensembl_gene_id)
  gene_lengths <- gene_lengths[gene_lengths$ensembl_gene_id %in% valid_genes, ]
  normalized_data <- normalized_data[valid_genes, ]
  
  # Convert lengths to kilobases
  gene_lengths <- setNames(gene_lengths$transcript_length / 1000, gene_lengths$ensembl_gene_id)
  
  # Calculate TPM
  tpm_matrix <- t(t(normalized_data) / gene_lengths[rownames(normalized_data)]) * 1e6
  
  # Calculate average expression for the cell type
  avg_expression <- rowMeans(tpm_matrix[, cells_in_type, drop = FALSE])
  
  # Identify the top N genes
  top_ensembl_ids <- names(sort(avg_expression, decreasing = TRUE))[1:n]
  
  # Map Ensembl IDs to gene names
  gene_mapping <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = top_ensembl_ids,
    mart = mart
  )
  
  # Remove genes with missing HGNC symbols
  top_gene_names <- gene_mapping$hgnc_symbol[gene_mapping$hgnc_symbol != ""]
  
  message("Top ", n, " genes expressed in cell type '", cell_type, "':")
  print(top_gene_names)
  
  return(top_gene_names)
}
