topExp <- function(seurat_object, cell_type, mart, n = 100) {
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
  
  # Fetch gene lengths from Ensembl
  gene_lengths <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "transcript_length"),
    filters = "ensembl_gene_id",
    values = rownames(normalized_data),
    mart = mart
  )
  
  gene_lengths <- gene_lengths$transcript_length / 1000
  normalized_data <- normalized_data[rownames(normalized_data) %in% names(gene_lengths), ]
  gene_lengths <- gene_lengths[rownames(normalized_data)]
  
  tpm_matrix <- t(t(normalized_data) / gene_lengths) * 1e6
  avg_expression <- rowMeans(tpm_matrix[, cells_in_type, drop = FALSE])
  top_ensembl_ids <- names(sort(avg_expression, decreasing = TRUE))[1:n]
  
  # Map Ensembl IDs to gene names
  gene_mapping <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = top_ensembl_ids,
    mart = mart
  )
  
  top_gene_names <- gene_mapping$hgnc_symbol[gene_mapping$hgnc_symbol != ""]
  
  message("Top ", n, " genes expressed in cell type '", cell_type, "':")
  print(top_gene_names)
  
  return(top_gene_names)
}

