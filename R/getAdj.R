# Function to create STRING adjacency matrices
get_string_adjacency <- function(top_genes, mart, score_threshold = 900, excluded_interactions = c("textmining")) {
  # Load the STRING database
  string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = score_threshold)
  
  # Map top genes to STRING IDs
  mapped_genes <- string_db$map(data.frame(gene = top_genes), "gene", removeUnmappedRows = TRUE)
  if (nrow(mapped_genes) == 0) {
    stop("No genes could be mapped to STRING IDs.")
  }
  
  # Retrieve interactions for mapped genes
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)
  
  # Filter interactions by type (exclude unwanted interaction types)
  if (!is.null(excluded_interactions)) {
    interactions <- interactions[!interactions$evidence %in% excluded_interactions, ]
  }
  
  # Filter interactions based on the score threshold
  interactions <- interactions[interactions$combined_score >= score_threshold, ]
  
  if (nrow(interactions) == 0) {
    stop("No interactions remain after filtering by type and score threshold.")
  }
  
  # Construct weighted adjacency matrix
  unique_proteins <- unique(c(interactions$from, interactions$to))
  weighted_adj_matrix <- matrix(0, nrow = length(unique_proteins), ncol = length(unique_proteins),
                                dimnames = list(unique_proteins, unique_proteins))
  for (i in seq_len(nrow(interactions))) {
    weighted_adj_matrix[interactions$from[i], interactions$to[i]] <- interactions$combined_score[i]
    weighted_adj_matrix[interactions$to[i], interactions$from[i]] <- interactions$combined_score[i]
  }
  
  # Create a binary adjacency matrix
  binary_adj_matrix <- ifelse(weighted_adj_matrix > 0, 1, 0)
  
  # Map STRING IDs (ENSP) to gene names
  ensp_ids <- rownames(weighted_adj_matrix)
  mapping <- biomaRt::getBM(
    attributes = c("ensembl_peptide_id", "external_gene_name"),
    filters = "ensembl_peptide_id",
    values = sub("9606\\.", "", ensp_ids),
    mart = mart
  )
  
  # Create a mapping dictionary
  id_to_gene <- setNames(mapping$external_gene_name, mapping$ensembl_peptide_id)
  
  # Replace STRING IDs with gene names in the adjacency matrices
  rownames(weighted_adj_matrix) <- id_to_gene[sub("9606\\.", "", rownames(weighted_adj_matrix))]
  colnames(weighted_adj_matrix) <- id_to_gene[sub("9606\\.", "", colnames(weighted_adj_matrix))]
  
  rownames(binary_adj_matrix) <- id_to_gene[sub("9606\\.", "", rownames(binary_adj_matrix))]
  colnames(binary_adj_matrix) <- id_to_gene[sub("9606\\.", "", colnames(binary_adj_matrix))]
  
  # Filter out rows and columns with missing gene names
  valid_genes <- !is.na(rownames(weighted_adj_matrix)) & !is.na(colnames(weighted_adj_matrix))
  weighted_adj_matrix <- weighted_adj_matrix[valid_genes, valid_genes]
  binary_adj_matrix <- binary_adj_matrix[valid_genes, valid_genes]
  
  # Return both matrices as a list
  return(list(weighted = weighted_adj_matrix, binary = binary_adj_matrix))
}

