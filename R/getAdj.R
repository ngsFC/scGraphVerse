getAdj <- function(top_genes, mart, score_threshold = 900) {
  string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = score_threshold)
  
  mapped_genes <- string_db$map(data.frame(gene = top_genes), "gene", removeUnmappedRows = TRUE)
  if (nrow(mapped_genes) == 0) {
    stop("No genes could be mapped to STRING IDs.")
  }
  
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)
  interactions <- interactions[interactions$combined_score >= score_threshold, ]
  
  if (nrow(interactions) == 0) {
    stop("No interactions remain after filtering by score threshold.")
  }
  
  unique_proteins <- unique(c(interactions$from, interactions$to))
  weighted_adj_matrix <- matrix(0, nrow = length(unique_proteins), ncol = length(unique_proteins),
                                dimnames = list(unique_proteins, unique_proteins))
  for (i in seq_len(nrow(interactions))) {
    weighted_adj_matrix[interactions$from[i], interactions$to[i]] <- interactions$combined_score[i]
    weighted_adj_matrix[interactions$to[i], interactions$from[i]] <- interactions$combined_score[i]
  }
  
  binary_adj_matrix <- ifelse(weighted_adj_matrix > 0, 1, 0)
  
  ensp_ids <- rownames(weighted_adj_matrix)
  mapping <- biomaRt::getBM(
    attributes = c("ensembl_peptide_id", "external_gene_name"),
    filters = "ensembl_peptide_id",
    values = sub("9606\\.", "", ensp_ids),
    mart = mart
  )
  
  id_to_gene <- setNames(mapping$external_gene_name, mapping$ensembl_peptide_id)
  
  rownames(weighted_adj_matrix) <- id_to_gene[sub("9606\\.", "", rownames(weighted_adj_matrix))]
  colnames(weighted_adj_matrix) <- id_to_gene[sub("9606\\.", "", colnames(weighted_adj_matrix))]
  rownames(binary_adj_matrix) <- id_to_gene[sub("9606\\.", "", rownames(binary_adj_matrix))]
  colnames(binary_adj_matrix) <- id_to_gene[sub("9606\\.", "", colnames(binary_adj_matrix))]
  
  valid_genes <- !is.na(rownames(weighted_adj_matrix)) & !is.na(colnames(weighted_adj_matrix))
  weighted_adj_matrix <- weighted_adj_matrix[valid_genes, valid_genes]
  binary_adj_matrix <- binary_adj_matrix[valid_genes, valid_genes]
  
  return(list(weighted = weighted_adj_matrix, binary = binary_adj_matrix))
}
