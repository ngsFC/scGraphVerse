#' Build adjacency matrices (weighted & binary) for physical interactions from STRING using POST
#'
#' This function:
#' 1) Uses the STRING API with a POST request (network_type=physical).
#' 2) Accepts a large set of genes without triggering "URI Too Large" (414) errors.
#' 3) Allows choosing which score column (e.g., "score", "escore", etc.) to use as weights.
#' 4) Returns both a weighted adjacency matrix (numerical) and a binary adjacency matrix (0/1).
#' 5) Renames row/column names from internal STRING IDs to the "preferredName" in the returned data.
#'
#' @param genes A character vector of gene symbols/identifiers (e.g. c("TP53", "BRCA1", ...)).
#' @param species NCBI taxon ID. Default is 9606 (human).
#' @param required_score Minimum confidence (0–1000). Default 400.
#' @param score_col Name of the column to use as edge weights. 
#'        Commonly "score" (combined) or "escore" (experimental).
#' @param remove_missing_score Drop rows where the score is missing/"-"? Default TRUE.
#' @param verbose If TRUE, prints messages. Default TRUE.
#'
#' @return A list with two matrices: \code{$weighted} (numeric) and \code{$binary} (0/1).
#'
#' @examples
#' \dontrun{
#' large_gene_set <- c("TP53","BRCA1", "MYC", ...) # possibly hundreds/thousands of genes
#' adjacency_list <- stringdb_adjacency(
#'   genes          = large_gene_set,
#'   species        = 9606,
#'   required_score = 700,
#'   score_col      = "escore"
#' )
#' }
stringdb_adjacency <- function(
    genes,
    species              = 9606,
    required_score       = 400,
    keep_all_genes       = TRUE,  # Include all genes or only STRING-mapped ones
    verbose              = TRUE
) {
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop("Package 'STRINGdb' is required. Please install it via Bioconductor.")
  }
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required. Please install it.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required. Please install it.")
  }
  
  if (length(genes) == 0) {
    stop("Please provide at least one gene in 'genes'.")
  }
  
  if (verbose) message("Initializing STRINGdb...")
  
  # Initialize STRINGdb
  string_db <- STRINGdb$new(
    version = "11.5",
    species = species,
    score_threshold = required_score,
    input_directory = ""
  )
  
  # Map gene symbols to STRING IDs
  if (verbose) message("Mapping genes to STRING IDs...")
  mapped_genes <- string_db$map(
    data.frame(genes, stringsAsFactors = FALSE), "genes", removeUnmappedRows = FALSE
  )
  
  # Extract mapped and unmapped genes
  mapped_genes <- mapped_genes[!is.na(mapped_genes$STRING_id), ]  # Keep only mapped genes
  unmapped_genes <- setdiff(genes, mapped_genes$genes)  # Genes not found in STRING
  
  if (verbose) {
    message("Mapped ", nrow(mapped_genes), " genes to STRING IDs.")
    if (length(unmapped_genes) > 0 && keep_all_genes) {
      message(length(unmapped_genes), " genes were not found in STRING but will be included as zero rows/columns.")
    }
  }
  
  if (nrow(mapped_genes) == 0) {
    stop("No valid STRING IDs found for the provided genes.")
  }
  
  # Prepare API request using STRING IDs
  if (verbose) message("Retrieving **physical** interactions from STRING API...")
  
  base_url <- "https://string-db.org/api/json/network"
  identifiers_str <- paste(mapped_genes$STRING_id, collapse = "\n")  # Pass STRING IDs, NOT gene names
  
  res <- httr::POST(
    url = base_url,
    body = list(
      identifiers    = identifiers_str,
      species        = species,
      required_score = required_score,  # API already filters based on score
      network_type   = "physical"
    ),
    encode = "form"
  )
  
  if (res$status_code != 200) {
    stop("STRING API query failed. Status code: ", res$status_code)
  }
  
  # Parse JSON response
  interactions <- jsonlite::fromJSON(httr::content(res, "text", encoding = "UTF-8"))
  
  if (!is.data.frame(interactions) || nrow(interactions) == 0) {
    if (verbose) message("No STRING physical interactions found.")
    return(list(weighted = matrix(0, length(genes), length(genes), dimnames = list(genes, genes)),
                binary = matrix(0, length(genes), length(genes), dimnames = list(genes, genes))))
  }
  
  if (verbose) message("Found ", nrow(interactions), " STRING physical interactions.")
  
  # Ensure we use "score" instead of "combined_score"
  interactions$interaction_score <- interactions$score  # Rename for clarity
  
  # Map STRING IDs to Gene Names
  id_to_gene <- setNames(mapped_genes$genes, mapped_genes$STRING_id)
  
  # Convert STRING IDs in interaction table to Gene Names
  interactions$gene_A <- id_to_gene[interactions$stringId_A]
  interactions$gene_B <- id_to_gene[interactions$stringId_B]
  
  # Remove interactions where gene names couldn't be mapped
  interactions <- interactions[!is.na(interactions$gene_A) & !is.na(interactions$gene_B), ]
  
  # Select genes to include in the adjacency matrix
  if (keep_all_genes) {
    final_gene_list <- genes  # Keep all original genes
  } else {
    final_gene_list <- unique(c(interactions$gene_A, interactions$gene_B))  # Only genes in STRING interactions
  }
  
  # Initialize p×p matrices (filled with 0s)
  p <- length(final_gene_list)
  weighted_mat <- matrix(0, nrow = p, ncol = p, dimnames = list(final_gene_list, final_gene_list))
  
  # Populate adjacency matrix with STRING interaction data
  if (nrow(interactions) > 0) {
    for (i in seq_len(nrow(interactions))) {
      a <- interactions$gene_A[i]
      b <- interactions$gene_B[i]
      s <- interactions$interaction_score[i]  # Use correct "score" column
      
      if (!is.na(a) && !is.na(b) && a %in% final_gene_list && b %in% final_gene_list) {
        weighted_mat[a, b] <- s
        weighted_mat[b, a] <- s
      }
    }
  }
  
  # Create binary adjacency matrix (0/1)
  binary_mat <- ifelse(weighted_mat > 0, 1, 0)
  
  if (verbose) message("Adjacency matrices constructed successfully.")
  
  # Return adjacency matrices
  list(weighted = weighted_mat, binary = binary_mat)
}

