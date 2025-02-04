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
    species            = 9606,
    required_score     = 400,
    score_col          = "escore",
    remove_missing_score = TRUE,
    verbose           = TRUE
) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required. Please install it.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required. Please install it.")
  }
  library(httr)
  library(jsonlite)
  
  if (length(genes) == 0) {
    stop("Please provide at least one gene in 'genes'.")
  }
  
  if (verbose) message("Querying STRING for physical interactions via POST...")
  
  # String's JSON API endpoint
  base_url <- "https://string-db.org/api/json/network"
  
  # The STRING API can accept identifiers via POST.  
  # We join them with newlines (the API also accepts %0A or %0D, but newline is simpler in POST).
  identifiers_str <- paste(genes, collapse = "\n")
  
  # We'll pass these as form-encoded body parameters:
  body_list <- list(
    identifiers    = identifiers_str,
    species        = species,
    required_score = required_score,
    network_type   = "physical"
  )
  
  # Perform POST request
  res <- POST(url = base_url, body = body_list, encode = "form")
  if (res$status_code != 200) {
    stop("STRING query failed. Status code: ", res$status_code)
  }
  
  # Parse JSON content
  parsed_data <- fromJSON(content(res, "text", encoding = "UTF-8"))
  
  if (!is.data.frame(parsed_data) || nrow(parsed_data) == 0) {
    if (verbose) message("No STRING physical interactions found for these genes.")
    return(list(weighted = matrix(0,0,0), binary = matrix(0,0,0)))
  }
  
  if (verbose) {
    message("Found ", nrow(parsed_data), " STRING interaction records.")
  }
  
  # Typical columns: stringId_A, stringId_B, preferredName_A, preferredName_B,
  # score, escore, nscore, pscore, ascore, etc.
  if (!score_col %in% colnames(parsed_data)) {
    stop(
      "The chosen score_col '", score_col, "' is not in the data.\n",
      "Available: ", paste(colnames(parsed_data), collapse=", ")
    )
  }
  
  # Rename chosen column to "Score"
  parsed_data$Score <- parsed_data[[score_col]]
  
  # Remove missing/"-" scores if requested
  if (remove_missing_score) {
    parsed_data <- parsed_data[parsed_data$Score != "-", ]
  }
  parsed_data$Score <- suppressWarnings(as.numeric(parsed_data$Score))
  parsed_data <- parsed_data[!is.na(parsed_data$Score), ]
  
  if (!nrow(parsed_data)) {
    stop("No valid interaction rows remain after filtering out missing scores.")
  }
  
  # Build adjacency
  unique_ids <- unique(c(parsed_data$stringId_A, parsed_data$stringId_B))
  weighted_mat <- matrix(
    0,
    nrow = length(unique_ids),
    ncol = length(unique_ids),
    dimnames = list(unique_ids, unique_ids)
  )
  
  for (i in seq_len(nrow(parsed_data))) {
    a <- parsed_data$stringId_A[i]
    b <- parsed_data$stringId_B[i]
    s <- parsed_data$Score[i]
    weighted_mat[a, b] <- s
    weighted_mat[b, a] <- s
  }
  
  # Binary adjacency
  binary_mat <- ifelse(weighted_mat > 0, 1, 0)
  
  # Prune isolated nodes
  keep_idx <- which(rowSums(weighted_mat) > 0 & colSums(weighted_mat) > 0)
  if (length(keep_idx) == 0) {
    stop("No connected nodes remain after pruning. All edges have zero score?")
  }
  weighted_mat <- weighted_mat[keep_idx, keep_idx, drop = FALSE]
  binary_mat   <- binary_mat[keep_idx, keep_idx, drop = FALSE]
  
  # OPTIONAL: rename row/col from stringId -> gene name using "preferredName"
  # Build a map from each stringId to the corresponding "preferredName"
  mapA <- setNames(parsed_data$preferredName_A, parsed_data$stringId_A)
  mapB <- setNames(parsed_data$preferredName_B, parsed_data$stringId_B)
  id_map <- c(mapA, mapB)
  id_map <- id_map[!duplicated(names(id_map))]  # remove duplicates
  # Limit map to the IDs actually present
  id_map <- id_map[rownames(weighted_mat)]
  
  # Replace row/col names
  rownames(weighted_mat) <- id_map
  colnames(weighted_mat) <- id_map
  rownames(binary_mat)   <- id_map
  colnames(binary_mat)   <- id_map
  
  # Return
  list(weighted = weighted_mat, binary = binary_mat)
}
