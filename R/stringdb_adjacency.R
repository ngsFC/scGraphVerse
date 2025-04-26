#' Build Adjacency Matrices for Physical Interactions from STRING (POST API)
#'
#' Constructs weighted and binary adjacency matrices for physical protein-protein interactions
#' using a POST request to the STRING database API. 
#'
#' @param genes A character vector of gene symbols or identifiers (e.g., \code{c("TP53", "BRCA1", ...)}).
#' @param species Integer. NCBI taxonomy ID of the species. Default is \code{9606} (human).
#' @param required_score Integer between 0 and 1000. Minimum confidence score for interactions. Default is \code{400}.
#' @param keep_all_genes Logical. If \code{TRUE} (default), includes all input genes in the final matrix even if unmapped.
#' @param verbose Logical. If \code{TRUE}, displays progress messages. Default is \code{TRUE}.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{weighted}: A square numeric adjacency matrix with interaction scores as weights.
#'   \item \code{binary}: A corresponding binary (0/1) adjacency matrix.
#' }
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Maps the input genes to STRING internal IDs.
#'   \item Uses a POST request to retrieve physical protein-protein interactions from STRING.
#'   \item Builds a weighted adjacency matrix using the STRING combined score.
#'   \item Builds a binary adjacency matrix indicating the presence or absence of interactions.
#' }
#'
#' Genes that could not be mapped to STRING are optionally retained as rows/columns of zeros if \code{keep_all_genes = TRUE}.
#'
#' @note
#' Requires the following packages: \pkg{STRINGdb}, \pkg{httr}, and \pkg{jsonlite}.
#'
#' @importFrom STRINGdb STRINGdb
#' @importFrom httr POST content
#' @importFrom jsonlite fromJSON
#' @export
#'
#' @examples
#' \dontrun{
#' # Define a large set of genes
#' large_gene_set <- c("TP53", "BRCA1", "MYC", "EGFR", "PTEN")
#'
#' # Retrieve adjacency matrices from STRING
#' adjacency_list <- stringdb_adjacency(
#'   genes = large_gene_set,
#'   species = 9606,
#'   required_score = 700
#' )
#'
#' # Access the weighted adjacency matrix
#' weighted_mat <- adjacency_list$weighted
#'
#' # Access the binary adjacency matrix
#' binary_mat <- adjacency_list$binary
#' }

stringdb_adjacency <- function(
    genes,
    species              = 9606,
    required_score       = 400,
    keep_all_genes       = TRUE, 
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
  mapped_genes <- mapped_genes[!is.na(mapped_genes$STRING_id), ]
  unmapped_genes <- setdiff(genes, mapped_genes$genes)
  
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
  identifiers_str <- paste(mapped_genes$STRING_id, collapse = "\n")
  
  res <- httr::POST(
    url = base_url,
    body = list(
      identifiers    = identifiers_str,
      species        = species,
      required_score = required_score,
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
    final_gene_list <- genes
  } else {
    final_gene_list <- unique(c(interactions$gene_A, interactions$gene_B)) 
  }
  
  # Initialize p×p matrices (filled with 0s)
  p <- length(final_gene_list)
  weighted_mat <- matrix(0, nrow = p, ncol = p, dimnames = list(final_gene_list, final_gene_list))
  
  # Populate adjacency matrix with STRING interaction data
  if (nrow(interactions) > 0) {
    for (i in seq_len(nrow(interactions))) {
      a <- interactions$gene_A[i]
      b <- interactions$gene_B[i]
      s <- interactions$interaction_score[i]  
      
      if (!is.na(a) && !is.na(b) && a %in% final_gene_list && b %in% final_gene_list) {
        weighted_mat[a, b] <- s
        weighted_mat[b, a] <- s
      }
    }
  }
  
  binary_mat <- ifelse(weighted_mat > 0, 1, 0)
  
  if (verbose) message("Adjacency matrices constructed successfully.")
  
  list(weighted = weighted_mat, binary = binary_mat)
}

