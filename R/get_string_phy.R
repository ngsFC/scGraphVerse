#' Get physical interactions from STRING for a set of genes
#'
#' This function queries the STRING API for physical (PPI) interactions
#' using the `network_type=physical` parameter.
#'
#' @param genes Character vector of gene names or identifiers (e.g. HGNC or Ensembl)
#' @param species NCBI Taxonomy ID (default 9606 for human)
#' @param required_score Minimum confidence score (0–1000, defaults to 400)
#' @param verbose Logical indicating whether to message progress
#' @return A data.frame of STRING interactions
#'
#' @examples
#' \dontrun{
#'   my_string_interactions <- get_physical_interactions_string(
#'     genes           = c("TP53", "BRCA1"),
#'     species         = 9606,
#'     required_score  = 400
#'   )
#' }
#'
get_string_phys <- function(
    genes,
    species        = 9606,
    required_score = 400,
    verbose        = TRUE
) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("The 'httr' package is required for this function. Please install it.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The 'jsonlite' package is required for this function. Please install it.")
  }
  
  if (length(genes) < 1) {
    stop("Please provide at least one gene name.")
  }
  
  # The STRING API expects identifiers joined by "%0d" or "%0A" (line break).
  gene_list_str <- paste(genes, collapse = "%0D")
  
  # Build URL. For details, see:
  # https://string-db.org/cgi/help?subpage=api
  base_url <- "https://string-db.org/api/json/network"
  
  
  params <- list(
    identifiers     = gene_list_str,
    species         = species,
    required_score  = required_score,
    network_type    = "physical"
  )
  
  if (verbose) message("Querying STRING for physical interactions...")
  
  res <- GET(base_url, query = params)
  if (res$status_code != 200) {
    stop("STRING query failed. Status code: ", res$status_code)
  }
  
  # Parse JSON
  parsed_data <- fromJSON(content(res, "text", encoding = "UTF-8"))
  
  # If the query returns no data or an error message
  if (!is.data.frame(parsed_data) || nrow(parsed_data) == 0) {
    if (verbose) message("No STRING physical interactions found.")
    return(data.frame())
  }
  
  if (verbose) {
    message("Found ", nrow(parsed_data), " STRING interaction records.")
  }
  
  parsed_data
}