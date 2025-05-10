#' Edge Mining of Gene Interactions Using PubMed
#'
#' Query PubMed for Literature Evidence Supporting Predicted Gene-Gene Interactions
#'
#' This function compares predicted gene adjacency matrices against a ground truth
#' matrix, identifies specific edge types (True Positives, False Positives, False Negatives),
#' and queries PubMed for each gene pair. It returns the number of PubMed hits, associated PMIDs, and query status.
#'
#' @param predicted_list A list of predicted gene adjacency matrices (with row and column names as gene symbols).
#' @param ground_truth A ground truth adjacency matrix (0/1) with row and column names.
#' @param delay Numeric. Seconds to wait between consecutive queries (default = 0.5).
#' @param query_field Character. The PubMed search field to use. Options: "Title/Abstract" (default), "Title", "Abstract".
#' @param query_edge_types Character vector. Which edge types to query: "TP", "FP", "FN" (default all).
#' @param max_retries Integer. Number of maximum retries for PubMed queries (default = 3).
#' @param BPPARAM A BiocParallel parameter object specifying how parallelization should be performed. Default = bpparam().
#'
#' @return A named list of data frames. Each data frame includes:
#' \describe{
#'   \item{gene1}{First gene in the interaction}
#'   \item{gene2}{Second gene}
#'   \item{edge_type}{"TP", "FP", or "FN"}
#'   \item{pubmed_hits}{Number of PubMed hits}
#'   \item{PMIDs}{Comma-separated PubMed IDs or NA}
#'   \item{query_status}{"hits_found", "no_hits", "error"}
#' }
#'
#' @import BiocParallel rentrez
#' @export
#'
#' @examples
#'
#' set.seed(123)
#' predicted <- matrix(rbinom(100, 1, 0.1), nrow = 10)
#' rownames(predicted) <- colnames(predicted) <- paste0("Gene", 1:10)
#' ground_truth <- matrix(rbinom(100, 1, 0.05), nrow = 10)
#' rownames(ground_truth) <- colnames(ground_truth) <- paste0("Gene", 1:10)
#'
#' results <- edge_mining(predicted_list = list(pred_net = predicted), ground_truth = ground_truth)
#' head(results$pred_net)
edge_mining <- function(predicted_list, ground_truth, delay = 1, query_field = "Title/Abstract",
                        query_edge_types = c("TP", "FP", "FN"), max_retries = 10,
                        BPPARAM = BiocParallel::bpparam()) {
  if (!is.list(predicted_list) || !is.matrix(ground_truth)) {
    stop("predicted_list must be a list of matrices and ground_truth must be a matrix")
  }

  field_map <- c("Title/Abstract" = "TIAB", "Title" = "TI", "Abstract" = "AB")
  if (query_field %in% names(field_map)) query_field <- field_map[[query_field]]

  results_list <- BiocParallel::bplapply(seq_along(predicted_list), function(i) {
    predicted <- predicted_list[[i]]
    if (!is.matrix(predicted) || is.null(rownames(predicted)) || is.null(colnames(predicted))) {
      stop(sprintf("Predicted matrix at index %d lacks row/column names.", i))
    }

    gene_pairs <- .identify_edges(predicted, ground_truth, query_edge_types)
    if (is.null(gene_pairs)) {
      return(data.frame(
        gene1 = character(0), gene2 = character(0), edge_type = character(0),
        pubmed_hits = integer(0), PMIDs = character(0), query_status = character(0)
      ))
    }

    .query_edge_pairs(gene_pairs, query_field, delay, max_retries, BPPARAM)
  }, BPPARAM = BPPARAM)

  if (!is.null(names(predicted_list))) {
    names(results_list) <- names(predicted_list)
  }

  return(results_list)
}
