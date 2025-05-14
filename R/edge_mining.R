#' Edge Mining of Gene Interactions Using PubMed
#'
#' Query PubMed for literature evidence supporting predicted
#' gene–gene interactions.
#'
#' This function compares predicted adjacency matrices against a ground truth
#' matrix, identifies edge types (TP, FP, FN), and queries PubMed for each
#' gene pair. Returns counts of hits, PMIDs, and query status.
#'
#' @param predicted_list A list of predicted adjacency matrices (row and
#'   column names are gene symbols).
#' @param ground_truth A 0/1 adjacency matrix with row and column names.
#' @param delay Numeric. Seconds to wait between consecutive queries
#'   (default = 1).
#' @param query_field Character. PubMed search field. Options:
#'   "Title/Abstract" (default), "Title", "Abstract".
#' @param query_edge_types Character vector. Edge types to query:
#'   c("TP", "FP", "FN") (default all).
#' @param max_retries Integer. Max retries for PubMed queries
#'   (default = 10).
#' @param BPPARAM A BiocParallel parameter object. Default = bpparam().
#'
#' @return A named list of data.frames. Each data.frame has columns:
#'   \describe{
#'     \item{gene1}{First gene in interaction}
#'     \item{gene2}{Second gene}
#'     \item{edge_type}{One of "TP", "FP", or "FN"}
#'     \item{pubmed_hits}{Number of PubMed hits}
#'     \item{PMIDs}{Comma-separated PubMed IDs or NA}
#'     \item{query_status}{One of "hits_found", "no_hits", or "error"}
#'   }
#'
#' @import BiocParallel rentrez
#' @export
#'
#' @examples
#' data(count_matrices)
#' data(adj_truth)
#'
#' networks <- infer_networks(
#'     count_matrices_list = count_matrices,
#'     method = "GENIE3",
#'     nCores = 1
#' )
#' head(networks[[1]])
#'
#' wadj_list <- generate_adjacency(networks)
#' swadj_list <- symmetrize(wadj_list, weight_function = "mean")
#'
#' binary_listj <- cutoff_adjacency(
#'     count_matrices = count_matrices,
#'     weighted_adjm_list = swadj_list,
#'     n = 2,
#'     method = "GENIE3",
#'     quantile_threshold = 0.99,
#'     nCores = 1,
#'     debug = TRUE
#' )
#' head(binary_listj[[1]])
#'
#' consensus <- create_consensus(binary_listj, method = "vote")
#' head(consensus)
#' em <- edge_mining(list(consensus), adj_truth, query_edge_types = "TP")
edge_mining <- function(
    predicted_list,
    ground_truth,
    delay = 1,
    query_field = "Title/Abstract",
    query_edge_types = c("TP", "FP", "FN"),
    max_retries = 10,
    BPPARAM = BiocParallel::bpparam()) {
    stopifnot(
        is.list(predicted_list),
        is.matrix(ground_truth)
    )

    field_map <- c(
        "Title/Abstract" = "TIAB",
        "Title"          = "TI",
        "Abstract"       = "AB"
    )
    if (query_field %in% names(field_map)) {
        query_field <- field_map[[query_field]]
    }

    results_list <- BiocParallel::bplapply(
        seq_along(predicted_list),
        function(i) {
            predicted <- predicted_list[[i]]
            if (!is.matrix(predicted) ||
                is.null(rownames(predicted)) ||
                is.null(colnames(predicted))) {
                stop(sprintf(
                    "Predicted matrix at index %d lacks row/column names.", i
                ))
            }

            gene_pairs <- .identify_edges(
                predicted,
                ground_truth,
                query_edge_types
            )
            if (is.null(gene_pairs)) {
                return(data.frame(
                    gene1        = character(0),
                    gene2        = character(0),
                    edge_type    = character(0),
                    pubmed_hits  = integer(0),
                    PMIDs        = character(0),
                    query_status = character(0)
                ))
            }

            .query_edge_pairs(
                gene_pairs,
                query_field,
                delay,
                max_retries,
                BPPARAM
            )
        },
        BPPARAM = BPPARAM
    )

    if (!is.null(names(predicted_list))) {
        names(results_list) <- names(predicted_list)
    }

    return(results_list)
}
