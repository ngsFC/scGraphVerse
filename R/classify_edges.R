#' Classify Edges as TP, FP, or FN
#'
#' Compares a consensus network to a reference network and classifies
#' edges as True Positives (TP), False Positives (FP), or False
#' Negatives (FN).
#'
#' @param consensus_matrix A \linkS4class{SummarizedExperiment} object
#'   representing the consensus network.
#' @param reference_matrix A \linkS4class{SummarizedExperiment} object
#'   representing the reference (ground truth) network. If \code{NULL},
#'   STRINGdb is used to generate a high-confidence physical network.
#' @param use_stringdb Logical. If TRUE and \code{reference_matrix} is NULL,
#'   queries STRINGdb for reference network. Default: TRUE.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{tp_edges}: Character vector of True Positive edges
#'     \item \code{fp_edges}: Character vector of False Positive edges
#'     \item \code{fn_edges}: Character vector of False Negative edges
#'     \item \code{consensus_graph}: igraph object of consensus network
#'     \item \code{reference_graph}: igraph object of reference network
#'     \item \code{edge_colors}: Color vector for TP (red) and FN (blue) edges
#'     \item \code{use_stringdb}: Logical indicating if STRINGdb was used
#'   }
#'
#' @details If \code{reference_matrix} is NULL and \code{use_stringdb = TRUE},
#'   this function queries STRINGdb to generate a human high-confidence
#'   (score > 900) physical interaction network.
#'
#' @importFrom igraph graph_from_adjacency_matrix as_edgelist
#' @export
#'
#' @examples
#' data(toy_counts)
#' data(toy_adj_matrix)
#'
#'
#' # Infer networks (toy_counts is already a MultiAssayExperiment)
#' networks <- infer_networks(
#'     count_matrices_list = toy_counts,
#'     method = "GENIE3",
#'     nCores = 1
#' )
#'
#' # Generate and symmetrize adjacency matrices (returns SummarizedExperiment)
#' wadj_se <- generate_adjacency(networks)
#' swadj_se <- symmetrize(wadj_se, weight_function = "mean")
#'
#' # Apply cutoff (returns SummarizedExperiment)
#' binary_se <- cutoff_adjacency(
#'     count_matrices = toy_counts,
#'     weighted_adjm_list = swadj_se,
#'     n = 1,
#'     method = "GENIE3",
#'     quantile_threshold = 0.95,
#'     nCores = 1
#' )
#'
#' # Create consensus (returns SummarizedExperiment)
#' consensus <- create_consensus(binary_se, method = "union")
#'
#' # Wrap reference matrix in SummarizedExperiment
#' ref_se <- build_network_se(list(reference = toy_adj_matrix))
#'
#' # classify_edges expects SummarizedExperiment objects
#' edge_class <- classify_edges(consensus, ref_se)
classify_edges <- function(
    consensus_matrix,
    reference_matrix = NULL,
    use_stringdb = TRUE) {
    # Check input types
    if (!inherits(consensus_matrix, "SummarizedExperiment")) {
        stop("consensus_matrix must be a SummarizedExperiment object.")
    }

    if (!is.null(reference_matrix) && !inherits(
        reference_matrix,
        "SummarizedExperiment"
    )) {
        stop("reference_matrix must be a SummarizedExperiment object or NULL.")
    }

    # Extract matrices from SummarizedExperiment
    consensus_list <- .extract_networks_from_se(consensus_matrix)
    consensus_matrix <- consensus_list[[1]]

    if (!is.null(reference_matrix)) {
        reference_list <- .extract_networks_from_se(reference_matrix)
        reference_matrix <- reference_list[[1]]
    }

    # Generate STRINGdb reference if needed
    stringdb_used <- FALSE
    if (is.null(reference_matrix) && use_stringdb) {
        BiocBaseUtils::checkInstalled("STRINGdb")

        if (is.null(rownames(consensus_matrix))) {
            stop("consensus_matrix must have row names to query STRINGdb.")
        }

        adj <- stringdb_adjacency(
            genes          = rownames(consensus_matrix),
            species        = 9606,
            required_score = 900,
            keep_all_genes = TRUE
        )$binary

        reference_matrix <- symmetrize(
            list(
                adj[
                    rownames(consensus_matrix),
                    rownames(consensus_matrix)
                ]
            ),
            "mean"
        )[[1]]
        stringdb_used <- TRUE
    }

    # Validate dimensions match
    if (!identical(dim(consensus_matrix), dim(reference_matrix))) {
        stop("Matrices must have the same dimensions.")
    }

    # Create graphs
    graph_ref <- igraph::graph_from_adjacency_matrix(
        reference_matrix,
        mode = "max",
        diag = FALSE
    )
    graph_cons <- igraph::graph_from_adjacency_matrix(
        consensus_matrix,
        mode = "max",
        diag = FALSE
    )

    # Get edge lists
    ref_edges <- .edge_to_str(igraph::as_edgelist(graph_ref))
    cons_edges <- .edge_to_str(igraph::as_edgelist(graph_cons))

    # Classify edges
    tp_edges <- intersect(ref_edges, cons_edges)
    fp_edges <- setdiff(cons_edges, ref_edges)
    fn_edges <- setdiff(ref_edges, cons_edges)

    # Edge colors for reference graph (TP=red, FN=blue)
    edge_colors <- ifelse(ref_edges %in% cons_edges, "red", "blue")

    list(
        tp_edges = tp_edges,
        fp_edges = fp_edges,
        fn_edges = fn_edges,
        consensus_graph = graph_cons,
        reference_graph = graph_ref,
        edge_colors = edge_colors,
        use_stringdb = stringdb_used
    )
}
