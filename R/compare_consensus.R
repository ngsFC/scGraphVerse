#' Compare Consensus and Reference Graphs or STRINGdb Networks
#'
#' Compares a consensus adjacency matrix to a reference network, either
#' provided manually or generated from STRINGdb. Visualizes True Positives
#' (TP), False Negatives (FN), and optionally False Positives (FP) edges.
#'
#' @param consensus_matrix A binary square adjacency matrix representing the
#'   consensus network. Row and column names should correspond to gene
#'   symbols.
#' @param reference_matrix Optional. A binary square adjacency matrix
#'   representing the reference (ground truth) network. If \code{NULL}, a
#'   STRINGdb high-confidence physical interaction network (human, score >
#'   900) is used.
#' @param false_plot Logical. If \code{TRUE}, an additional plot of False
#'   Positives (FP) is generated. Default is \code{FALSE}.
#'
#' @return A \code{ggplot} object visualizing the comparison. If
#'   \code{false_plot = TRUE}, a combined plot of True Positives / False
#'   Negatives and False Positives is returned.
#'
#' @details If no \code{reference_matrix} is provided, the function
#'   automatically queries STRINGdb to generate a high-confidence physical
#'   interaction network.
#'
#'   The plots differentiate:
#'     \itemize{
#'       \item Confirmed Edges (TP or CE): Present in both consensus and
#'         reference.
#'       \item Missing Edges (FN or ME): Present in reference but absent in
#'         consensus.
#'       \item Extra Edges (FP or EE): Present in consensus but absent in
#'         reference (only if \code{false_plot = TRUE}).
#'     }
#'
#' @note Requires the \pkg{igraph}, \pkg{ggraph}, \pkg{patchwork},
#'   \pkg{Matrix}, and \pkg{STRINGdb} packages.
#'
#' @importFrom igraph degree V graph_from_adjacency_matrix as_edgelist
#'   graph_from_edgelist delete_vertices
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @importFrom patchwork wrap_plots
#' @importFrom Matrix Matrix
#' @importFrom STRINGdb STRINGdb
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' # Simulate small example matrices
#' original <- matrix(
#'     sample(0:1, 25, replace = TRUE, prob = c(0.8, 0.2)),
#'     5, 5
#' )
#' consensus <- matrix(
#'     sample(0:1, 25, replace = TRUE, prob = c(0.8, 0.2)),
#'     5, 5
#' )
#' diag(original) <- diag(consensus) <- 0
#' rownames(original) <- colnames(original) <- paste0("Gene", 1:5)
#' rownames(consensus) <- colnames(consensus) <- paste0("Gene", 1:5)
#'
#' # Compare consensus network to original network
#' compare_consensus(
#'     consensus,
#'     reference_matrix = original,
#'     false_plot       = TRUE
#' )
compare_consensus <- function(
    consensus_matrix,
    reference_matrix = NULL,
    false_plot = FALSE) {
    if (!is.matrix(consensus_matrix)) {
        stop("consensus_matrix must be a binary adjacency matrix.")
    }

    use_STRINGdb <- is.null(reference_matrix)

    if (use_STRINGdb) {
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
    }

    if (!is.matrix(reference_matrix)) {
        stop("reference_matrix must be a binary adjacency matrix.")
    }
    if (!identical(
        dim(consensus_matrix),
        dim(reference_matrix)
    )) {
        stop("Matrices must have the same dimensions.")
    }

    graph_ref <- igraph::graph_from_adjacency_matrix(
        reference_matrix,
        mode = "undirected",
        diag = FALSE
    )
    graph_cons <- igraph::graph_from_adjacency_matrix(
        consensus_matrix,
        mode = "undirected",
        diag = FALSE
    )

    ref_edges <- .edge_to_str(igraph::as_edgelist(graph_ref))
    cons_edges <- .edge_to_str(igraph::as_edgelist(graph_cons))

    if (use_STRINGdb) {
        TP_label <- "CE (Confirmed Edges)"
        FN_label <- "ME (Missing Edges)"
        FP_label <- "EE (Extra Edges)"
    } else {
        TP_label <- "TP (True Positives)"
        FN_label <- "FN (False Negatives)"
        FP_label <- "FP (False Positives)"
    }

    edge_colors <- ifelse(
        ref_edges %in% cons_edges,
        "red",
        "blue"
    )
    plot_tp_fn <- .plot_tp_fn_graph(
        graph_ref,
        edge_colors,
        TP_label,
        FN_label
    )

    if (!false_plot) {
        return(plot_tp_fn)
    }

    fp_edges <- setdiff(cons_edges, ref_edges)
    if (length(fp_edges) == 0) {
        return(plot_tp_fn)
    }

    plot_fp <- .plot_fp_graph(fp_edges, FP_label)
    return(patchwork::wrap_plots(
        plot_tp_fn,
        plot_fp,
        nrow = 1
    ))
}
