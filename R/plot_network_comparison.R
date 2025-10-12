#' Visualize Network Comparison
#'
#' Creates visualization plots comparing consensus and reference networks,
#' showing True Positives (TP), False Negatives (FN), and optionally
#' False Positives (FP) edges.
#'
#' @param edge_classification A list output from `classify_edges()`
#'   containing edge classifications and graph objects.
#' @param show_fp Logical. If TRUE, displays False Positive edges in a
#'   separate plot. Default: FALSE.
#'
#' @return A \code{ggplot} object visualizing the comparison. If
#'   \code{show_fp = TRUE}, a combined plot using patchwork is returned.
#'
#' @details This function requires the \strong{ggraph} and \strong{ggplot2}
#'   packages. If \code{show_fp = TRUE}, the \strong{patchwork} package is
#'   also required.
#'
#'   The plots differentiate:
#'   \itemize{
#'     \item TP/CE (True Positives/Confirmed Edges): Red edges present in both
#'     \item FN/ME (False Negatives/Missing Edges): Blue edges in reference only
#'     \item FP/EE (False Positives/Extra Edges): Edges in consensus only
#'   }
#'
#'   If STRINGdb was used for reference, labels are CE/ME/EE. Otherwise,
#'   TP/FN/FP.
#'
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
#' # Classify edges
#' edge_class <- classify_edges(consensus, ref_se)
#'
#' # Plot comparison
#' plot_network_comparison(edge_class, show_fp = FALSE)
plot_network_comparison <- function(edge_classification, show_fp = FALSE) {
    BiocBaseUtils::checkInstalled("ggraph")
    BiocBaseUtils::checkInstalled("ggplot2")
    if (show_fp) {
        BiocBaseUtils::checkInstalled("patchwork")
    }

    # Determine labels based on whether STRINGdb was used
    if (edge_classification$use_stringdb) {
        TP_label <- "CE (Confirmed Edges)"
        FN_label <- "ME (Missing Edges)"
        FP_label <- "EE (Extra Edges)"
    } else {
        TP_label <- "TP (True Positives)"
        FN_label <- "FN (False Negatives)"
        FP_label <- "FP (False Positives)"
    }

    # Create TP/FN plot
    plot_tp_fn <- .plot_tp_fn_graph(
        edge_classification$reference_graph,
        edge_classification$edge_colors,
        TP_label,
        FN_label
    )

    if (!show_fp) {
        return(plot_tp_fn)
    }

    # Create FP plot if requested
    fp_edges <- edge_classification$fp_edges
    if (length(fp_edges) == 0) {
        return(plot_tp_fn)
    }

    plot_fp <- .plot_fp_graph(fp_edges, FP_label)
    return(patchwork::wrap_plots(plot_tp_fn, plot_fp, nrow = 1))
}
