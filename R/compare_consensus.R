#' Compare Consensus and Reference Graphs or STRINGdb Networks
#'
#' Convenience wrapper that classifies edges and visualizes the comparison
#' between consensus and reference networks. For more control, use the
#' individual functions: \code{\link{classify_edges}} and
#' \code{\link{plot_network_comparison}}.
#'
#' @param consensus_matrix A \linkS4class{SummarizedExperiment} object
#'   representing the consensus network.
#' @param reference_matrix Optional. A \linkS4class{SummarizedExperiment} obj
#'   representing the reference network. If \code{NULL}, STRINGdb is queried.
#' @param false_plot Logical. If \code{TRUE}, displays False Positives plot.
#'   Default is \code{FALSE}.
#'
#' @return A \code{ggplot} object visualizing the comparison.
#'
#' @details If no \code{reference_matrix} is provided, STRINGdb is queried
#'   to generate a high-confidence physical interaction network.
#'
#' @seealso \code{\link{classify_edges}}, \code{\link{plot_network_comparison}}
#'
#' @note Requires \pkg{ggraph} and \pkg{ggplot2}. If \code{reference_matrix}
#'   is NULL, also requires \pkg{STRINGdb}. If \code{false_plot = TRUE},
#'   requires \pkg{patchwork}.
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
#' head(networks[[1]])
#'
#' # Generate adjacency matrices
#' wadj_se <- generate_adjacency(networks)
#' swadj_se <- symmetrize(wadj_se, weight_function = "mean")
#'
#' # Apply cutoff
#' binary_se <- cutoff_adjacency(
#'     count_matrices = toy_counts,
#'     weighted_adjm_list = swadj_se,
#'     n = 1,
#'     method = "GENIE3",
#'     quantile_threshold = 0.95,
#'     nCores = 1,
#'     debug = TRUE
#' )
#' head(binary_se[[1]])
#'
#' consensus <- create_consensus(binary_se, method = "union")
#'
#' # Wrap reference matrix in SummarizedExperiment
#' ref_se <- build_network_se(list(reference = toy_adj_matrix))
#'
#' # Compare consensus to reference
#' compare_consensus(
#'     consensus,
#'     reference_matrix = ref_se,
#'     false_plot = FALSE
#' )
#'
compare_consensus <- function(
    consensus_matrix,
    reference_matrix = NULL,
    false_plot = FALSE) {
    # Classify edges
    edge_class <- classify_edges(
        consensus_matrix,
        reference_matrix,
        use_stringdb = is.null(reference_matrix)
    )

    # Create visualization
    plot_network_comparison(edge_class, show_fp = false_plot)
}
