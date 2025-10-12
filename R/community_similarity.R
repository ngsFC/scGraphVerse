#' Compare Community Assignments and Topological Properties
#'
#' Convenience wrapper that computes community assignment metrics,
#' topological properties, and optionally visualizes the comparison.
#' For more control, use the individual functions:
#' \code{\link{compute_community_metrics}},
#' \code{\link{compute_topology_metrics}}, and
#' \code{\link{plot_community_comparison}}.
#'
#' @param control_output A list output from `community_path()` representing the
#'   ground truth network. Must contain a `graph` (igraph object) and
#'   `communities$membership`.
#' @param predicted_list A list of lists, each output from `community_path()`
#'   representing predicted networks to compare.
#' @param plot Logical. If TRUE, displays plots immediately. If FALSE, no plots
#'   are displayed. Default: TRUE.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{community_metrics}: A data frame with VI, NMI, and ARI
#'       scores for each prediction.
#'     \item \code{topology_measures}: A data frame with raw topological
#'       metrics for each prediction.
#'     \item \code{control_topology}: A list of raw topological metrics for
#'       the ground truth network.
#'   }
#'
#' @details This function requires the \strong{igraph} package. If
#'   \code{plot = TRUE}, the \strong{fmsb} package is also required.
#'   Community similarity is measured using variation of information (VI),
#'   normalized mutual information (NMI), and adjusted Rand index (ARI).
#'
#' @seealso \code{\link{compute_community_metrics}},
#'   \code{\link{compute_topology_metrics}},
#'   \code{\link{plot_community_comparison}}
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
#' comm_cons <- community_path(consensus)
#' comm_truth <- community_path(toy_adj_matrix)
#'
#' sim_score <- community_similarity(comm_truth, list(comm_cons))
community_similarity <- function(
    control_output,
    predicted_list,
    plot = TRUE) {
    # Compute community metrics
    comm_df <- compute_community_metrics(control_output, predicted_list)

    # Compute topology metrics
    topo_result <- compute_topology_metrics(control_output, predicted_list)

    # Display plots if requested
    if (plot) {
        plot_community_comparison(
            comm_df,
            topo_result$topology_measures,
            topo_result$control_topology
        )
    }

    list(
        community_metrics = comm_df,
        topology_measures = topo_result$topology_measures,
        control_topology  = topo_result$control_topology
    )
}
