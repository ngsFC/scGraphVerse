#' Visualize Community and Topology Comparison
#'
#' Creates visualization plots for community assignment metrics and
#' topological properties comparison.
#'
#' @param community_metrics A data frame with VI, NMI, and ARI scores
#'   (output from `compute_community_metrics()`).
#' @param topology_measures A data frame with Modularity, Communities,
#'   Density, and Transitivity (from `compute_topology_metrics()`).
#' @param control_topology A named numeric vector of control network
#'   topology metrics (from `compute_topology_metrics()`).
#'
#' @return Invisibly returns NULL. Displays a radar plot for community
#'   metrics and bar plots for topology comparison.
#'
#' @details This function requires the \strong{fmsb} package for radar
#'   chart visualization. The radar plot shows normalized community
#'   similarity metrics. The bar plots compare raw topological properties
#'   between predicted and control networks.
#'
#' @importFrom graphics barplot par legend
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
#'     nCores = 1
#' )
#'
#' consensus <- create_consensus(binary_se, method = "union")
#' comm_cons <- community_path(consensus)
#' comm_truth <- community_path(toy_adj_matrix)
#'
#' comm_metrics <- compute_community_metrics(comm_truth, list(comm_cons))
#' topo_metrics <- compute_topology_metrics(comm_truth, list(comm_cons))
#'
#' plot_community_comparison(
#'     comm_metrics,
#'     topo_metrics$topology_measures,
#'     topo_metrics$control_topology
#' )
plot_community_comparison <- function(
    community_metrics,
    topology_measures,
    control_topology) {
    BiocBaseUtils::checkInstalled("fmsb")

    # Plot radar chart for community metrics
    .plot_radar_communities(community_metrics)

    # Plot bar charts for topology
    .plot_topo_barplots(topology_measures, control_topology)

    invisible(NULL)
}
