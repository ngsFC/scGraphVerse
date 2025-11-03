#' Compute Network Topological Properties
#'
#' Calculates topological properties (modularity, number of communities,
#' density, transitivity) for one or more networks.
#'
#' @param control_output A list output from `community_path()` representing the
#'   ground truth network. Must contain a `graph` (igraph object) and
#'   `communities$membership`.
#' @param predicted_list A list of lists, each output from `community_path()`
#'   representing predicted networks.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{topology_measures}: A data frame with Modularity,Communities,
#'       Density, and Transitivity for each prediction.
#'     \item \code{control_topology}: A numeric vector of the same metrics
#'       for the control network.
#'   }
#'
#' @details This function requires the \strong{igraph} package. Topological
#'   metrics include:
#'   \itemize{
#'     \item Modularity: Quality of community division
#'     \item Communities: Number of detected communities
#'     \item Density: Proportion of actual vs. possible edges
#'     \item Transitivity: Clustering coefficient
#'   }
#'
#' @importFrom igraph modularity edge_density transitivity is_igraph
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
#' topo <- compute_topology_metrics(comm_truth, list(comm_cons))
compute_topology_metrics <- function(control_output, predicted_list) {
    control_comm <- control_output$communities$membership
    control_graph <- control_output$graph
    control_topo <- .compute_topo_metrics(control_graph, control_comm)

    topology_comparison <- list()

    for (i in seq_along(predicted_list)) {
        pred <- predicted_list[[i]]
        pred_comm <- pred$communities$membership
        pred_graph <- pred$graph

        if (is.null(pred_graph) || !igraph::is_igraph(pred_graph)) {
            warning("Prediction ", i, " has no valid graph. Skipping topology.")
            topology_comparison[[paste0("Predicted_", i)]] <- rep(NA, 4)
            next
        }

        topology_comparison[[paste0("Predicted_", i)]] <-
            .compute_topo_metrics(pred_graph, pred_comm)
    }

    topo_df <- as.data.frame(do.call(rbind, topology_comparison))
    colnames(topo_df) <- c(
        "Modularity",
        "Communities",
        "Density",
        "Transitivity"
    )

    list(
        topology_measures = topo_df,
        control_topology  = control_topo
    )
}
