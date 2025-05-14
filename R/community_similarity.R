#' Compare Community Assignments and Topological Properties
#'
#' Evaluates similarity between a ground truth community structure and one
#' or more predicted community structures. Computes community assignment
#' metrics (VI, NMI, ARI) and raw topological properties (Modularity,
#' Number of Communities, Density, Transitivity). Visualizes results via a
#' radar plot for community assignment and bar plots for topology.
#'
#' @param control_output A list output from `community_path()` representing the
#'   ground truth network. Must contain a `graph` (igraph object) and
#'   `communities$membership`.
#' @param predicted_list A list of lists, each output from `community_path()`
#'   representing predicted networks to compare.
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
#' @details This function requires the \strong{igraph} and \strong{fmsb}
#'   packages. Community similarity is measured using variation of
#'   information (VI), normalized mutual information (NMI), and adjusted
#'   Rand index (ARI). Topological properties are compared by directly
#'   plotting raw values without normalization.
#'
#' @importFrom igraph modularity edge_density transitivity compare is_igraph
#' @importFrom fmsb radarchart
#' @importFrom graphics barplot par legend
#' @export
#'
#' @examples
#' data(count_matrices)
#' data(adj_truth)
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
#' comm_cons <- community_path(consensus)
#' comm_truth <- community_path(adj_truth)
#'
#' sim_score <- community_similarity(comm_truth, list(comm_cons))
community_similarity <- function(
    control_output,
    predicted_list) {
    required_pkgs <- c("igraph", "fmsb")
    missing_pkgs <- required_pkgs[
        !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing_pkgs) > 0) {
        stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
    }

    control_comm <- control_output$communities$membership
    control_graph <- control_output$graph
    control_topo <- .compute_topo_metrics(control_graph, control_comm)

    community_metrics <- list()
    topology_comparison <- list()

    for (i in seq_along(predicted_list)) {
        pred <- predicted_list[[i]]
        pred_comm <- pred$communities$membership
        pred_graph <- pred$graph

        community_metrics[[paste0("Predicted_", i)]] <-
            .compare_communities(control_comm, pred_comm)

        if (is.null(pred_graph) || !igraph::is_igraph(pred_graph)) {
            warning("Prediction ", i, " has no valid graph. Skipping topology.")
            topology_comparison[[paste0("Predicted_", i)]] <- rep(NA, 4)
            next
        }

        topology_comparison[[paste0("Predicted_", i)]] <-
            .compute_topo_metrics(pred_graph, pred_comm)
    }

    comm_df <- as.data.frame(do.call(rbind, community_metrics))
    topo_df <- as.data.frame(do.call(rbind, topology_comparison))
    colnames(topo_df) <- c(
        "Modularity",
        "Communities",
        "Density",
        "Transitivity"
    )

    .plot_radar_communities(comm_df)
    .plot_topo_barplots(topo_df, control_topo)

    list(
        community_metrics = comm_df,
        topology_measures = topo_df,
        control_topology  = control_topo
    )
}
