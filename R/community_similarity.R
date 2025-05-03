#' Compare Community Assignments and Topological Properties
#'
#' This function evaluates the similarity between a ground truth community structure and one or more predicted community structures.
#' It computes community assignment metrics (VI, NMI, ARI) and raw topological properties (Modularity, Number of Communities, Density, Transitivity).
#' Results are visualized through a radar plot for community assignment and bar plots for topology.
#'
#' @param control_output A list output from `community_path()` representing the ground truth network. Must contain a `graph` (igraph object) and `communities$membership`.
#' @param predicted_list A list of lists, each output from `community_path()` representing predicted networks to compare.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{community_metrics}: A data frame with VI, NMI, and ARI scores for each prediction.
#'   \item \code{topology_measures}: A data frame with raw topological metrics for each prediction.
#'   \item \code{control_topology}: A list of raw topological metrics for the ground truth network.
#' }
#'
#' @details
#' This function requires the \strong{igraph} and \strong{fmsb} packages.
#' Community similarity is measured using variation of information (VI), normalized mutual information (NMI), and adjusted Rand index (ARI).
#' Topological properties are compared by directly plotting raw values without normalization.
#'
#' @importFrom igraph modularity edge_density transitivity compare is_igraph
#' @importFrom fmsb radarchart
#' @importFrom graphics barplot par legend
#' @export
#'
#' @examples
#' # Example usage:
#' # Assuming you have outputs from community_path()
#' control <- community_path(example_graph_control)
#' predictions <- list(
#'   community_path(example_graph_pred1),
#'   community_path(example_graph_pred2)
#' )
#' result <- community_similarity(control, predictions)
#' print(result$community_metrics)
#' print(result$topology_measures)
community_similarity <- function(control_output, predicted_list) {
  required_pkgs <- c("igraph", "fmsb")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
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
    
    community_metrics[[paste0("Predicted_", i)]] <- .compare_communities(control_comm, pred_comm)
    
    if (is.null(pred_graph) || !igraph::is_igraph(pred_graph)) {
      warning("Prediction ", i, " has no valid graph. Skipping topology comparison.")
      topology_comparison[[paste0("Predicted_", i)]] <- rep(NA, 4)
      next
    }
    
    topology_comparison[[paste0("Predicted_", i)]] <- .compute_topo_metrics(pred_graph, pred_comm)
  }
  
  comm_df <- as.data.frame(do.call(rbind, community_metrics))
  topo_df <- as.data.frame(do.call(rbind, topology_comparison))
  colnames(topo_df) <- c("Modularity", "Communities", "Density", "Transitivity")
  
  .plot_radar_communities(comm_df)
  .plot_topo_barplots(topo_df, control_topo)
  
  return(list(
    community_metrics = comm_df,
    topology_measures = topo_df,
    control_topology = control_topo
  ))
}

