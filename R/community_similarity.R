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
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
  }

  control_comm <- control_output$communities$membership
  control_graph <- control_output$graph

  # Compute control graph topological metrics
  control_topo <- list(
    Modularity = igraph::modularity(control_graph, membership = control_comm),
    Communities = length(unique(control_comm)),
    Density = igraph::edge_density(control_graph),
    Transitivity = igraph::transitivity(control_graph)
  )

  community_metrics <- list()
  topology_comparison <- list()

  for (i in seq_along(predicted_list)) {
    pred <- predicted_list[[i]]
    pred_comm <- pred$communities$membership
    pred_graph <- pred$graph

    # --- Community similarity metrics ---
    vi <- igraph::compare(control_comm, pred_comm, method = "vi")
    nmi <- igraph::compare(control_comm, pred_comm, method = "nmi")
    ari <- igraph::compare(control_comm, pred_comm, method = "adjusted.rand")

    community_metrics[[paste0("Predicted_", i)]] <- c(VI = vi, NMI = nmi, ARI = ari)

    if (is.null(pred_graph) || !igraph::is_igraph(pred_graph)) {
      warning("Prediction ", i, " has no valid graph. Skipping topology comparison.")
      topology_comparison[[paste0("Predicted_", i)]] <- rep(NA, 4)
      next
    }

    pred_topo <- c(
      Modularity = igraph::modularity(pred_graph, pred_comm),
      Communities = length(unique(pred_comm)),
      Density = igraph::edge_density(pred_graph),
      Transitivity = igraph::transitivity(pred_graph)
    )

    topology_comparison[[paste0("Predicted_", i)]] <- pred_topo
  }

  comm_df <- as.data.frame(do.call(rbind, community_metrics))
  topo_df <- as.data.frame(do.call(rbind, topology_comparison))
  colnames(topo_df) <- c("Modularity", "Communities", "Density", "Transitivity")

  # Radar plot for Community Similarity Metrics (one plot)
  max_val <- ceiling(max(comm_df, na.rm = TRUE))
  axis_steps <- pretty(c(0, max_val), n = 5)
  radar_comm <- rbind(rep(max_val, ncol(comm_df)), rep(0, ncol(comm_df)), comm_df)

  colors <- grDevices::rainbow(nrow(comm_df))
  graphics::par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))

  fmsb::radarchart(radar_comm,
    axistype = 2, pcol = colors, plwd = 2, plty = 1,
    cglcol = "grey", axislabcol = "black",
    caxislabels = axis_steps, vlcex = 1.1,
    title = "Community Similarity Metrics"
  )

  graphics::legend("topright", legend = rownames(comm_df), col = colors, lty = 1, lwd = 2)

  # Bar plots for Topological Measures (Raw Values)
  for (i in seq_len(nrow(topo_df))) {
    pred_name <- rownames(topo_df)[i]
    pred_topo <- topo_df[i, ]

    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 4, 2))
    for (metric in colnames(pred_topo)) {
      control_value <- control_topo[[metric]]
      pred_value <- pred_topo[[metric]]
      graphics::barplot(
        height = c(control_value, pred_value),
        names.arg = c("Control", "Predicted"),
        main = paste0(metric, " Comparison\n", pred_name),
        ylab = metric,
        col = c("lightblue", "salmon")
      )
    }
  }

  return(list(
    community_metrics = comm_df,
    topology_measures = topo_df,
    control_topology = control_topo
  ))
}
