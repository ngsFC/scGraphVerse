#' Compare Community Assignments and Topology with Radar Plots
#'
#' Evaluates predicted communities against a ground truth using similarity metrics
#' and scaled topological similarity scores, visualized with dual radar plots.
#'
#' @param control_output Output from `community_path()` for the ground truth.
#' @param predicted_list List of predicted `community_path()` outputs to compare.
#'
#' @return A list with `community_metrics`, `topology_similarity` (scaled), and `control_topology`.
#' @export
community_similarity <- function(control_output, predicted_list) {
  required_pkgs <- c("igraph", "fmsb")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
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
  topology_similarity <- list()
  
  for (i in seq_along(predicted_list)) {
    pred <- predicted_list[[i]]
    pred_comm <- pred$communities$membership
    pred_graph <- pred$graph
    
    # --- Community similarity metrics (no Split-Join) ---
    vi  <- igraph::compare(control_comm, pred_comm, method = "vi")
    nmi <- igraph::compare(control_comm, pred_comm, method = "nmi")
    ari <- igraph::compare(control_comm, pred_comm, method = "adjusted.rand")
    
    community_metrics[[paste0("Predicted_", i)]] <- c(VI = vi, NMI = nmi, ARI = ari)
    
    # --- Scaled topological similarity metrics ---
    if (is.null(pred_graph) || !igraph::is_igraph(pred_graph)) {
      warning("Prediction ", i, " has no valid graph. Skipping topology comparison.")
      topology_similarity[[paste0("Predicted_", i)]] <- rep(NA, 4)
      next
    }
    
    mod_score <- 1 - abs(igraph::modularity(pred_graph, pred_comm) - control_topo$Modularity) / abs(control_topo$Modularity)
    dens_score <- 1 - abs(igraph::edge_density(pred_graph) - control_topo$Density) / abs(control_topo$Density)
    trans_score <- 1 - abs(igraph::transitivity(pred_graph) - control_topo$Transitivity) / abs(control_topo$Transitivity)
    
    pred_n_comms <- length(unique(pred_comm))
    comms_score <- 1 - abs(pred_n_comms - control_topo$Communities) / control_topo$Communities
    
    scores <- pmax(0, pmin(1, c(
      Modularity = mod_score,
      Communities = comms_score,
      Density = dens_score,
      Transitivity = trans_score
    )))
    
    topology_similarity[[paste0("Predicted_", i)]] <- scores
  }
  
  comm_df <- as.data.frame(do.call(rbind, community_metrics))
  topo_df <- as.data.frame(do.call(rbind, topology_similarity))
  colnames(topo_df) <- c("Modularity", "Communities", "Density", "Transitivity")
  
  # === Radar chart for community similarity (no Split-Join) ===
  max_val <- ceiling(max(comm_df, na.rm = TRUE))
  axis_steps <- pretty(c(0, max_val), n = 5)  # Clean numeric scale
  radar_comm <- rbind(rep(max_val, ncol(comm_df)), rep(0, ncol(comm_df)), comm_df)
  
  # === Radar chart for topological similarity (scaled 0–1) ===
  radar_topo <- rbind(rep(1, ncol(topo_df)), rep(0, ncol(topo_df)), topo_df)
  
  colors <- grDevices::rainbow(nrow(comm_df))
  graphics::par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
  
  # Community Similarity Plot (VI, NMI, ARI)
  fmsb::radarchart(radar_comm, axistype = 2, pcol = colors, plwd = 2, plty = 1,
                   cglcol = "grey", axislabcol = "black",
                   caxislabels = axis_steps,
                   vlcex = 1.1, title = "Community Similarity Metrics")
  
  # Topological Similarity Plot
  fmsb::radarchart(radar_topo, axistype = 2, pcol = colors, plwd = 2, plty = 1,
                   cglcol = "grey", axislabcol = "black", caxislabels = seq(0, 1, 0.2),
                   vlcex = 1.1, title = "Topological Similarity Scores")
  
  graphics::legend("topright", legend = rownames(comm_df), col = colors, lty = 1, lwd = 2)
  
  return(list(
    community_metrics = comm_df,
    topology_similarity = topo_df,
    control_topology = control_topo
  ))
}

