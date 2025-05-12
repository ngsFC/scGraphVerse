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
#' genes <- c(
#'     "TP53",   "BRCA1", "EGFR",   "MYC",     "CDKN1A",
#'     "BCL2",   "MDM2",  "PTEN",   "AKT1",    "MAPK1"
#' )
#'
#' adj <- matrix(
#'     0,
#'     nrow = length(genes),
#'     ncol = length(genes),
#'     dimnames = list(genes, genes)
#' )
#'
#' edge_list <- list(
#'     c("TP53", "MDM2"),
#'     c("TP53", "CDKN1A"),
#'     c("BRCA1", "BCL2"),
#'     c("PTEN", "AKT1"),
#'     c("EGFR", "MAPK1"),
#'     c("MYC", "CDKN1A"),
#'     c("MYC", "BCL2"),
#'     c("AKT1", "MAPK1")
#' )
#'
#' for (e in edge_list) {
#'     adj[e[1], e[2]] <- 1
#'     adj[e[2], e[1]] <- 1
#' }
#'
#' diag(adj) <- 0
#' control <- community_path(
#'     adj_matrix = adj,
#'     methods    = "louvain",
#'     pathway_db = "KEGG",
#'     genes_path = 2,
#'     plot       = FALSE,
#'     verbose    = FALSE
#' )
#'
#' pred1_adj <- adj
#' pred1_adj["TP53", "MDM2"] <- pred1_adj["MDM2", "TP53"] <- 0
#'
#' pred2_adj <- adj
#' pred2_adj["MYC", "PTEN"] <- pred2_adj["PTEN", "MYC"] <- 1
#' pred1 <- community_path(
#'     adj_matrix = pred1_adj,
#'     methods = "louvain",
#'     pathway_db = "KEGG",
#'     genes_path = 2,
#'     plot = FALSE,
#'     verbose = FALSE
#' )
#'
#' pred2 <- community_path(
#'     adj_matrix = pred2_adj,
#'     methods    = "louvain",
#'     pathway_db = "KEGG",
#'     genes_path = 2,
#'     plot       = FALSE,
#'     verbose    = FALSE
#' )
#'
#' comparison <- community_similarity(
#'     control_output = control,
#'     predicted_list = list(pred1, pred2)
#' )
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
