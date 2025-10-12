#' Compute Community Assignment Similarity Metrics
#'
#' Calculates community assignment similarity between a reference
#' community structure and one or more predicted structures using
#' variation of information (VI), normalized mutual information (NMI),
#' and adjusted Rand index (ARI).
#'
#' @param control_output A list output from `community_path()` representing the
#'   ground truth network. Must contain `communities$membership`.
#' @param predicted_list A list of lists, each output from `community_path()`
#'   representing predicted networks to compare.
#'
#' @return A data frame with columns VI, NMI, and ARI for each prediction.
#'   Row names indicate which prediction (e.g., "Predicted_1").
#'
#' @details This function requires the \strong{igraph} package. Lower VI
#'   values indicate better similarity (VI = 0 is perfect match). Higher
#'   NMI and ARI values indicate better similarity (both range 0-1).
#'
#' @importFrom igraph compare
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
#' metrics <- compute_community_metrics(comm_truth, list(comm_cons))
compute_community_metrics <- function(control_output, predicted_list) {
    control_comm <- control_output$communities$membership

    community_metrics <- list()

    for (i in seq_along(predicted_list)) {
        pred <- predicted_list[[i]]
        pred_comm <- pred$communities$membership

        community_metrics[[paste0("Predicted_", i)]] <-
            .compare_communities(control_comm, pred_comm)
    }

    comm_df <- as.data.frame(do.call(rbind, community_metrics))
    return(comm_df)
}
