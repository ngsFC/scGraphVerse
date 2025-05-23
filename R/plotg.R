#' Visualize Graphs from Adjacency Matrices
#'
#' Generates and arranges multiple graph visualizations from a list of
#' adjacency matrices. Each matrix is converted into an undirected
#' \pkg{igraph} object and visualized using a force-directed layout via
#' \pkg{ggraph}.
#'
#' @param adj_matrix_list A list of square, symmetric adjacency matrices
#'   with zeros on the diagonal (no self-loops). Each matrix represents an
#'   undirected graph.
#'
#' @return A grid of plots displaying all valid graphs in the input list.
#'
#' @details Each adjacency matrix is validated to ensure it is square and
#'   symmetric. Disconnected nodes (degree zero) are removed prior to
#'   visualization. Graphs are visualized with a force-directed layout
#'   using \pkg{ggraph}, and multiple plots are arranged into a grid with
#'   \pkg{gridExtra}.
#'
#'   Each subplot title includes the graph index, number of nodes, and
#'   number of edges.
#'
#' @note This function requires the following packages:
#'   \pkg{igraph}, \pkg{ggraph}, and \pkg{gridExtra}. If any are missing,
#'   an informative error will be thrown.
#'
#' @importFrom igraph graph_from_adjacency_matrix delete_vertices V
#'   degree vcount ecount
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 labs theme_minimal theme element_text
#' @export
#'
#' @examples
#' data(count_matrices)
#'
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
#' plotg(binary_listj)
plotg <- function(
    adj_matrix_list) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("The 'igraph' package is required but not installed.")
    }
    if (!requireNamespace("ggraph", quietly = TRUE)) {
        stop("The 'ggraph' package is required but not installed.")
    }
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
        stop("The 'gridExtra' package is required but not installed.")
    }

    if (!is.list(adj_matrix_list)) {
        stop("Input must be a list of adjacency matrices.")
    }

    plot_list <- list()

    for (i in seq_along(adj_matrix_list)) {
        mat <- adj_matrix_list[[i]]

        # Inline validation
        if (!is.matrix(mat) ||
            nrow(mat) != ncol(mat) ||
            !all(mat == t(mat))) {
            warning(sprintf(
                "Skipping graph %d: Adjacency matrix must symmetric.",
                i
            ))
            next
        }

        p <- .create_igraph_plot(mat, i)
        if (is.null(p)) {
            warning(sprintf(
                "Skipping graph %d: No edges remaining after filtering.",
                i
            ))
            next
        }

        plot_list[[length(plot_list) + 1]] <- p
    }

    if (length(plot_list) == 0) {
        stop("No valid graphs were generated.")
    }

    gridExtra::grid.arrange(
        grobs = plot_list,
        ncol  = ceiling(sqrt(length(plot_list)))
    )
}
