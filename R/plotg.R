#' Visualize Graphs from Adjacency Matrices
#'
#' Generates and arranges multiple graph visualizations from a list of
#' adjacency matrices. Each matrix is converted into an undirected
#' \pkg{igraph} object and visualized using a force-directed layout via
#' \pkg{ggraph}.
#'
#' @param adj_matrix_list A list of square, symmetric adjacency matrices
#'   with zeros on the diagonal (no self-loops), or a
#'   \linkS4class{SummarizedExperiment} object containing such matrices as
#'   assays. Each matrix represents an undirected graph.
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
#' @export
#'
#' @examples
#' data(toy_counts)
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
#' plotg(binary_se)
plotg <- function(
    adj_matrix_list) {
    BiocBaseUtils::checkInstalled("ggraph")
    BiocBaseUtils::checkInstalled("ggplot2")
    BiocBaseUtils::checkInstalled("gridExtra")

    # Accept SummarizedExperiment and extract assays
    if (inherits(adj_matrix_list, "SummarizedExperiment")) {
        adj_matrix_list <- as.list(
            SummarizedExperiment::assays(adj_matrix_list)
        )
    }

    if (!is.list(adj_matrix_list)) {
        stop("Input must be a list of matrices or a SummarizedExperiment.")
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
