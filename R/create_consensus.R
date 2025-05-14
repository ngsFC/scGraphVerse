#' Create a Consensus Adjacency Matrix from Multiple Networks
#'
#' Builds a consensus adjacency matrix from a list of networks using one
#' of three methods: \code{"vote"}, \code{"union"}, or \code{"INet"}.
#'
#' @param adj_matrix_list A list of binary adjacency matrices (square,
#'   0/1) with identical dimensions and matching row/column names.
#' @param method Character string specifying the consensus strategy. One of:
#'   \itemize{
#'     \item \code{"vote"} (default): An edge is included if supported
#'       by at least \code{threshold} fraction of matrices.
#'     \item \code{"union"}: An edge is included if present in any
#'       matrix.
#'     \item \code{"INet"}: Combines normalized weighted matrices using
#'       \code{\link[INetTool]{consensusNet}}.
#'   }
#' @param weighted_list A list of weighted adjacency matrices (required if
#'   \code{method = "INet"}).
#' @param theta Numeric. Tuning parameter passed to \code{consensusNet}
#'   (default: \code{0.04}).
#' @param threshold Numeric between 0 and 1. Threshold for "vote" and
#'   "INet" methods. Default is \code{0.5}.
#' @param ncores Integer. Number of CPU cores to use when \code{method =
#'   "INet"}. Default is \code{1}.
#'
#' @return A square consensus adjacency matrix (binary or weighted,
#'   depending on the method).
#'
#' @details Consensus construction depends on the selected method:
#'   \describe{
#'     \item{\strong{vote}}{Counts the presence of each edge across all
#'       matrices and includes edges supported by at least
#'       \code{threshold × N} matrices.}
#'     \item{\strong{union}}{Includes any edge that appears in any
#'       matrix.}
#'     \item{\strong{INet}}{Multiplies binary matrices by corresponding
#'       weighted matrices, normalizes the results, and applies
#'       \code{consensusNet} to generate a consensus network.}
#'   }
#'
#' For "INet", both binary and weighted adjacency matrices must be
#' provided with matching dimensions.
#'
#' @importFrom igraph as_adjacency_matrix
#' @importFrom INetTool consensusNet
#' @export
#'
#' @examples
#' data(count_matrices)
#'
#' networks <- infer_networks(
#'     count_matrices_list = count_matrices,
#'     method = "GENIE3",
#'     nCores = 15
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
#'     nCores = 15,
#'     debug = TRUE
#' )
#' head(binary_listj[[1]])
#'
#' consensus <- create_consensus(binary_listj, method = "vote")
#' head(consensus)
create_consensus <- function(
    adj_matrix_list,
    method = "vote",
    weighted_list = NULL,
    theta = 0.04,
    threshold = 0.5,
    ncores = 1) {
    consensus_matrix <- Reduce("+", adj_matrix_list)

    if (method == "vote") {
        threshold_value <- round(length(adj_matrix_list) * threshold)
        consensus_matrix[
            consensus_matrix < threshold_value
        ] <- 0
        consensus_matrix[
            consensus_matrix >= threshold_value
        ] <- 1
    } else if (method == "union") {
        consensus_matrix[
            consensus_matrix >= 1
        ] <- 1
    } else if (method == "INet") {
        if (!requireNamespace("INetTool", quietly = TRUE)) {
            stop(
                "Package 'INetTool' is required for method = 'INet'.\n",
                "Please install it via:",
                "    remotes::install_github('ValeriaPolicastro/INet-Tool')",
                call. = FALSE
            )
        }
        if (is.null(weighted_list)) {
            stop(
                "For method 'INet', a weighted_list must be provided."
            )
        }

        list_weighted <- mapply(
            function(bin_mat, wt_mat) {
                if (!is.matrix(bin_mat) || !is.matrix(wt_mat)) {
                    stop(
                        "Both binary and weighted elements must be matrices."
                    )
                }
                if (!all(dim(bin_mat) == dim(wt_mat))) {
                    stop(
                        "Dimension mismatch between binary and weighted ",
                        "adjacency matrices."
                    )
                }
                bin_mat * wt_mat
            },
            adj_matrix_list,
            weighted_list,
            SIMPLIFY = FALSE
        )

        list_norm <- lapply(
            list_weighted,
            function(mat) mat / max(mat, na.rm = TRUE)
        )

        Con <- consensusNet(
            list_norm,
            theta     = theta,
            ncores    = ncores,
            threshold = threshold
        )
        consensus_matrix <- as.matrix(
            as_adjacency_matrix(Con$graphConsensus)
        )
    } else {
        stop(
            "Invalid method. Choose 'vote', 'union', or 'INet'."
        )
    }

    consensus_matrix
}
