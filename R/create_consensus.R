#' Create a Consensus Adjacency Matrix from Multiple Networks
#'
#' Builds a consensus adjacency matrix from networks stored in a
#' \linkS4class{SummarizedExperiment} using one of three methods:
#' \code{"vote"}, \code{"union"}, or \code{"INet"}.
#'
#' @param adj_matrix_list A \linkS4class{SummarizedExperiment} object
#'   containing binary adjacency matrices (square, 0/1) with identical
#'   dimensions and matching row/column names, or a list of such matrices.
#' @param method Character string specifying the consensus strategy. One of:
#'   \itemize{
#'     \item \code{"vote"} (default): An edge is included if supported
#'       by at least \code{threshold} fraction of matrices.
#'     \item \code{"union"}: An edge is included if present in any
#'       matrix.
#'     \item \code{"INet"}: Combines normalized weighted matrices using
#'       \code{\link[INetTool]{consensusNet}}.
#'   }
#' @param weighted_list A \linkS4class{SummarizedExperiment} object containing
#'   weighted adjacency matrices (required if \code{method = "INet"}), or a
#'   list of such matrices.
#' @param theta Numeric. Tuning parameter passed to \code{consensusNet}
#'   (default: \code{0.04}).
#' @param threshold Numeric between 0 and 1. Threshold for "vote" and
#'   "INet" methods. Default is \code{0.5}.
#' @param ncores Integer. Number of CPU cores to use when \code{method =
#'   "INet"}. Default is \code{1}.
#' @param tolerance Numeric. Tolerance for differences between similar graphs
#'   in INet method. Default is \code{0.1}.
#' @param nitermax Integer. Maximum number of iterations for INet algorithm.
#'   Default is \code{50}.
#' @param verbose Logical. If TRUE, display verbose output for INet method.
#'   Default is \code{FALSE}.
#'
#' @return A \linkS4class{SummarizedExperiment} object with a single assay
#'   containing the consensus adjacency matrix (binary or weighted, depending
#'   on the method). Metadata includes consensus method and parameters.
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
#'
#' consensus <- create_consensus(binary_se, method = "union")
#' head(consensus)
create_consensus <- function(
    adj_matrix_list,
    method = c("vote", "union", "INet"),
    weighted_list = NULL,
    theta = 0.04,
    threshold = 0.5,
    ncores = 1,
    tolerance = 0.1,
    nitermax = 50,
    verbose = FALSE) {
    method <- match.arg(method)

    if (inherits(adj_matrix_list, "SummarizedExperiment")) {
        adj_matrix_list <- .extract_networks_from_se(adj_matrix_list)
    }

    if (!is.null(weighted_list) && inherits(
        weighted_list,
        "SummarizedExperiment"
    )) {
        weighted_list <- .extract_networks_from_se(weighted_list)
    }

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

        if (!requireNamespace("INetTool", quietly = TRUE)) {
            stop(
                "Package 'INetTool' is required for method = 'INet'. ",
                "Install it with: BiocManager::install('INetTool')"
            )
        }

        Con <- INetTool::consensusNet(
            list_norm,
            theta     = theta,
            ncores    = ncores,
            threshold = threshold,
            tolerance = tolerance,
            nitermax  = nitermax,
            verbose   = verbose
        )
        consensus_matrix <- as.matrix(
            igraph::as_adjacency_matrix(Con$graphConsensus)
        )
    }

    # Return as SummarizedExperiment
    build_network_se(
        networks = list(consensus = consensus_matrix),
        networkData = S4Vectors::DataFrame(
            network = "consensus",
            n_edges = sum(consensus_matrix > 0),
            n_input_networks = length(adj_matrix_list),
            row.names = "consensus"
        ),
        metadata = list(
            type = if (method == "INet") "weighted" else "binary",
            method = method,
            threshold = threshold
        )
    )
}
