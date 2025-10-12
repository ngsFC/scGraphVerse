#' Symmetrize Adjacency Matrices in a SummarizedExperiment
#'
#' Symmetrizes each adjacency matrix in a \linkS4class{SummarizedExperiment}
#' by ensuring entries (i, j) and (j, i) are identical, using a specified
#' combination function.
#'
#' @param matrix_list A \linkS4class{SummarizedExperiment} object containing
#'   adjacency matrices to symmetrize.
#' @param weight_function Character string or function. Method to combine
#'   entries (i, j) and (j, i). Options include \code{"mean"},
#'   \code{"max"}, \code{"min"}, or a user-defined function.
#' @param nCores Integer. Number of CPU cores to use for parallel
#'   processing. Defaults to the number of available workers in the current
#'   \pkg{BiocParallel} backend.
#'
#' @return A \linkS4class{SummarizedExperiment} object with symmetrized
#'   adjacency matrices, where for each matrix \eqn{A[i, j] = A[j, i]} for
#'   all \eqn{i \neq j}.
#'
#' @details For each pair of off-diagonal elements (i, j) and (j, i):
#'   \itemize{
#'     \item If one value is zero, the non-zero value is used.
#'     \item If both are non-zero, they are combined using the specified
#'       \code{weight_function}.
#'   }
#'   Diagonal entries are preserved as-is and not modified.
#'
#'   Parallelization is managed via \pkg{BiocParallel} for improved
#'   performance.
#'
#' @importFrom BiocParallel bplapply bpworkers bpparam MulticoreParam
#' @export
#'
#' @examples
#' data("toy_counts")
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
symmetrize <- function(
    matrix_list,
    weight_function = "mean",
    nCores = 1) {
    if (!inherits(matrix_list, "SummarizedExperiment")) {
        stop("matrix_list must be a SummarizedExperiment object")
    }

    metadata_orig <- S4Vectors::metadata(matrix_list)
    matrix_list <- .extract_networks_from_se(matrix_list)

    weight_function <- match.fun(weight_function)

    symmetrized_matrices <- BiocParallel::bplapply(
        matrix_list,
        function(mat) {
            p <- nrow(mat)
            sym_mat <- mat

            for (i in seq_len(p - 1)) {
                for (j in seq(i + 1, p)) {
                    val_ij <- mat[i, j]
                    val_ji <- mat[j, i]

                    if (val_ij == 0 || val_ji == 0) {
                        symmetric_value <- max(val_ij, val_ji)
                    } else {
                        symmetric_value <- weight_function(c(val_ij, val_ji))
                    }

                    sym_mat[i, j] <- symmetric_value
                    sym_mat[j, i] <- symmetric_value
                }
            }

            sym_mat
        },
        BPPARAM = BiocParallel::MulticoreParam(nCores)
    )

    build_network_se(
        networks = symmetrized_matrices,
        metadata = c(metadata_orig, list(symmetrized = TRUE))
    )
}
