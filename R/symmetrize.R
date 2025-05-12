#' Symmetrize a List of Square Matrices
#'
#' Symmetrizes each square matrix in a list by ensuring entries (i, j) and
#' (j, i) are identical, using a specified combination function.
#'
#' @param matrix_list A list of square numeric matrices to symmetrize.
#' @param weight_function Character string or function. Method to combine
#'   entries (i, j) and (j, i). Options include \code{"mean"},
#'   \code{"max"}, \code{"min"}, or a user-defined function.
#' @param nCores Integer. Number of CPU cores to use for parallel
#'   processing. Defaults to the number of available workers in the current
#'   \pkg{BiocParallel} backend.
#'
#' @return A list of symmetric matrices, where for each matrix
#'   A[i, j] = A[j, i] for all i ≠ j.
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
#' mat1 <- matrix(c(0, 2, 3, 4), nrow = 2)
#' mat2 <- matrix(c(0, 5, 6, 0), nrow = 2)
#' matrix_list <- list(mat1, mat2)
#'
#' sym_list <- symmetrize(
#'     matrix_list,
#'     weight_function = "mean"
#' )
#'
#' sym_list[[1]]
symmetrize <- function(
    matrix_list,
    weight_function = "mean",
    nCores = 1) {
    if (!is.list(matrix_list) ||
        !all(vapply(matrix_list, is.matrix, logical(1)))) {
        stop("matrix_list must be a list of matrices")
    }

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

    symmetrized_matrices
}
