#' Symmetrize a List of Square Matrices
#'
#' Symmetrizes each square matrix in a list by ensuring that entries \eqn{(i,j)} and \eqn{(j,i)}
#' are identical, based on a specified combination function.
#'
#' @param matrix_list A list of square numeric matrices to symmetrize.
#' @param weight_function Character string or function. Method to combine \eqn{(i,j)} and \eqn{(j,i)} entries.
#'   Options include \code{"mean"}, \code{"max"}, \code{"min"}, or a user-defined function.
#' @param nCores Integer. Number of CPU cores to use for parallel processing.
#'   Defaults to the number of available workers in the current \pkg{BiocParallel} backend.
#'
#' @return
#' A list of symmetric matrices, where for each matrix \eqn{A[i,j] = A[j,i]} for all \eqn{i ≠ j}.
#'
#' @details
#' For each pair of off-diagonal elements \eqn{(i,j)} and \eqn{(j,i)}:
#' \itemize{
#'   \item If one value is zero, the non-zero value is used.
#'   \item If both values are non-zero, they are combined using the specified \code{weight_function}.
#' }
#' Diagonal entries are preserved as-is and not modified.
#'
#' Parallelization is managed via \pkg{BiocParallel} for improved performance.
#'
#' @importFrom BiocParallel bplapply bpworkers bpparam MulticoreParam
#' @export
#'
#' @examples
#' # Create two small asymmetric matrices
#' mat1 <- matrix(c(0, 2, 3, 4), nrow = 2)
#' mat2 <- matrix(c(0, 5, 6, 0), nrow = 2)
#' matrix_list <- list(mat1, mat2)
#'
#' # Symmetrize using the mean function
#' sym_list <- symmetrize(matrix_list, weight_function = "mean")
#'
#' # View the first symmetrized matrix
#' sym_list[[1]]

symmetrize <- function(matrix_list, weight_function = "mean", nCores = BiocParallel::bpworkers(BiocParallel::bpparam())) {
  if (!is.list(matrix_list) || !all(sapply(matrix_list, is.matrix))) {
    stop("matrix_list must be a list of matrices")
  }
  
  weight_function <- match.fun(weight_function) 
  
  symmetrized_matrices <- BiocParallel::bplapply(matrix_list, function(mat) {
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
    return(sym_mat)
  }, BPPARAM = BiocParallel::MulticoreParam(nCores))
  
  return(symmetrized_matrices)
}

