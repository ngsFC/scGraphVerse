#' Symmetrize a List of Square Matrices
#'
#' This function symmetrizes each square matrix in a list by ensuring that the values at 
#' positions \eqn{(i, j)} and \eqn{(j, i)} are identical. The symmetry is enforced by combining 
#' the two values using a specified function, such as the mean or maximum.
#'
#' @param matrix_list A list of square numeric matrices to be symmetrized.
#' @param weight_function Character or function. A function (or name of a function) 
#'   used to combine \eqn{(i, j)} and \eqn{(j, i)} values. Common options include 
#'   \code{"mean"}, \code{"max"}, \code{"min"}, or a custom function.
#' @param nCores Integer. Number of CPU cores to use for parallel processing. 
#'   Defaults to the number of available workers from \pkg{BiocParallel}.
#'
#' @return A list of symmetric matrices, where each output matrix satisfies 
#'   \eqn{A[i,j] = A[j,i]} for all \eqn{i ≠ j}.
#'
#' @details
#' For each off-diagonal pair of entries \eqn{(i, j)} and \eqn{(j, i)} in the matrix, 
#' the function:
#' \itemize{
#'   \item Uses the non-zero value if one of the two is zero.
#'   \item Applies the specified \code{weight_function} if both values are non-zero.
#' }
#' The diagonal values are not modified.
#'
#' Parallel execution is handled via \pkg{BiocParallel} for improved scalability.
#'
#' @examples
#' mat1 <- matrix(c(0, 2, 3, 4), nrow = 2)
#' mat2 <- matrix(c(0, 5, 6, 0), nrow = 2)
#' matrix_list <- list(mat1, mat2)
#' 
#' sym_list <- symmetrize(matrix_list, weight_function = "mean")
#'
#' @importFrom BiocParallel bplapply bpworkers bpparam MulticoreParam
#' @export

symmetrize <- function(matrix_list, weight_function = "mean", nCores = BiocParallel::bpworkers(BiocParallel::bpparam())) {
  if (!is.list(matrix_list) || !all(sapply(matrix_list, is.matrix))) {
    stop("matrix_list must be a list of matrices")
  }
  
  weight_function <- match.fun(weight_function)  # Ensure valid function input
  
  # Parallel processing of matrix symmetrization
  symmetrized_matrices <- BiocParallel::bplapply(matrix_list, function(mat) {
    p <- nrow(mat)
    sym_mat <- mat  # Copy structure
    
    for (i in seq_len(p - 1)) {
      for (j in seq(i + 1, p)) {
        val_ij <- mat[i, j]
        val_ji <- mat[j, i]
        
        if (val_ij == 0 || val_ji == 0) {
          symmetric_value <- max(val_ij, val_ji)  # Use the non-zero value
        } else {
          symmetric_value <- weight_function(c(val_ij, val_ji))  # Apply function
        }
        
        sym_mat[i, j] <- symmetric_value
        sym_mat[j, i] <- symmetric_value
      }
    }
    return(sym_mat)
  }, BPPARAM = BiocParallel::MulticoreParam(nCores))
  
  return(symmetrized_matrices)
}

