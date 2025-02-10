#' Symmetrize a List of Matrices
#'
#' This function takes a list of matrices and symmetrizes each matrix by ensuring that 
#' for each pair of off-diagonal elements (i, j) and (j, i), their values are the same. 
#' The function allows you to specify how to combine the values at (i, j) and (j, i) 
#' using a weight function, such as the "mean" or "max" of the two values.
#'
#' @param matrix_list A list of matrices to be symmetrized. Each matrix should be square 
#'   and contain numeric values.
#' @param weight_function A string specifying the function to combine values at (i, j) and (j, i). 
#'   The default is "mean". Other options include "max", or any other function that can be 
#'   applied to a vector of two numeric values (e.g., "min").
#' @param nCores Number of cores for parallel processing. Default is `BiocParallel::bpworkers(BiocParallel::bpparam())`.
#' 
#' @return A list of symmetrized matrices. Each matrix will have identical values at positions 
#'   (i, j) and (j, i), computed using the specified weight function.
#'
#' @details 
#' For each matrix in the input list, the function processes every pair of off-diagonal elements 
#' and applies the chosen `weight_function` to symmetrize the matrix. If one of the values at 
#' (i, j) or (j, i) is zero, the function will use the non-zero value. If both are non-zero, 
#' the function will apply the `weight_function` to the two values.
#'
#' @examples
#' 
#' # Example list of matrices
#' mat1 <- matrix(c(0, 2, 3, 4), nrow = 2, ncol = 2)
#' mat2 <- matrix(c(0, 5, 6, 0), nrow = 2, ncol = 2)
#' matrix_list <- list(mat1, mat2)
#' 
#' # Symmetrize with "mean" function
#' symmetrized_matrices <- symmetrize(matrix_list, weight_function = "mean")
#'
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

