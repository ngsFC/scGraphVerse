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
#' \dontrun{
#' # Example list of matrices
#' mat1 <- matrix(c(0, 2, 3, 4), nrow = 2, ncol = 2)
#' mat2 <- matrix(c(0, 5, 6, 0), nrow = 2, ncol = 2)
#' matrix_list <- list(mat1, mat2)
#' 
#' # Symmetrize with "mean" function
#' symmetrized_matrices <- symmetrize(matrix_list, weight_function = "mean")
#' }
#'
#' @export
symmetrize <- function(matrix_list, weight_function = "mean") {
  
  # Helper function to process and symmetrize a single matrix
  process_matrix <- function(mat, weight_function) {
    p <- nrow(mat)
    sym_mat <- mat  
    
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        # Get values at (i, j) and (j, i)
        val_ij <- mat[i, j]
        val_ji <- mat[j, i]
        
        # Check for zeros and apply logic
        if (val_ij == 0 || val_ji == 0) {
          # Use the non-zero or the higher value if one is zero
          symmetric_value <- max(val_ij, val_ji)
        } else {
          # Use the specified weight function if both are non-zero
          symmetric_value <- match.fun(weight_function)(c(val_ij, val_ji))
        }
        
        # Set both (i, j) and (j, i) to the symmetric value
        sym_mat[i, j] <- symmetric_value
        sym_mat[j, i] <- symmetric_value
      }
    }
    
    return(sym_mat)
  }
  
  # Apply the process_matrix function to each matrix in the list
  symmetrized_matrices <- lapply(matrix_list, process_matrix, weight_function = weight_function)
  
  return(symmetrized_matrices)
}

