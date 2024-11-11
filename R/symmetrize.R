symmetrize <- function(matrix_list, weight_function = "mean") {
  
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
  
  symmetrized_matrices <- lapply(matrix_list, process_matrix, weight_function = weight_function)
  
  return(symmetrized_matrices)
}
