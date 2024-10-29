symmetrize <- function(matrix_list, weight_function = "mean") {
  process_matrix <- function(mat, weight_function) {
    # Ensure the matrix is symmetric by averaging corresponding pairs
    p <- nrow(mat)
    sym_mat <- mat  # Initialize with the original matrix
    
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        # Calculate the average of mat[i, j] and mat[j, i] using the weight_function
        symmetric_value <- match.fun(weight_function)(c(mat[i, j], mat[j, i]))
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