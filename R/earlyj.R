earlyj  <- function(matrix_list) {
  # Initialize an empty list to store the modified matrices
  modified_matrices <- list()
  
  # Loop through each matrix in the list and modify its row names
  for (i in seq_along(matrix_list)) {
    # Get the current matrix
    current_matrix <- matrix_list[[i]]
    
    # Modify row names to include the matrix index (m1, m2, ...)
    new_rownames <- paste0(rownames(current_matrix), "_m", i)
    rownames(current_matrix) <- new_rownames
    
    # Store the modified matrix in the list
    modified_matrices[[i]] <- current_matrix
  }
  
  # Combine all modified matrices using rbind
  combined_matrix <- do.call(rbind, modified_matrices)
  
  return(combined_matrix)
}
