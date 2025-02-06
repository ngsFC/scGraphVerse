#' Modify Row Names of Matrices and Combine Them
#'
#' This function modifies the row names of each matrix in a list by appending 
#' the matrix index (e.g., "_m1", "_m2", ...) to the row names. It then combines 
#' all the modified matrices into a single matrix using `rbind`.
#'
#' @param matrix_list A list of matrices. Each matrix will have its row names modified.
#' @return A combined matrix with modified row names.
#' @examples
#' \dontrun{
#' combined_matrix <- earlyj(list(matrix1, matrix2, matrix3))
#' }
#' @export
earlyj <- function(input_list) {
  # Ensure input is a list
  if (!is.list(input_list) || length(input_list) == 0) {
    stop("Input must be a non-empty list of matrices, Seurat objects, or SingleCellExperiment objects.")
  }
  
  # Detect the type of input
  first_element_type <- class(input_list[[1]])
  valid_types <- c("matrix", "Seurat", "SingleCellExperiment")
  
  if (!first_element_type %in% valid_types) {
    stop("All elements must be matrices, Seurat objects, or SingleCellExperiment objects.")
  }
  
  # If input is a list of matrices, standardize them before combining
  if (first_element_type == "matrix") {
    standardized_matrices <- lapply(seq_along(input_list), function(i) {
      mat <- input_list[[i]]
      
      # Check if cells are in columns (genes in rows), and transpose if needed
      if (ncol(mat) > nrow(mat)) {
        mat <- t(mat)
      }
      
      # Modify cell names (row names)
      if (is.null(rownames(mat))) {
        rownames(mat) <- paste0("cell", seq_len(nrow(mat)))
      }
      rownames(mat) <- paste0(rownames(mat), "-m", i)
      
      return(mat)
    })
    
    # Combine matrices using rbind (since cells are in rows)
    combined_matrix <- do.call(rbind, standardized_matrices)
    return(combined_matrix)
  }
  
  # If input is a list of Seurat objects, combine them
  if (first_element_type == "Seurat") {
    # Ensure all Seurat objects use the same feature (gene) set
    all_features <- Reduce(intersect, lapply(input_list, rownames))
    
    if (length(all_features) == 0) {
      stop("Seurat objects do not share common features. Cannot merge.")
    }
    
    modified_seurat_list <- lapply(seq_along(input_list), function(i) {
      obj <- input_list[[i]]
      
      # Subset to common features to avoid errors
      obj <- subset(obj, features = all_features)
      
      # Modify cell names
      obj <- RenameCells(obj, new.names = paste0(Cells(obj), "-m", i))
      
      return(obj)
    })
    
    # Combine Seurat objects using merge
    combined_seurat <- merge(modified_seurat_list[[1]], y = modified_seurat_list[-1])
    return(combined_seurat)
  }
  
  # If input is a list of SingleCellExperiment objects, combine them
  if (first_element_type == "SingleCellExperiment") {
    # Ensure all SCE objects have the same features (genes)
    all_features <- Reduce(intersect, lapply(input_list, rownames))
    
    if (length(all_features) == 0) {
      stop("SingleCellExperiment objects do not share common features. Cannot merge.")
    }
    
    modified_sce_list <- lapply(seq_along(input_list), function(i) {
      obj <- input_list[[i]]
      
      # Subset to common features
      obj <- obj[all_features, ]
      
      # Modify cell names
      colnames(obj) <- paste0(colnames(obj), "-m", i)
      
      return(obj)
    })
    
    # Combine SCE objects using cbind
    combined_sce <- do.call(cbind, modified_sce_list)
    return(combined_sce)
  }
}
