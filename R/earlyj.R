#' Modify Cell Names and Combine Datasets
#'
#' This function modifies the cell identifiers of each element in a list by appending 
#' the matrix index (e.g., "-m1", "-m2", ...) to the cell names. For matrices, if 
#' `rowg` is TRUE, it assumes that genes are in the rows (and cells are in the columns) 
#' and attaches the suffix to the column names. If `rowg` is FALSE, it first transposes the 
#' matrix so that genes become rows and cells become columns, then modifies the column names.
#'
#' For Seurat or SingleCellExperiment objects, the function standardizes features and renames 
#' the cells similarly before merging them.
#'
#' @param input_list A list of matrices, Seurat objects, or SingleCellExperiment objects.
#' @param rowg Logical. If TRUE (default), the function assumes that genes are in the rows 
#' and cells are in the columns. If FALSE, it assumes that genes are in the columns and transposes 
#' the matrices so that genes are on the rows.
#' @return A combined matrix, Seurat object, or SingleCellExperiment object with modified cell names.
#' @examples
#' \dontrun{
#' # When matrices have genes on the rows (default):
#' combined_matrix <- earlyj(list(matrix1, matrix2, matrix3))
#'
#' # When matrices have genes on the columns:
#' combined_matrix <- earlyj(list(matrix1, matrix2, matrix3), rowg = FALSE)
#' }
#' @importFrom Seurat RenameCells merge Cells
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods is
#' @export
earlyj <- function(input_list, rowg = TRUE) {
  # Ensure input is a non-empty list
  if (!is.list(input_list) || length(input_list) == 0) {
    stop("Input must be a non-empty list of matrices, Seurat objects, or SingleCellExperiment objects.")
  }
  
  # Detect object types and ensure all elements are of the same class
  object_classes <- unique(sapply(input_list, function(x) class(x)[1]))
  first_element_type <- object_classes[1]
  
  if (length(object_classes) > 1) {
    stop("All elements in input_list must be of the same type.")
  }
  
  valid_types <- c("matrix", "Seurat", "SingleCellExperiment")
  if (!(first_element_type %in% valid_types)) {
    stop("All elements must be matrices, Seurat objects, or SingleCellExperiment objects.")
  }
  
  ### Handling Matrices ###
  if (first_element_type == "matrix") {
    standardized_matrices <- lapply(seq_along(input_list), function(i) {
      mat <- input_list[[i]]
      
      # Validate matrix structure
      if (!is.matrix(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
        stop("Each matrix must be a non-empty numeric matrix.")
      }
      
      # If rowg is FALSE, assume genes are in columns and transpose so that genes become rows
      if (!rowg) {
        mat <- t(mat)
      }
      
      # Now, genes are in rows and cells are in columns.
      # Check if colnames (cell identifiers) exist; if not, assign default names.
      if (is.null(colnames(mat))) {
        colnames(mat) <- paste0("cell", seq_len(ncol(mat)))
      }
      
      # Append a unique suffix to each cell identifier
      colnames(mat) <- paste0(colnames(mat), "-m", i)
      return(mat)
    })
    
    # Combine matrices column-wise (cells as columns, genes as rows)
    combined_matrix <- do.call(cbind, standardized_matrices)
    return(combined_matrix)
  }
  
  ### Handling Seurat Objects ###
  if (first_element_type == "Seurat") {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat package is required but not installed.")
    }
    
    # Determine common features across all Seurat objects
    all_features <- Reduce(intersect, lapply(input_list, rownames))
    if (length(all_features) == 0) {
      stop("Seurat objects do not share common features. Cannot merge.")
    }
    
    modified_seurat_list <- lapply(seq_along(input_list), function(i) {
      obj <- input_list[[i]]
      obj <- subset(obj, features = all_features)  # Standardize features
      obj <- Seurat::RenameCells(obj, new.names = paste0(Seurat::Cells(obj), "-m", i))  # Rename cells
      return(obj)
    })
    
    # Merge Seurat objects
    combined_seurat <- do.call(Seurat::merge, modified_seurat_list)
    return(combined_seurat)
  }
  
  ### Handling SingleCellExperiment Objects ###
  if (first_element_type == "SingleCellExperiment") {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("SingleCellExperiment package is required but not installed.")
    }
    
    # Determine common features across all SingleCellExperiment objects
    all_features <- Reduce(intersect, lapply(input_list, rownames))
    if (length(all_features) == 0) {
      stop("SingleCellExperiment objects do not share common features. Cannot merge.")
    }
    
    modified_sce_list <- lapply(seq_along(input_list), function(i) {
      obj <- input_list[[i]]
      obj <- obj[all_features, ]  # Standardize features
      colnames(obj) <- paste0(colnames(obj), "-m", i)  # Rename cells
      return(obj)
    })
    
    # Combine SingleCellExperiment objects column-wise
    combined_sce <- do.call(cbind, modified_sce_list)
    return(combined_sce)
  }
}
