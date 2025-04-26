#' Modify Cell Names and Combine Datasets
#'
#' Modifies cell identifiers of each element in a list of matrices, \linkS4class{Seurat} objects,
#' or \linkS4class{SingleCellExperiment} objects by appending a unique matrix index (e.g., "-m1", "-m2", etc.).
#' After renaming, the datasets are merged into a single object.
#'
#' @param input_list A list of matrices, \linkS4class{Seurat} objects, or \linkS4class{SingleCellExperiment} objects.
#'   All elements must be of the same class.
#' @param rowg Logical. If \code{TRUE} (default), genes are assumed to be rows and cells columns.
#'   If \code{FALSE}, matrices are transposed before renaming and combining.
#'
#' @return
#' A combined matrix, Seurat object, or SingleCellExperiment object with modified (unique) cell names.
#'
#' @details
#' For matrices, this function optionally transposes the input before combining. For \code{Seurat}
#' and \code{SingleCellExperiment} objects, only features (genes) common across all input datasets are retained
#' before merging. The cell names are suffixed with "-m1", "-m2", etc., according to their original list position.
#'
#' @importFrom Seurat RenameCells merge Cells
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods is
#' @export
#'
#' @examples
#' \dontrun{
#' # Example with matrices where genes are rows (default behavior)
#' mat1 <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
#' mat2 <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
#' rownames(mat1) <- paste0("Gene", 1:10)
#' rownames(mat2) <- paste0("Gene", 1:10)
#'
#' combined_matrix <- earlyj(list(mat1, mat2))
#'
#' # Example with matrices where genes are columns
#' mat3 <- t(mat1)
#' mat4 <- t(mat2)
#'
#' combined_matrix2 <- earlyj(list(mat3, mat4), rowg = FALSE)
#'
#' # Example with Seurat objects (requires Seurat installed)
#' # combined_seurat <- earlyj(list(seurat_obj1, seurat_obj2))
#'
#' # Example with SingleCellExperiment objects (requires SingleCellExperiment installed)
#' # combined_sce <- earlyj(list(sce_obj1, sce_obj2))
#' }

earlyj <- function(input_list, rowg = TRUE) {
  if (!is.list(input_list) || length(input_list) == 0) {
    stop("Input must be a non-empty list of matrices, Seurat objects, or SingleCellExperiment objects.")
  }
  
  object_classes <- unique(sapply(input_list, function(x) class(x)[1]))
  first_element_type <- object_classes[1]
  
  if (length(object_classes) > 1) {
    stop("All elements in input_list must be of the same type.")
  }
  
  valid_types <- c("matrix", "Seurat", "SingleCellExperiment")
  if (!(first_element_type %in% valid_types)) {
    stop("All elements must be matrices, Seurat objects, or SingleCellExperiment objects.")
  }
  
  if (first_element_type == "matrix") {
    standardized_matrices <- lapply(seq_along(input_list), function(i) {
      mat <- input_list[[i]]
      
      if (!is.matrix(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
        stop("Each matrix must be a non-empty numeric matrix.")
      }
      
      if (!rowg) {
        mat <- t(mat)
      }
      
      if (is.null(colnames(mat))) {
        colnames(mat) <- paste0("cell", seq_len(ncol(mat)))
      }
      
      colnames(mat) <- paste0(colnames(mat), "-m", i)
      return(mat)
    })
    
    combined_matrix <- do.call(cbind, standardized_matrices)
    return(combined_matrix)
  }
  
  if (first_element_type == "Seurat") {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat package is required but not installed.")
    }
    
    all_features <- Reduce(intersect, lapply(input_list, rownames))
    if (length(all_features) == 0) {
      stop("Seurat objects do not share common features. Cannot merge.")
    }
    
    modified_seurat_list <- lapply(seq_along(input_list), function(i) {
      obj <- input_list[[i]]
      obj <- subset(obj, features = all_features)
      obj <- Seurat::RenameCells(obj, new.names = paste0(Seurat::Cells(obj), "-m", i))
      return(obj)
    })
    
    combined_seurat <- do.call(Seurat::merge, modified_seurat_list)
    return(combined_seurat)
  }
  
  if (first_element_type == "SingleCellExperiment") {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("SingleCellExperiment package is required but not installed.")
    }
    
    all_features <- Reduce(intersect, lapply(input_list, rownames))
    if (length(all_features) == 0) {
      stop("SingleCellExperiment objects do not share common features. Cannot merge.")
    }
    
    modified_sce_list <- lapply(seq_along(input_list), function(i) {
      obj <- input_list[[i]]
      obj <- obj[all_features, ]
      colnames(obj) <- paste0(colnames(obj), "-m", i)
      return(obj)
    })
    
    combined_sce <- do.call(cbind, modified_sce_list)
    return(combined_sce)
  }
}
