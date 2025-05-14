#' Modify Cell Names and Combine Datasets
#'
#' Modifies cell identifiers of each element in a list of matrices,
#' \linkS4class{Seurat} objects, or \linkS4class{SingleCellExperiment} objects
#' by appending a unique matrix index (e.g., "-m1", "-m2", etc.).
#' After renaming, the datasets are merged into a single object.
#'
#' @param input_list A list of matrices, \linkS4class{Seurat} objects,
#' or \linkS4class{SingleCellExperiment} objects. All elements must be of the
#' same class.
#' @param rowg Logical. If \code{TRUE} (default), genes are assumed to be rows
#'  and cells columns.
#' If \code{FALSE}, matrices are transposed before renaming and combining.
#'
#' @return
#' A combined matrix, Seurat object, or SingleCellExperiment object with
#' modified (unique) cell names.
#'
#' @details
#' For matrices, this function optionally transposes the input before combining.
#'  For \code{Seurat} and \code{SingleCellExperiment} objects, only features
#'  (genes) common across all input datasets are retained before merging.
#'  The cell names are suffixed with "-m1", "-m2", etc., according to their
#'  original list position.
#'  
#' @importFrom SeuratObject merge.Assay subset.Assay
#' @importFrom Seurat RenameCells Cells
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods is
#' @export
#'
#' @examples
#' data(count_matrices)
#'
#' earlym <- earlyj(count_matrices)
#' head(earlym)
#'
earlyj <- function(input_list, rowg = TRUE) {
    if (!is.list(input_list) || length(input_list) == 0) {
        stop("Input must be a list of matrices, Seurat obj, or SCE obj")
    }

    object_classes <- unique(
        vapply(input_list, function(x) class(x)[1], character(1))
    )

    if (length(object_classes) > 1) {
        stop("All elements in input_list must be of the same type.")
    }

    first_type <- object_classes[1]
    valid_types <- c("matrix", "Seurat", "SingleCellExperiment")
    if (!(first_type %in% valid_types)) {
        stop("All elements must be matrices, Seurat obj, or SCE obj")
    }

    if (first_type == "matrix") {
        return(.merge_matrix_list(input_list, rowg))
    }

    if (first_type == "Seurat") {
        return(.merge_seurat_list(input_list))
    }

    if (first_type == "SingleCellExperiment") {
        return(.merge_sce_list(input_list))
    }
}
