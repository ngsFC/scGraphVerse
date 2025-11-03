#' Modify Cell Names and Combine Datasets
#'
#' Extracts expression data from a \linkS4class{MultiAssayExperiment} object,
#' modifies cell identifiers by appending a unique experiment index
#' (e.g., "-m1", "-m2", etc.), and merges the datasets into a single object.
#'
#' @param input_list A \linkS4class{MultiAssayExperiment} object containing
#'   expression data from multiple experiments or conditions.
#' @param rowg Logical. If \code{TRUE} (default), genes are assumed to be rows
#'  and cells columns.
#' If \code{FALSE}, matrices are transposed before renaming and combining.
#'
#' @return
#' A \linkS4class{MultiAssayExperiment} object containing a single merged
#' experiment with modified (unique) cell names.
#'
#' @details
#' For matrices, this function optionally transposes the input before combining.
#'  For \code{Seurat} and \code{SingleCellExperiment} objects, only features
#'  (genes) common across all input datasets are retained before merging.
#'  The cell names are suffixed with "-m1", "-m2", etc., according to their
#'  original list position. The result is returned as a MultiAssayExperiment
#'  with a single merged experiment.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom methods is
#' @export
#'
#' @examples
#' data(toy_counts)
#'
#' merged_mae <- earlyj(toy_counts)
#' merged_mae
#'
earlyj <- function(input_list, rowg = TRUE) {
    if (!inherits(input_list, "MultiAssayExperiment")) {
        stop("input_list must be a MultiAssayExperiment object")
    }

    input_list <- .extract_from_mae(input_list)

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

    merged_data <- if (first_type == "matrix") {
        .merge_matrix_list(input_list, rowg)
    } else if (first_type == "Seurat") {
        .merge_seurat_list(input_list)
    } else if (first_type == "SingleCellExperiment") {
        .merge_sce_list(input_list)
    }

    # Return as MultiAssayExperiment with single merged experiment
    create_mae(list(merged = merged_data))
}
