#' Create MultiAssayExperiment from Multiple Single-Cell Datasets
#'
#' Converts a list of count matrices, Seurat objects, or SingleCellExperiment
#' objects into a MultiAssayExperiment for integrated network inference.
#'
#' @param datasets A named list of datasets. Each element can be:
#'   \itemize{
#'     \item A matrix (genes × cells)
#'     \item A Seurat object
#'     \item A SingleCellExperiment object
#'   }
#' @param colData Optional. A DataFrame with metadata for each experiment.
#'   If NULL, automatically generated from list names.
#' @param ... Additional arguments (currently unused)
#'
#' @return A MultiAssayExperiment object with:
#'   \itemize{
#'     \item experiments: List of SingleCellExperiment objects
#'     \item colData: Metadata for each experiment/condition
#'   }
#'
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame
#' @export
#'
#' @examples
#' # Load the example MAE
#' data("toy_counts")
#'
#' # Extract the list of SingleCellExperiment objects
#' sce_list <- MultiAssayExperiment::experiments(toy_counts)
#' sce_list <- as.list(sce_list)
#'
#' # Create a new MAE from the SCE list
#' mae <- create_mae(sce_list)
#' mae
create_mae <- function(datasets, colData = NULL, ...) {
    if (!is.list(datasets)) {
        stop("datasets must be a named list")
    }

    if (is.null(names(datasets))) {
        names(datasets) <- paste0("experiment_", seq_along(datasets))
    }

    sce_list <- .convert_to_sce_list(datasets)

    # For network inference with different cells per experiment,
    # we use ExperimentList without explicit sample mapping
    MultiAssayExperiment::MultiAssayExperiment(
        experiments = sce_list
    )
}
