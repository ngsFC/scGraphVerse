#' Load an AnnData H5AD File into Seurat or SingleCellExperiment
#'
#' @title Load H5AD into Seurat or SingleCellExperiment Object
#'
#' @description
#' Download an \*.h5ad file from a URL to a temporary file, read it via Python’s
#' \pkg{anndata} (through \pkg{reticulate}), and return either a
#' \code{\link[Seurat]{Seurat-class}} or a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment-class}} object.
#'
#' @details
#' This function hides all Python-side dependencies (auto-installs
#' \pkg{anndata} if missing), downloads the specified URL via
#' \code{\link[utils]{download.file}} into a temp file, reads it with
#' \code{anndata::read_h5ad()}, and constructs your chosen R object.
#'
#' @param url          Character(1). URL pointing to the \*.h5ad file.
#' @param output       Character(1). Either \code{"seurat"} (default) for a
#'                     \pkg{Seurat} object, or \code{"sce"} for a
#'                     \pkg{SingleCellExperiment} object.
#' @param min.features Integer(1). Minimum number of features detected per cell
#'                     (applies only when \code{output = "seurat"}). Default: 500.
#' @param min.cells    Integer(1). Minimum number of cells expressing each feature
#'                     (applies only when \code{output = "seurat"}). Default: 30.
#'
#' @return
#' A \code{\link[Seurat]{Seurat-class}} object if \code{output = "seurat"}, or
#' a \code{\link[SingleCellExperiment]{SingleCellExperiment-class}} object if
#' \code{output = "sce"}.
#'
#' @author
#' Your Name
#'
#' @examples
#' # Seurat object
#' pbmc_seurat <- load_h5ad(
#'   url    = "https://datasets.cellxgene.cziscience.com/42e24d96-69c6-4f05-8237-b5d305c49b45.h5ad",
#'   output = "seurat",
#'   min.features = 500,
#'   min.cells    = 30
#' )
#'
#' # SingleCellExperiment object
#' pbmc_sce <- load_h5ad(
#'   url    = "https://datasets.cellxgene.cziscience.com/42e24d96-69c6-4f05-8237-b5d305c49b45.h5ad",
#'   output = "sce"
#' )
#'
#' @import reticulate
#' @importFrom utils download.file
#' @importFrom Seurat CreateSeuratObject
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom Matrix Matrix
#' @export
download_Atlas <- function(
    url,
    output       = c("seurat", "sce"),
    min.features = 500,
    min.cells    = 30
) {
  output <- match.arg(output)
  
  # Ensure Python anndata is installed
  if (!reticulate::py_module_available("anndata")) {
    reticulate::py_install("anndata", pip = TRUE)
  }
  adata_mod <- reticulate::import("anndata")
  
  # Download to a temp file
  tmpfile <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmpfile), add = TRUE)
  utils::download.file(url, destfile = tmpfile, mode = "wb", quiet = TRUE)
  
  adata <- adata_mod$read_h5ad(tmpfile)
  
  counts_mat <- t(as.matrix(adata$X))
  meta_df    <- reticulate::py_to_r(adata$obs)
  
  if (output == "seurat") {
    obj <- Seurat::CreateSeuratObject(
      counts       = counts_mat,
      meta.data    = meta_df,
      min.features = min.features,
      min.cells    = min.cells
    )
  } else {
    sparse_mat <- Matrix::Matrix(counts_mat, sparse = TRUE)
    coldata    <- S4Vectors::DataFrame(meta_df)
    obj <- SingleCellExperiment::SingleCellExperiment(
      assays  = list(counts = sparse_mat),
      colData = coldata
    )
  }
  
  obj
}
