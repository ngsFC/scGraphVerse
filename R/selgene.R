#' Select Top Expressed Genes from Single-Cell Data
#'
#' Identifies and returns the top \code{n} most highly expressed genes across
#' all cells or within a specific cell type. Supports objects of class
#' \linkS4class{Seurat}, \linkS4class{SingleCellExperiment}, or a numeric
#' expression matrix (genes × cells).
#'
#' The function assumes that log-normalized values are available in the
#' \code{"data"} slot (for Seurat objects) or the \code{"logcounts"} assay
#' (for SingleCellExperiment). If raw counts are provided as a matrix, no
#' transformation is applied.
#'
#' Optional filtering is available to exclude mitochondrial genes
#' (\code{"^MT-"}) and ribosomal genes (\code{"^RP[SL]"}), which may otherwise
#' dominate the top expressed genes.
#'
#' @param object A Seurat object, SingleCellExperiment object, or numeric
#'   matrix (genes × cells).
#' @param top_n Integer. Number of top expressed genes to return.
#' @param cell_type Optional string. If provided, filters the expression
#'   matrix to only include cells of this type.
#' @param cell_type_col Character. Name of the column in metadata (Seurat
#'   \code{meta.data} or SCE \code{colData}) containing cell type annotations.
#'   Default is \code{"cell_type"}.
#' @param assay Character. For SingleCellExperiment objects only. Name of the
#'   assay to use. If \code{NULL}, defaults to \code{"logcounts"}.
#' @param remove_mt Logical. If \code{TRUE}, remove mitochondrial genes matching
#'   \code{"^MT-"} (case-insensitive).
#' @param remove_rib Logical. If \code{TRUE}, remove ribosomal genes matching
#'   \code{"^RP[SL]"} (case-insensitive).
#'
#' @return A character vector of the top \code{n} most highly expressed gene
#'   names.
#'
#' @section Details:
#' When using a Seurat object, the function retrieves the log-normalized data
#' from the default assay's \code{"data"} slot. For SingleCellExperiment, it
#' uses the specified assay (default is \code{"logcounts"}). For matrices, no
#' checks or transformations are applied, and subsetting by cell type is not
#' supported.
#'
#' Mitochondrial and ribosomal gene removal is based on regular expressions
#' matching gene names. These should follow standard naming conventions (e.g.,
#' \code{MT-ND1}, \code{RPL13A}, \code{RPS6}).
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' expr_mat <- matrix(rnorm(1000, mean = 3), nrow = 100, ncol = 10)
#' rownames(expr_mat) <- c(
#'     paste0("MT-", 1:5),
#'     paste0("RPL", 6:15),
#'     paste0("Gene", 16:100)
#' )
#' colnames(expr_mat) <- paste0("Cell", 1:10)
#'
#' cell_md <- DataFrame(celltype = rep(c("A", "B"), each = 5))
#'
#' sce <- SingleCellExperiment(
#'     assays   = list(logcounts = expr_mat),
#'     colData  = cell_md
#' )
#'
#' top_genes <- selgene(
#'     object         = sce,
#'     top_n          = 10,
#'     cell_type      = "A",
#'     cell_type_col  = "celltype",
#'     remove_mt      = TRUE,
#'     remove_rib     = TRUE
#' )
#'
#' @seealso \linkS4class{Seurat}, \linkS4class{SingleCellExperiment}
#' @importFrom Seurat DefaultAssay
#' @importFrom methods slotNames
#' @importFrom SummarizedExperiment assay assayNames
#' @export
selgene <- function(
    object,
    top_n,
    cell_type = NULL,
    cell_type_col = "cell_type",
    assay = NULL,
    remove_mt = FALSE,
    remove_rib = FALSE) {
    if (missing(top_n) ||
        length(top_n) != 1 ||
        !is.numeric(top_n) ||
        top_n <= 0) {
        stop(
            "Please provide a valid 'top_n' (a single positive integer)."
        )
    }

    expr <- .extract_expression(object, assay)

    if (!is.null(cell_type)) {
        expr <- .filter_by_cell_type(
            expr,
            object,
            cell_type,
            cell_type_col
        )
    }

    expr <- .filter_genes(expr, remove_mt, remove_rib)
    selected_genes <- .select_top_genes(expr, top_n)
    message(
        "Top ", top_n,
        " genes selected based on mean expression."
    )
    return(selected_genes)
}
