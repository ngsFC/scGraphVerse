#' Select Top Expressed Genes from Single-Cell Data
#'
#' Identifies and returns the top \code{n} most highly expressed genes across all cells
#' or within a specific cell type. Supports objects of class \linkS4class{Seurat},
#' \linkS4class{SingleCellExperiment}, or a numeric expression matrix (genes × cells).
#'
#' The function assumes that log-normalized values are available in the
#' \code{"data"} slot (for Seurat objects) or the \code{"logcounts"} assay (for
#' SingleCellExperiment). If raw counts are provided as a matrix, no transformation
#' is applied.
#'
#' Optional filtering is available to exclude mitochondrial genes (\code{"^MT-"})
#' and ribosomal genes (\code{"^RP[SL]"}), which may otherwise dominate the top expressed genes.
#'
#' @param object A Seurat object, SingleCellExperiment object, or numeric matrix (genes × cells).
#' @param top_n Integer. Number of top expressed genes to return.
#' @param cell_type Optional string. If provided, filters the expression matrix to only include cells of this type.
#' @param cell_type_col Character. Name of the column in metadata (Seurat \code{meta.data} or SCE \code{colData})
#'                      containing cell type annotations. Default is \code{"cell_type"}.
#' @param assay Character. For SingleCellExperiment objects only. Name of the assay to use.
#'              If NULL, defaults to \code{"logcounts"}.
#' @param remove_mt Logical. If TRUE, remove mitochondrial genes matching \code{"^MT-"} (case-insensitive).
#' @param remove_rib Logical. If TRUE, remove ribosomal genes matching \code{"^RP[SL]"} (case-insensitive).
#'
#' @return A character vector of the top \code{n} most highly expressed gene names.
#'
#' @section Details:
#' When using a Seurat object, the function retrieves the log-normalized data
#' from the default assay's \code{"data"} slot. For SingleCellExperiment, it uses
#' the specified assay (default is \code{"logcounts"}). For matrices, no checks or
#' transformations are applied, and subsetting by cell type is not supported.
#'
#' Mitochondrial and ribosomal gene removal is based on regular expressions matching
#' gene names. These should follow standard naming conventions (e.g., \code{MT-ND1},
#' \code{RPL13A}, \code{RPS6}).
#'
#' @examples
#' # Example with mock SingleCellExperiment
#' library(SingleCellExperiment)
#'
#' # Simulated expression matrix (log-normalized)
#' expr_mat <- matrix(rnorm(1000, mean = 3), nrow = 100, ncol = 10)
#' rownames(expr_mat) <- c(paste0("MT-", 1:5), paste0("RPL", 6:15), paste0("Gene", 16:100))
#' colnames(expr_mat) <- paste0("Cell", 1:10)
#'
#' # Simulated metadata with cell type
#' cell_md <- DataFrame(celltype = rep(c("A", "B"), each = 5))
#'
#' # Create SCE object
#' sce <- SingleCellExperiment(assays = list(logcounts = expr_mat), colData = cell_md)
#'
#' # Run function with mitochondrial and ribosomal gene removal
#' top_genes <- selgene(sce,
#'   top_n = 10, cell_type = "A", cell_type_col = "celltype",
#'   remove_mt = TRUE, remove_rib = TRUE
#' )
#' @seealso \linkS4class{Seurat}, \linkS4class{SingleCellExperiment}
#' @importFrom Seurat DefaultAssay
#' @importFrom methods slotNames
#' @importFrom SummarizedExperiment assay assayNames
#' @export

selgene <- function(object, top_n,
                    cell_type = NULL, cell_type_col = "cell_type", assay = NULL,
                    remove_mt = FALSE, remove_rib = FALSE) {
  if (missing(top_n) || length(top_n) != 1 || !is.numeric(top_n) || top_n <= 0) {
    stop("Please provide a valid 'top_n' (a single positive integer).")
  }

  if (inherits(object, "Seurat")) {
    assay_name <- Seurat::DefaultAssay(object)
    seurat_assay <- object[[assay_name]]
    slots_avail <- methods::slotNames(seurat_assay)

    if ("data" %in% slots_avail) {
      expr <- seurat_assay@data
      message("Using Seurat assay '", assay_name, "' slot 'data' (log-normalized).")
    } else {
      stop(
        "Assay '", assay_name, "' has no 'data' slot. Available slots: ",
        paste(slots_avail, collapse = ", ")
      )
    }

    if (!is.null(cell_type)) {
      meta <- object@meta.data
      if (!cell_type_col %in% colnames(meta)) {
        stop("Seurat metadata must contain a '", cell_type_col, "' column.")
      }
      keep_cells <- rownames(meta)[meta[[cell_type_col]] == cell_type]
      expr <- expr[, colnames(expr) %in% keep_cells, drop = FALSE]
      message(
        "Subsetted to ", length(keep_cells), " cells where ", cell_type_col,
        " = '", cell_type, "'."
      )
    }
  } else if (inherits(object, "SingleCellExperiment")) {
    available_assays <- SummarizedExperiment::assayNames(object)
    assay_to_use <- if (!is.null(assay)) assay else "logcounts"

    if (!assay_to_use %in% available_assays) {
      stop(
        "Requested assay '", assay_to_use, "' not found. Available assays: ",
        paste(available_assays, collapse = ", ")
      )
    }
    expr <- SummarizedExperiment::assay(object, assay_to_use)
    message("Using SCE assay '", assay_to_use, "' (assumed log-normalized).")

    if (!is.null(cell_type)) {
      meta <- as.data.frame(SummarizedExperiment::colData(object))
      if (!cell_type_col %in% colnames(meta)) {
        stop("SCE colData must contain a '", cell_type_col, "' column.")
      }
      keep_cells <- rownames(meta)[meta[[cell_type_col]] == cell_type]
      expr <- expr[, colnames(expr) %in% keep_cells, drop = FALSE]
      message(
        "Subsetted to ", length(keep_cells), " cells where ", cell_type_col,
        " = '", cell_type, "'."
      )
    }
  } else if (is.matrix(object)) {
    if (!is.null(cell_type)) {
      stop("'cell_type' filtering is not supported when input is a raw matrix.")
    }
    expr <- object
  } else {
    stop("Input must be a Seurat object, SingleCellExperiment object, or numeric matrix.")
  }

  expr <- as.matrix(expr)
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]

  # Optional gene filtering
  gene_names <- rownames(expr)
  keep_genes <- rep(TRUE, length(gene_names))

  if (remove_mt) {
    keep_genes <- keep_genes & !grepl("^MT-", gene_names, ignore.case = TRUE)
    message("Removed mitochondrial genes matching '^MT-'.")
  }

  if (remove_rib) {
    keep_genes <- keep_genes & !grepl("^RP[SL]", gene_names, ignore.case = TRUE)
    message("Removed ribosomal genes matching '^RP[SL]'.")
  }

  expr <- expr[keep_genes, , drop = FALSE]

  avg_expression <- rowMeans(expr, na.rm = TRUE)
  sorted_genes <- names(sort(avg_expression, decreasing = TRUE))
  selected_genes <- head(sorted_genes, top_n)
  message("Top ", top_n, " genes selected based on mean expression.")

  return(selected_genes)
}
