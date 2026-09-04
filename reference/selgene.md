# Select Top Expressed Genes from Single-Cell Data

Identifies and returns the top `n` most highly expressed genes across
all cells or within a specific cell type. Supports objects of class
Seurat,
[SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html),
or a numeric expression matrix (genes × cells).

## Usage

``` r
selgene(
  object,
  top_n,
  cell_type = NULL,
  cell_type_col = "cell_type",
  assay = NULL,
  remove_mt = FALSE,
  remove_rib = FALSE
)
```

## Arguments

- object:

  A Seurat object, SingleCellExperiment object, or numeric matrix (genes
  × cells).

- top_n:

  Integer. Number of top expressed genes to return.

- cell_type:

  Optional string. If provided, filters the expression matrix to only
  include cells of this type.

- cell_type_col:

  Character. Name of the column in metadata (Seurat `meta.data` or SCE
  `colData`) containing cell type annotations. Default is `"cell_type"`.

- assay:

  Character. For SingleCellExperiment objects only. Name of the assay to
  use. If `NULL`, defaults to `"logcounts"`.

- remove_mt:

  Logical. If `TRUE`, remove mitochondrial genes matching `"^MT-"`
  (case-insensitive).

- remove_rib:

  Logical. If `TRUE`, remove ribosomal genes matching `"^RP[SL]"`
  (case-insensitive).

## Value

A character vector of the top `n` most highly expressed gene names.

## Details

The function assumes that log-normalized values are available in the
`"data"` slot (for Seurat objects) or the `"logcounts"` assay (for
SingleCellExperiment). If raw counts are provided as a matrix, no
transformation is applied.

Optional filtering is available to exclude mitochondrial genes
(`"^MT-"`) and ribosomal genes (`"^RP[SL]"`), which may otherwise
dominate the top expressed genes.

## Details

When using a Seurat object, the function retrieves the log-normalized
data from the default assay's `"data"` slot. For SingleCellExperiment,
it uses the specified assay (default is `"logcounts"`). For matrices, no
checks or transformations are applied, and subsetting by cell type is
not supported.

Mitochondrial and ribosomal gene removal is based on regular expressions
matching gene names. These should follow standard naming conventions
(e.g., `MT-ND1`, `RPL13A`, `RPS6`).

## See also

[SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)

## Examples

``` r

data(toy_counts)
genes <- selgene(
    object = toy_counts[[1]],
    top_n = 5,
    cell_type = "T_cells",
    cell_type_col = "CELL_TYPE",
    remove_rib = TRUE,
    remove_mt = TRUE,
    assay = "counts"
)
#> Using SCE assay 'counts' (log-normalized).
#> Subsetted to 40 cells where CELL_TYPE = 'T_cells'.
#> Removed mitochondrial genes matching '^MT-'.
#> Removed ribosomal genes matching '^RP[SL]'.
#> Top 5 genes selected based on mean expression.
```
