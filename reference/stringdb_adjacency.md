# Build Adjacency Matrices for Physical Interactions from STRING (POST API)

Constructs weighted and binary adj matrices for physical protein-protein
interactions using a POST request to the STRING database API.

## Usage

``` r
stringdb_adjacency(
  genes,
  species = 9606,
  required_score = 400,
  keep_all_genes = TRUE,
  verbose = TRUE
)
```

## Arguments

- genes:

  A character vector of gene symbols or identifiers, e.g.,
  `c("TP53", "BRCA1", ...)`.

- species:

  Integer. NCBI taxonomy ID of the species. Default is `9606` (human).

- required_score:

  Integer in \\0,1000\\. Minimum confidence score for interactions.
  Default is `400`.

- keep_all_genes:

  Logical. If `TRUE` (default), includes all input genes in the final
  matrix even if unmapped.

- verbose:

  Logical. If `TRUE`, displays progress messages. Default is `TRUE`.

## Value

A list containing:

- `weighted`: A square numeric adjacency matrix with scores as weights.

- `binary`: A corresponding binary (0/1) adjacency matrix.

## Details

This function:

1.  Maps input genes to STRING internal IDs.

2.  Uses a POST request to retrieve physical protein-protein
    interactions from STRING.

3.  Builds a weighted adjacency matrix using the STRING combined score.

4.  Builds a binary adjacency matrix indicating presence/absence.

Genes not mapped to STRING are optionally retained as zero rows/columns
if `keep_all_genes = TRUE`.

## Note

Requires packages: STRINGdb, httr, jsonlite.

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

str_res <- stringdb_adjacency(
    genes = genes,
    species = 9606,
    required_score = 900,
    keep_all_genes = FALSE
)
#> Initializing STRINGdb...
#> Registered S3 method overwritten by 'gplots':
#>   method         from     
#>   reorder.factor DescTools
#> Mapping genes to STRING IDs...
#> Mapped 5 genes to STRING IDs.
#> Retrieving **physical** interactions from STRING API...
#> Found 3 STRING physical interactions.
#> Adjacency matrices constructed successfully.
wadj_truth <- str_res$weighted
toy_adj_matrix <- str_res$binary
```
