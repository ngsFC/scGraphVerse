# Modify Cell Names and Combine Datasets

Extracts expression data from a MultiAssayExperiment object, modifies
cell identifiers by appending a unique experiment index (e.g., "-m1",
"-m2", etc.), and merges the datasets into a single object.

## Usage

``` r
earlyj(input_list, rowg = TRUE)
```

## Arguments

- input_list:

  A MultiAssayExperiment object containing expression data from multiple
  experiments or conditions.

- rowg:

  Logical. If `TRUE` (default), genes are assumed to be rows and cells
  columns. If `FALSE`, matrices are transposed before renaming and
  combining.

## Value

A MultiAssayExperiment object containing a single merged experiment with
modified (unique) cell names.

## Details

For matrices, this function optionally transposes the input before
combining. For `Seurat` and `SingleCellExperiment` objects, only
features (genes) common across all input datasets are retained before
merging. The cell names are suffixed with "-m1", "-m2", etc., according
to their original list position. The result is returned as a
MultiAssayExperiment with a single merged experiment.

## Examples

``` r
data(toy_counts)

merged_mae <- earlyj(toy_counts)
merged_mae
#> A MultiAssayExperiment object of 1 listed
#>  experiment with a user-defined name and respective class.
#>  Containing an ExperimentList class object of length 1:
#>  [1] merged: SingleCellExperiment with 35 rows and 120 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```
