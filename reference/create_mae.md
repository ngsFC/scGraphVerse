# Create MultiAssayExperiment from Multiple Single-Cell Datasets

Converts a list of count matrices, Seurat objects, or
SingleCellExperiment objects into a MultiAssayExperiment for integrated
network inference.

## Usage

``` r
create_mae(datasets, colData = NULL, ...)
```

## Arguments

- datasets:

  A named list of datasets. Each element can be:

  - A matrix (genes × cells)

  - A Seurat object

  - A SingleCellExperiment object

- colData:

  Optional. A DataFrame with metadata for each experiment. If NULL,
  automatically generated from list names.

- ...:

  Additional arguments (currently unused)

## Value

A MultiAssayExperiment object with:

- experiments: List of SingleCellExperiment objects

- colData: Metadata for each experiment/condition

## Examples

``` r
# Load the example MAE
data("toy_counts")

# Extract the list of SingleCellExperiment objects
sce_list <- MultiAssayExperiment::experiments(toy_counts)
sce_list <- as.list(sce_list)

# Create a new MAE from the SCE list
mae <- create_mae(sce_list)
mae
#> A MultiAssayExperiment object of 3 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 3:
#>  [1] experiment_1: SingleCellExperiment with 35 rows and 40 columns
#>  [2] experiment_2: SingleCellExperiment with 35 rows and 40 columns
#>  [3] experiment_3: SingleCellExperiment with 35 rows and 40 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```
