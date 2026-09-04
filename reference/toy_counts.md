# Toy MultiAssayExperiment for Network Inference

A simulated dataset generated using zinb_simdata() for demonstrating
network inference functions in the scGraphVerse package.

## Usage

``` r
data(toy_counts)
```

## Format

A MultiAssayExperiment object containing 3 SingleCellExperiment objects
(experiments), each with 10 genes x 40 cells. The data was simulated
using a known ground truth network (toy_adj_matrix) with zero-inflated
negative binomial distributions.

## Examples

``` r
data(toy_counts)
toy_counts
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

# Access individual experiments
MultiAssayExperiment::experiments(toy_counts)
#> ExperimentList class object of length 3:
#>  [1] experiment_1: SingleCellExperiment with 35 rows and 40 columns
#>  [2] experiment_2: SingleCellExperiment with 35 rows and 40 columns
#>  [3] experiment_3: SingleCellExperiment with 35 rows and 40 columns

# Use directly with infer_networks
networks <- infer_networks(
    count_matrices_list = toy_counts,
    method = "GENIE3",
    nCores = 1
)
```
