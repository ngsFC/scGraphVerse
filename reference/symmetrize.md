# Symmetrize Adjacency Matrices in a SummarizedExperiment

Symmetrizes each adjacency matrix in a
[SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
by ensuring entries (i, j) and (j, i) are identical, using a specified
combination function.

## Usage

``` r
symmetrize(matrix_list, weight_function = "mean", nCores = 1)
```

## Arguments

- matrix_list:

  A
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object containing adjacency matrices to symmetrize.

- weight_function:

  Character string or function. Method to combine entries (i, j) and (j,
  i). Options include `"mean"`, `"max"`, `"min"`, or a user-defined
  function.

- nCores:

  Integer. Number of CPU cores to use for parallel processing. Defaults
  to the number of available workers in the current BiocParallel
  backend.

## Value

A
[SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object with symmetrized adjacency matrices, where for each matrix
\\A\[i, j\] = A\[j, i\]\\ for all \\i \neq j\\.

## Details

For each pair of off-diagonal elements (i, j) and (j, i):

- If one value is zero, the non-zero value is used.

- If both are non-zero, they are combined using the specified
  `weight_function`.

Diagonal entries are preserved as-is and not modified.

Parallelization is managed via BiocParallel for improved performance.

## Examples

``` r
data("toy_counts")


# Infer networks (toy_counts is already a MultiAssayExperiment)
networks <- infer_networks(
    count_matrices_list = toy_counts,
    method = "GENIE3",
    nCores = 1
)
head(networks[[1]])
#>   regulatoryGene targetGene    weight
#> 1          HLA-B        FTL 0.1943541
#> 2          HLA-A      HLA-B 0.1647698
#> 3           CD74      CXCR4 0.1615967
#> 4          HLA-B      HLA-A 0.1468035
#> 5            FTL       FTH1 0.1412889
#> 6           FTH1        FTL 0.1398331

# Generate adjacency matrices
wadj_se <- generate_adjacency(networks)
swadj_se <- symmetrize(wadj_se, weight_function = "mean")
```
